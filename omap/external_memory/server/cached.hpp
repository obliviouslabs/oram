#pragma once
#include <algorithm>
#include <list>
#include <unordered_map>
#include <vector>

#include "common/defs.hpp"
#include "common/dmcache.hpp"
#include "common/lrucache.hpp"
#include "common/tracing/perf.hpp"

namespace CACHE {
enum CacheType {
  LRU,
  DM  // direct map
};
template <typename T, typename S, uint64_t cache_size = SERVER__CACHE_SIZE,
          uint64_t tlb_size = 2, CacheType cache_type = CacheType::DM>
struct Cached : S {
 private:
  using S::S;
  using typename S::IndexType;

  static constexpr bool use_lru_flag = cache_type == LRU;
  using Cache =
      typename std::conditional<use_lru_flag,
                                CACHE::LRUCache<T, cache_size, tlb_size>,
                                CACHE::DMCache<T, cache_size>>::type;
  Cache cache;

  template <bool dirty = true, bool skip_read = false, bool writeBack = true>
  T& AccessLRU(const IndexType i) {
    PERFCTR_INCREMENT(accessCount);
    T* TLB_result = cache.AccessTLB(i, dirty, writeBack);
    if (TLB_result) {
      PERFCTR_INCREMENT(tlbHitCount);
      return *TLB_result;
    }

    if (!cache.CheckContains(i)) {
      if (cache.IsFull()) {
        IndexType evictedIndex;
        {
          auto& evicted = cache.GetNextToEvict(evictedIndex);
          if (evicted.dirty) {
            Write(evictedIndex, evicted.val);
          }
        }
        cache.EvictLRU(evictedIndex);
      }

      T ret;
      if (!skip_read) {
        Read(i, ret);
      }
      cache.Insert(i, ret);
    }

    return cache.Access(i, dirty, writeBack);
  }

  template <bool dirty = true, bool skip_read = false, bool writeBack = true>
  T& AccessDM(const IndexType i) {
    PERFCTR_INCREMENT(accessCount);
    auto& mappedSlot = cache.GetMappedSlot(i);
    if (mappedSlot.idx != i) {
      if (mappedSlot.dirty) {
        Write(mappedSlot.idx, mappedSlot.val);
      }
      if constexpr (!skip_read) {
        Read(i, mappedSlot.val);
        mappedSlot.dirty = false;
      }
      mappedSlot.idx = i;
    }
    mappedSlot.dirty = (mappedSlot.dirty || dirty) && writeBack;
    return mappedSlot.val;
  }

  void ReadLazy(const IndexType idx, T& out) {
    PERFCTR_INCREMENT(accessCount);
    S::ReadLazy(idx, out);
  }

  void WriteLazy(const IndexType idx, const T& in) {
    PERFCTR_INCREMENT(accessCount);
    S::WriteLazy(idx, in);
  }

  void flushRead() { S::flushRead(); }

  void flushWrite() { S::flushWrite(); }

  void ReadThruDMLazy(const IndexType idx, T& out) {
    const auto& mappedSlot = cache.GetMappedSlot(idx);
    if (mappedSlot.idx != idx) {
      ReadLazy(idx, out);
    } else {
      out = mappedSlot.val;
      // std::memcpy(&out[i], &mappedSlot.val, sizeof(T));
    }
  }

  void WriteBackDMLazy(const IndexType idx, const T& in) {
    auto& mappedSlot = cache.GetMappedSlot(idx);
    if (mappedSlot.idx != idx) {
      if (mappedSlot.dirty) {
        WriteLazy(mappedSlot.idx, mappedSlot.val);
      }
    }
    mappedSlot.val = in;
    mappedSlot.idx = idx;
    mappedSlot.dirty = true;
  }

 public:
  static uint64_t GetMemoryUsage() {
    return cache_size * sizeof(typename Cache::CacheEntry);
  }

  template <bool dirty = true, bool skip_read = false, bool writeBack = true>
  T& Access(const IndexType i) {
    if constexpr (use_lru_flag) {
      return AccessLRU<dirty, skip_read, writeBack>(i);
    } else {
      return AccessDM<dirty, skip_read, writeBack>(i);
    }
  }

  const T& AccessReadOnly(const IndexType i) { return Access<false>(i); }

  const T& AccessNoWriteBack(const IndexType i) {
    return Access<false, false, false>(i);
  }

  T& AccessWriteOnly(const IndexType i) { return Access<true, true>(i); }

  void Write(const IndexType i, const T& in) { S::Write(i, in); }

  void Read(const IndexType i, T& out) { S::Read(i, out); }

  struct BatchAccessor {
   private:
    Cached& parent;
    std::vector<std::pair<IndexType, uint32_t>> indices;
    std::vector<uint32_t> prefetchCachePosMap;
    std::vector<T> prefetchCache;

    void readBatch() {
      if (indices.empty()) {
        return;
      }
      prefetchCachePosMap.resize(indices.size());
      std::sort(indices.begin(), indices.end(),
                [](const auto& a, const auto& b) { return a.first < b.first; });
      for (uint32_t i = 0; i < indices.size(); ++i) {
        IndexType idx = indices[i].first;
        if (i == 0 || idx != indices[i - 1].first) {
          prefetchCachePosMap[indices[i].second] = prefetchCache.size();
          if (prefetchCache.size() == prefetchCache.capacity()) {
            parent.flushRead();  // we have to read everything to the prefetch
                                 // cache before the reference is invalidated
          }
          prefetchCache.emplace_back();
          // TODO change to batch read
          parent.ReadThruDMLazy(idx, prefetchCache.back());
        } else {
          prefetchCachePosMap[indices[i].second] = prefetchCache.size() - 1;
        }
      }
      parent.flushRead();
    }

    void writeBackBatch() {
      for (uint32_t i = 0, cacheIdx = 0; i < indices.size(); ++i) {
        IndexType idx = indices[i].first;
        if (i == 0 || idx != indices[i - 1].first) {
          parent.WriteBackDMLazy(idx, prefetchCache[cacheIdx++]);
        }
      }
      parent.flushWrite();
    }

   public:
    BatchAccessor(Cached& parent, uint32_t prefetchCapacity) : parent(parent) {
      indices.reserve(prefetchCapacity * 2);
      prefetchCache.reserve(prefetchCapacity);
    }

    uint32_t Prefetch(const IndexType idx) {
      indices.emplace_back(idx, indices.size());
      return indices.size() - 1;
    }

    void FlushRead() { readBatch(); }

    T& Access(uint32_t prefetchReceipt) {
      return prefetchCache[prefetchCachePosMap[prefetchReceipt]];
    }

    void FlushWrite() {
      writeBackBatch();
      indices.clear();
      prefetchCache.clear();
    }

    ~BatchAccessor() {
      if (!indices.empty()) {
        FlushWrite();
      }
    }
  };
};
}  // namespace CACHE