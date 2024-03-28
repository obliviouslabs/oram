#pragma once
#include <algorithm>
#include <list>
#include <unordered_map>

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
                                CACHE::DMCache<T, cache_size> >::type;
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
};
}  // namespace CACHE