#pragma once
#include <algorithm>
#include <list>
#include <unordered_map>

#include "common/defs.hpp"
#include "common/dmcache.hpp"
#include "common/tracing/perf.hpp"

#define USE_LRU_FLAG false

namespace CACHE {
enum CacheType {
  LRU,
  DM  // direct map
};
template <typename T, typename S, uint64_t CACHE_SIZE = SERVER__CACHE_SIZE,
          uint64_t TLB_SIZE = 2>
struct Cached : S {
  using Cache = typename std::conditional<USE_LRU_FLAG,
                                          CACHE::Cache<T, CACHE_SIZE, TLB_SIZE>,
                                          CACHE::DMCache<T, CACHE_SIZE> >::type;
  Cache cache;
  using typename S::IndexType;

  using S::S;

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

  template <bool dirty = true, bool skip_read = false, bool writeBack = true>
  T& Access(const IndexType i) {
    if constexpr (USE_LRU_FLAG) {
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