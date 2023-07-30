#pragma once
#include "common/defs.hpp"
#include <list>
#include <unordered_map>
#include <algorithm>

namespace CACHE {
template<typename T, uint64_t CACHE_SIZE=SERVER__CACHE_SIZE, uint64_t TLB_SIZE = 2>
struct Cache {
  struct CacheEntry {
    T val;
    bool dirty;
  };
  typedef uint64_t IndexType;
  uint64_t size;
  std::unordered_map<IndexType, CacheEntry> data;
  std::unordered_map<IndexType, std::list<IndexType>::iterator> dataRefs;
  std::list<IndexType> LRU;

  struct TLB_entry {
    IndexType idx;
    uint64_t counter = 0; // 0 stands for invalid
    CacheEntry* pentry; // rehashing should not invalidate pointer to value
  };

  uint64_t TLB_global_counter;
  std::vector<TLB_entry> TLB;

  void sortTLBByCounter() {
    auto cmp = [](const TLB_entry& entry1, const TLB_entry& entry2) {
      return entry1.counter > entry2.counter;
    };
    std::sort(TLB.begin(), TLB.end(), cmp);
  }

  void resetTLBCounters() {
    sortTLBByCounter();
    TLB_global_counter = TLB_SIZE + 1;
    uint64_t counter = TLB_SIZE + 1;
    for (TLB_entry& entry: TLB) {
      entry.counter = --counter;
    }
  }

  // if write back is false, reset dirty to false
  T* AccessTLB(const IndexType& accessedIndex, bool dirty = true, bool writeBack = true) {
    if (++TLB_global_counter == UINT64_MAX) {
      // to prevent overflow
      resetTLBCounters();
    }
    for (TLB_entry& entry: TLB) {
      if (entry.idx == accessedIndex && entry.counter) {
        entry.counter = TLB_global_counter;
        entry.pentry->dirty = (entry.pentry->dirty || dirty) && writeBack;
        T* pval = &entry.pentry->val;
        return pval;
      }
    }
    return NULL;
  }

  const T* AccessTLBReadOnly(const IndexType& accessedIndex) {
    return AccessTLB(accessedIndex, false);
  }

  const T* AccessTLBNoWriteBack(const IndexType& accessedIndex) {
    return AccessTLB(accessedIndex, false, false);
  }

  T& UpdateTLB(const IndexType& accessedIndex, CacheEntry& entry) {
    uint64_t minCounter = UINT64_MAX;
    auto minIt = TLB.begin();
    for (auto it = TLB.begin(); it != TLB.end(); ++it) {
      if (it->counter < minCounter) {
        minCounter = it->counter;
        minIt = it;
      }
    }
    *minIt = {accessedIndex, TLB_global_counter, &entry};
    sortTLBByCounter(); // speed up future lookups
    return entry.val;
  }

  Cache() {
    size = 0;
    data.reserve(CACHE_SIZE);
    dataRefs.reserve(CACHE_SIZE);
    TLB_global_counter = 0;
    TLB.resize(TLB_SIZE);
  }

  bool CheckContains(const IndexType& rootIdx) {
    return data.count(rootIdx) > 0;
  }

  bool IsFull() {
    return size == CACHE_SIZE;
  }

  CacheEntry& GetNextToEvict(IndexType& indexToEvict) {
    Assert(size > 0);
    Assert(size == CACHE_SIZE);
    indexToEvict = LRU.back();
    Assert(data.count(indexToEvict) > 0);
    return data[indexToEvict];
  }

  void EvictLRU(const IndexType& evictedIndex) {
    #ifndef NDEBUG
    Assert(size > 0);
    Assert(size == CACHE_SIZE);
    Assert(evictedIndex == LRU.back());
    Assert(data.count(evictedIndex) > 0);
    #endif
    size--;
    LRU.pop_back();
    data.erase(evictedIndex);
    dataRefs.erase(evictedIndex);
  }

  void Insert(const IndexType& newIndex, const T& val) {
    Assert(size < CACHE_SIZE);
    Assert(data.count(newIndex) == 0);
    size++;
    LRU.push_front(newIndex);

    data[newIndex] = { val, false };
    dataRefs[newIndex] = (LRU.begin());

    // return data[newIndex].val;
  }

  T& Access(const IndexType& accessedIndex, bool dirty = true, bool writeBack = true) {
    Assert(size > 0);
    Assert(data.count(accessedIndex) > 0);
    
    LRU.erase(dataRefs[accessedIndex]);
    LRU.push_front(accessedIndex);
    dataRefs[accessedIndex] = LRU.begin();
    auto& data_to_access = data[accessedIndex];
    data_to_access.dirty = (data_to_access.dirty || dirty) && writeBack;

    return UpdateTLB(accessedIndex, data_to_access);
  }

  const T& AccessReadOnly(const IndexType& accessedIndex) {
    return Access(accessedIndex, false);
  }

  const T& AccessNoWriteBack(const IndexType& accessedIndex) {
    return Access(accessedIndex, false, false);
  }
};
}