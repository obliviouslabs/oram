#pragma once
#include "common/defs.hpp"
#include <list>
#include <unordered_map>
#include <algorithm>

namespace CACHE {
template<typename T, uint64_t CACHE_SIZE=SERVER__CACHE_SIZE>
struct DMCache {
  typedef uint64_t IndexType;
  struct CacheEntry {
    T val;
    IndexType idx;
    bool dirty;
  };
  
  uint64_t size;
  std::vector<CacheEntry> data;


  DMCache() {
    size = CACHE_SIZE;
    data.resize(CACHE_SIZE);
  }

  bool CheckContains(const IndexType& rootIdx) {
    return data[rootIdx % CACHE_SIZE].idx == rootIdx;
  }

  void Insert(const IndexType& newIndex, const T& val) {
    data[newIndex % CACHE_SIZE] = { val, newIndex, false };
  }

  T& Access(const IndexType& accessedIndex, bool dirty = true, bool writeBack = true) {

    auto& data_to_access = data[accessedIndex % CACHE_SIZE];
    data_to_access.dirty = (data_to_access.dirty || dirty) && writeBack;

    return data_to_access.val;
  }

  const T& AccessReadOnly(const IndexType& accessedIndex) {
    return Access(accessedIndex, false);
  }

  const T& AccessNoWriteBack(const IndexType& accessedIndex) {
    return Access(accessedIndex, false, false);
  }

  CacheEntry& GetMappedSlot(IndexType indexToEvict) {
    return data[indexToEvict % CACHE_SIZE];
  }
};
}