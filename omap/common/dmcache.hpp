#pragma once
#include <algorithm>
#include <list>
#include <unordered_map>

#include "common/defs.hpp"

namespace CACHE {
template <typename T, uint64_t cache_size = SERVER__CACHE_SIZE>
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
    size = cache_size;
    data.resize(cache_size);
  }

  bool CheckContains(const IndexType& rootIdx) {
    return data[rootIdx % cache_size].idx == rootIdx;
  }

  void Insert(const IndexType& newIndex, const T& val) {
    data[newIndex % cache_size] = {val, newIndex, false};
  }

  T& Access(const IndexType& accessedIndex, bool dirty = true,
            bool writeBack = true) {
    auto& data_to_access = data[accessedIndex % cache_size];
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
    return data[indexToEvict % cache_size];
  }

  const CacheEntry& GetMappedSlot(IndexType indexToEvict) const {
    return data[indexToEvict % cache_size];
  }
};
}  // namespace CACHE