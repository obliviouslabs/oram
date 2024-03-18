#pragma once

#include <functional>
#include <vector>

#include "bucket.hpp"
#include "external_memory/algorithm/bitonic.hpp"
#include "external_memory/algorithm/or_compact_shuffle.hpp"
#include "external_memory/stdvector.hpp"
#include "external_memory/virtualvector.hpp"
#include "heap_tree.hpp"

namespace ODSL {

template <typename PositionType>
int commonSuffixLength(PositionType a, PositionType b) {
  return std::countr_zero(a ^ b);
}

template <class PathVec, typename UidType, typename T>
void ReadElementAndRemoveFromPath(PathVec& path, const UidType& uid, T& out) {
  using Block_ = PathVec::value_type;
  for (Block_& b : path) {
    b.ReadAndRemove(uid, out);
  }
}

template <class Iterator, typename UidType, typename T>
void ReadElementAndRemoveFromPath(Iterator begin, Iterator end,
                                  const UidType& uid, T& out) {
  for (auto it = begin; it != end; ++it) {
    it->ReadAndRemove(uid, out);
  }
}

// Write a new block to the path, return false if the path is full
template <class PathVec, typename Block_>
bool WriteNewBlockToPath(PathVec& path, const Block_& block) {
  int endIdx = path.size();
  bool cond = true;
  // fill the first slot that's empty
  for (int i = 0; i < endIdx; i++) {
    cond &= !path[i].Insert(cond, block);
  }
  return !cond;
}

// Write a new block to the top of the tree, return false if the top is full
template <class PathVec, typename Block_>
bool WriteNewBlockToTreeTop(PathVec& path, const Block_& block, int topSize) {
  bool cond = true;
  // fill the first slot that's empty
  for (int i = 0; i < topSize; i++) {
    cond &= !path[i].Insert(cond, block);
  }
  return !cond;
}

template <typename T, const int Z, const int stashSize,
          typename PositionType = uint64_t, typename UidType = uint64_t>
int GetMaxCacheLevel(PositionType size, size_t cacheBytes = 1UL << 62) {
  using Bucket_ = Bucket<T, Z, PositionType, UidType>;
  using Stash = Bucket<T, stashSize, PositionType, UidType>;
  using HeapTree_ = HeapTree<Bucket_, PositionType>;
  int maxCacheLevel = 0;
  if (cacheBytes >= HeapTree_::GetMemoryUsage(size, 62)) {
    // default value
    return 62;
  }
  if (cacheBytes < sizeof(Stash)) {
    return -1;
  }
  cacheBytes -= sizeof(Stash);
  for (;; ++maxCacheLevel) {
    size_t treeUsage = HeapTree_::GetMemoryUsage(size, maxCacheLevel);
    if (cacheBytes < treeUsage) {
      break;
    }
  }
  return maxCacheLevel - 1;
}

}  // namespace ODSL
