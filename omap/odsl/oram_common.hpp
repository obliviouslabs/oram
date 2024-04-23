#pragma once

#include <concepts>
#include <vector>

/// @brief This file defines some methods and concepts common to binary-tree
/// based ORAMs.

#include "algorithm/bitonic.hpp"
#include "algorithm/or_compact_shuffle.hpp"
#include "bucket.hpp"
#include "external_memory/rw_concepts.hpp"
#include "external_memory/stdvector.hpp"
#include "external_memory/virtualvector.hpp"
#include "heap_tree.hpp"

namespace ODSL {

/**
 * @brief An update or remove function should be able to update an object of
 * type T and return a boolean value indicating whether the object should be
 * kept in the ORAM.
 *
 * @tparam Func The type of the update function.
 * @tparam T The type of the object to be updated.
 */
template <class Func, typename T>
concept UpdateOrRemoveFunction = requires(Func f, T& t) {
  { f(t) } -> std::same_as<bool>;
};

/**
 * @brief An update function should be able to update an object of
 * type T.
 *
 * @tparam Func The type of the update function.
 * @tparam T The type of the object to be updated.
 */
template <class Func, typename T>
concept UpdateFunction = requires(Func f, T& t) {
  { f(t) };
};

/**
 * @brief A batch update or remove function should be able to update a batch of
 * object of type T and return a vector of boolean value indicating whether each
 * object should be kept in the ORAM.
 *
 * @tparam Func The type of the update function.
 * @tparam T The type of the object to be updated.
 */
template <class Func, typename T>
concept BatchUpdateOrRemoveFunction =
    requires(Func f, T* t, uint64_t batchSize) {
      { f(batchSize, t) } -> std::same_as<std::vector<bool>>;
    };

/**
 * @brief A batch update or remove function should be able to update a batch of
 * object of type T.
 *
 * @tparam Func The type of the update function.
 * @tparam T The type of the object to be updated.
 */
template <class Func, typename T>
concept BatchUpdateFunction = requires(Func f, std::vector<T>& t) {
  { f(t) };
};

/**
 * @brief Return the length of the common suffix of two positions.
 *
 * @tparam PositionType
 * @param a First position
 * @param b Second position
 * @return int The length of the common suffix
 */
template <typename PositionType>
int CommonSuffixLength(PositionType a, PositionType b) {
  return std::countr_zero(a ^ b);
}

/**
 * @brief Read an element from the path and remove it from the path.
 *
 * @tparam Iterator
 * @tparam UidType
 * @tparam T
 * @param begin The begin iterator of the path
 * @param end The end iterator of the path
 * @param uid The uid of the element to be read
 * @param out The output of the element
 * @return true if the element is read successfully, false otherwise
 */
template <class Iterator, typename UidType, typename T>
bool ReadElementAndRemoveFromPath(Iterator begin, Iterator end,
                                  const UidType& uid, T& out) {
  bool findFlag = false;
  for (auto it = begin; it != end; ++it) {
    findFlag |= it->ReadAndRemove(uid, out);
  }
  return findFlag;
}

/**
 * @brief Write a new block to the path.
 *
 * @tparam Iterator
 * @tparam Block_
 * @param begin The begin iterator of the path
 * @param end The end iterator of the path
 * @param block The block to be written
 * @return true if the block is written successfully, false otherwise
 */
template <class Iterator, typename Block_>
bool WriteNewBlockToPath(Iterator begin, Iterator end, const Block_& block) {
  bool cond = block.IsDummy();
  // fill the first slot that's empty
  for (auto it = begin; it != end; ++it) {
    cond |= it->Insert(!cond, block);
  }
  return cond;
}

/**
 * @brief Calculate the maximum number of levels to cache in the heap tree.
 *
 * @tparam T The type of the data stored in the ORAM
 * @tparam Z The bucket size of the ORAM
 * @tparam stashSize The size of the ORAM stash
 * @tparam PositionType The type of the position
 * @tparam UidType The type of the unique id
 * @tparam check_freshness Whether to check freshness of the data
 */
template <typename T, const int Z, const int stashSize,
          typename PositionType = uint64_t, typename UidType = uint64_t,
          const bool check_freshness = true>
int GetMaxCacheLevel(PositionType size, size_t cacheBytes = MAX_CACHE_SIZE) {
  using Bucket_ = Bucket<T, Z, PositionType, UidType>;
  using FreshORAMNode_ = FreshORAMNode<Bucket_>;

  using TreeNode_ =
      std::conditional_t<check_freshness, FreshORAMNode_, Bucket_>;
  using Stash = Bucket<T, stashSize, PositionType, UidType>;
  using HeapTree_ = HeapTree<TreeNode_, PositionType>;
  int maxCacheLevel = 0;
  if (cacheBytes >= HeapTree_::GetMemoryUsage(size, MAX_CACHE_LEVEL)) {
    // default value
    return MAX_CACHE_LEVEL;
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
