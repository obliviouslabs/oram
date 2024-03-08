#pragma once
#include "heap_tree.hpp"

/**
 * @brief A tree stored in a heap structure, but pack several levels together,
 * and optimize for reverse lexico order eviction.
 * Below shows an exmpale of a tree with 6 leaves and 4 levels, the number in
 * the brackets are the indices of the nodes in the data array.
 *
 *  Cached:                 ***[0]
 *                **0[1]               **1[2]
 * ---------------------------------------------------------
 *      *00[3]    |   *01[6]      |    *10[9]    |      *11[10]
 * 000[4]  100[5] | 001[7] 101[8] |              |
 *
 * @tparam T
 * @tparam PositionType
 * @tparam page_size The size of the page, default to 4096 bytes
 * @tparam evict_freq Number of reverse lexico evictions per random access
 */
template <typename T, typename PositionType = uint64_t,
          const size_t page_size = 4096, const int evict_freq = 2>
struct HeapTreeSub {
  HeapTree<T, PositionType, page_size, evict_freq> topTree;
  std::vector<HeapTree<T, PositionType, page_size, evict_freq>> subTrees;
  int cacheLevel = 0;
  int topLevel = 0;
  int totalLevel = 0;
  PositionType numSubTree = 0;
  PositionType leafCount = 0;

  HeapTreeSub() {}

  HeapTreeSub(PositionType _size, int _topLevel = 1, int _cacheLevel = 62,
              PositionType realCacheSize = -1) {
    Init(_size, _topLevel, _cacheLevel, realCacheSize);
  }

  /**
   * @brief Initialize the tree with a default value
   *
   * @param _size Number of leaves in the tree
   * @param defaultVal
   * @param _cacheLevel The number of top levels to cache, which determines the
   * logical structure of the tree
   * @param realCacheSize Allow caching more/fewer items physically.
   */
  void InitWithDefault(PositionType _size, int _topLevel, const T& defaultVal,
                       int _cacheLevel = 62, PositionType realCacheSize = -1) {
    if (leafCount != 0) {
      throw std::runtime_error("Init called on non-empty tree");
    }

    int remainCacheLevel = _cacheLevel - _topLevel;
    if (remainCacheLevel < 0) {
      throw std::runtime_error("cacheLevel must be >= topLevel");
    }
    topLevel = _topLevel;
    numSubTree = 1UL << _topLevel;
    _size = divRoundUp(_size, numSubTree) << _topLevel;
    leafCount = _size;
    cacheLevel = _cacheLevel;

    totalLevel = GetLogBaseTwo(_size - 1) + 2;
    topTree.InitWithDefault(numSubTree, defaultVal, _cacheLevel, realCacheSize);
    subTrees.resize(numSubTree);
    size_t subtreeRealCacheSize =
        realCacheSize == -1
            ? -1
            : (realCacheSize - (1UL << topLevel) + 1) / numSubTree;
    for (PositionType i = 0; i < numSubTree; ++i) {
      subTrees[i].InitWithDefault(_size >> _topLevel, defaultVal,
                                  remainCacheLevel, subtreeRealCacheSize);
    }
  }

  void Init(PositionType _size, int _cacheLevel = 62,
            PositionType realCacheSize = -1) {
    InitWithDefault(_size, T(), _cacheLevel, realCacheSize);
  }

  int GetCacheLevel() const { return cacheLevel; }

  size_t GetLeafCount() const { return leafCount; }
};
