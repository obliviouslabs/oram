#pragma once
#include "heap_tree.hpp"

/**
 * @brief Separate the tree into multiple subtrees, each with its own cache to
 * avoid contention
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
struct ParHeapTree {
  HeapTree<T, PositionType, page_size, evict_freq> topTree;
  std::vector<HeapTree<T, PositionType, page_size, evict_freq>> subTrees;
  int cacheLevel = 0;
  int topLevel = 0;
  int totalLevel = 0;
  PositionType numSubTree = 0;
  PositionType subtreeMask = 0;
  PositionType topMask = 0;
  PositionType leafCount = 0;

  ParHeapTree() {}

  ParHeapTree(PositionType _size, int _topLevel = 1, int _cacheLevel = 62,
              PositionType realCacheSize = -1) {
    Init(_size, _topLevel, _cacheLevel, realCacheSize);
  }

  /**
   * @brief Initialize the tree with a default value
   *
   * @param _size Number of leaves in the tree
   * @param defaultVal
   * @param _topLevel The number of top levels shared by all subtrees
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
      throw std::runtime_error("cacheLevel must >= topLevel");
    }
    topLevel = _topLevel;
    numSubTree = 1UL << _topLevel;
    subtreeMask = numSubTree - 1;
    topMask = numSubTree / 2 - 1;
    _size = divRoundUp(_size, numSubTree) << _topLevel;
    leafCount = _size;
    cacheLevel = _cacheLevel;

    totalLevel = GetLogBaseTwo(_size - 1) + 2;
    topTree.InitWithDefault(numSubTree / 2, defaultVal, _cacheLevel,
                            realCacheSize);
    subTrees.resize(numSubTree);
    PositionType topTreeSize = (1UL << topLevel) - 1;
    PositionType subtreeRealCacheSize =
        (realCacheSize - topTreeSize) / numSubTree;
    if (realCacheSize == -1) {
      subtreeRealCacheSize = -1;
    } else if (realCacheSize < topTreeSize) {
      subtreeRealCacheSize = 0;
    }
    for (PositionType i = 0; i < numSubTree; ++i) {
      subTrees[i].InitWithDefault(_size >> _topLevel, defaultVal,
                                  remainCacheLevel, subtreeRealCacheSize);
    }
    printf("ParHeapTree: topLevel %d, cacheLevel %d, realCacheSize %lu\n",
           topLevel, cacheLevel, realCacheSize);
  }

  void Init(PositionType _size, int _cacheLevel = 62,
            PositionType realCacheSize = -1) {
    InitWithDefault(_size, T(), _cacheLevel, realCacheSize);
  }

  int GetCacheLevel() const { return cacheLevel; }

  size_t GetLeafCount() const { return leafCount; }

  static uint64_t GetMemoryUsage(PositionType _size, int _topLevel,
                                 int _cacheLevel = 62,
                                 PositionType realCacheSize = -1) {
    PositionType numSubTree = 1UL << _topLevel;
    _size = divRoundUp(_size, numSubTree) << _topLevel;
    PositionType topTreeSize = numSubTree - 1;
    PositionType subtreeRealCacheSize =
        (realCacheSize - topTreeSize) / numSubTree;
    if (realCacheSize == -1) {
      subtreeRealCacheSize = -1;
    } else if (realCacheSize < topTreeSize) {
      subtreeRealCacheSize = 0;
    }

    uint64_t subTreeMemoryUsage =
        HeapTree<T, PositionType, page_size, evict_freq>::GetMemoryUsage(
            _size >> _topLevel, _cacheLevel - _topLevel, subtreeRealCacheSize);
    uint64_t topTreeMemoryUsage =
        HeapTree<T, PositionType, page_size, evict_freq>::GetMemoryUsage(
            numSubTree, _cacheLevel, realCacheSize);
    return subTreeMemoryUsage * numSubTree + topTreeMemoryUsage;
  }

  uint64_t GetMemoryUsage() const {
    return GetMemoryUsage(leafCount, topLevel, cacheLevel);
  }

  template <typename Iterator>
  int ReadPath(PositionType pos, Iterator pathBegin) {
    int topPathLevel = topTree.ReadPath(pos & topMask, pathBegin);
    Assert(topPathLevel == topLevel);
    return topLevel + subTrees[pos & subtreeMask].ReadPath(
                          pos >> topLevel, pathBegin + topLevel);
  }

  template <typename Iterator>
  int ReadSubPath(PositionType pos, Iterator pathBegin, int k) {
    if (k == topLevel) {
      return topLevel +
             subTrees[pos & subtreeMask].ReadPath(pos >> topLevel, pathBegin);
    } else if (k < topLevel) {
      topTree.ReadSubPath(pos & topMask, pathBegin, k);
      return topLevel + subTrees[pos & subtreeMask].ReadPath(
                            pos >> topLevel, pathBegin + (topLevel - k));
    } else {
      return topLevel + subTrees[pos & subtreeMask].ReadSubPath(
                            pos >> topLevel, pathBegin, k - topLevel);
    }
  }

  template <typename Iterator>
  int WritePath(PositionType pos, const Iterator pathBegin) {
    topTree.WritePath(pos & topMask, pathBegin);
    return topLevel + subTrees[pos & subtreeMask].WritePath(
                          pos >> topLevel, pathBegin + topLevel);
  }

  template <typename Iterator>
  int WriteSubPath(PositionType pos, const Iterator pathBegin, int k) {
    if (k == topLevel) {
      return topLevel +
             subTrees[pos & subtreeMask].WritePath(pos >> topLevel, pathBegin);
    } else if (k < topLevel) {
      topTree.WriteSubPath(pos & topMask, pathBegin, k);
      return topLevel + subTrees[pos & subtreeMask].WritePath(
                            pos >> topLevel, pathBegin + (topLevel - k));
    } else {
      return topLevel + subTrees[pos & subtreeMask].WriteSubPath(
                            pos >> topLevel, pathBegin, k - topLevel);
    }
  }

  void getTopKLevel(int k, T* output) {
    if (k <= topLevel) {
      topTree.getTopKLevel(k, output);
    } else {
      topTree.getTopKLevel(k, output);
      PositionType topSize = (1UL << topLevel) - 1;
      PositionType subSize = (1UL << (k - topLevel)) - 1;
      for (PositionType i = 0; i < numSubTree; ++i) {
        subTrees[i].getTopKLevel(k - topLevel, output + topSize + i * subSize);
      }
    }
  }

  void SetByPathAndLevel(PositionType path, int level, const T& val) {
    if (level < topLevel) {
      topTree.SetByPathAndLevel(path, level, val);
    } else {
      subTrees[path & subtreeMask].SetByPathAndLevel(path >> topLevel,
                                                     level - topLevel, val);
    }
  }
};
