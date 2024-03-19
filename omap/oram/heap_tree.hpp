#pragma once
#include "common/utils.hpp"
#include "external_memory/cachefrontvector.hpp"

/**
 * @brief A tree stored in a heap structure, but pack several levels together,
 * and optimize for reverse lexico order eviction.
 * Below shows an exmaple of a tree with 6 leaves and 4 levels, the number in
 * the brackets are the indices of the nodes in the data array.
 *
 *  Cached:                 ***[0]
 *                **0[1]               **1[2]
 * ---------------------------------------------------------
 *      *00[3]    |   *01[6]      |    *10[9]    |      *11[10]
 * 000[4]  100[5] | 001[7] 101[8] |              |
 *
 * @tparam T The type of node in the tree
 * @tparam PositionType
 * @tparam page_size The size of the page, default to 4096 bytes
 * @tparam evict_freq Number of reverse lexico evictions per random access
 */
template <typename T, typename PositionType = uint64_t,
          const size_t page_size = 4096, const int evict_freq = 2>
struct HeapTree {
 private:
  static const int evict_freq_ = evict_freq;
  static constexpr size_t max_node_per_page = divRoundUp(page_size, sizeof(T));
  static constexpr int node_per_page_log2 =
      GetLogBaseTwoConstexpr(max_node_per_page) + 1;

  /** By reducing the size of each packed subtree and store the packed tree in
   * each level following the reverse lexico order. We can trade off locality on
   * each path to locality in eviction*/
  static constexpr int findBestLexicoGroupLevel() {
    double minMissRate = 1.0 + double(evict_freq);
    int bestLevel = 0;
    for (int i = 0; i < node_per_page_log2; ++i) {
      double missRate =
          (1.0 + double(evict_freq) / (1 << i)) / (node_per_page_log2 - i);
      if (missRate < minMissRate) {
        minMissRate = missRate;
        bestLevel = i;
      }
    }
    return bestLevel;
  }
  static constexpr int bestLexicoGroupLevel = findBestLexicoGroupLevel();
  static constexpr int packLevel = node_per_page_log2 - bestLexicoGroupLevel;
  static constexpr PositionType packed_size = (1UL << packLevel) - 1;
  // put several packs to the same page
  static constexpr size_t node_per_page = packed_size << bestLexicoGroupLevel;
  static constexpr size_t actual_page_size = node_per_page * sizeof(T);
  using Vec = EM::CacheFrontVector::Vector<T, actual_page_size, true, true>;
  Vec arr;  // the actual data
  int cacheLevel = 0;
  int totalLevel = 0;
  PositionType leafCount = 0;
  PositionType cacheSize = 0;
  PositionType extSize = 0;
  PositionType totalSize = 0;  // total number of nodes

 public:
  HeapTree() {}

  HeapTree(PositionType _size, int _cacheLevel = 62) {
    Init(_size, _cacheLevel);
  }

  int GetDepth() const { return totalLevel; }

  /**
   * @brief Initialize the tree with a default value
   *
   * @param _size Number of leaves in the tree
   * @param defaultVal
   * @param _cacheLevel The number of top levels to cache, which affects the
   * layout of the tree.
   */
  void InitWithDefault(PositionType _size, const T& defaultVal,
                       int _cacheLevel = 62) {
    if (totalSize != 0) {
      throw std::runtime_error("Init called on non-empty tree");
    }
    cacheLevel = _cacheLevel;
    leafCount = _size;
    totalLevel = (int)GetLogBaseTwo(_size - 1) + 2;
    if (totalLevel > 64) {
      throw std::runtime_error("ORAM size too large.");
    }
    totalSize = 2 * _size - 1;
    cacheSize =
        (PositionType)std::min((uint64_t)totalSize, (2UL << _cacheLevel) - 1);
    arr.SetSize(totalSize, cacheSize, defaultVal);
    extSize = totalSize - cacheSize;
  }

  void Init(PositionType _size, int _cacheLevel = 62) {
    InitWithDefault(_size, T(), _cacheLevel);
  }

  /**
   * @brief Get the memory usage of the tree
   *
   * @param _size Number of leaves in the tree
   * @param _cacheLevel The number of top levels to cache.
   * @return uint64_t The memory usage in bytes
   */
  static uint64_t GetMemoryUsage(PositionType _size, int _cacheLevel = 62) {
    size_t totalSize = 2 * _size - 1;
    size_t cacheSize = std::min((uint64_t)totalSize, (2UL << _cacheLevel) - 1);
    return Vec::GetMemoryUsage(totalSize, cacheSize);
  }

  /**
   * @brief Get the memory usage of the tree
   *
   * @return uint64_t The memory usage in bytes of this tree
   */
  uint64_t GetMemoryUsage() const { return arr.GetMemoryUsage(); }

  /**
   * @brief Retrieve the node indices on a path idx
   *
   * @tparam Iterator
   * @param outputBegin
   * @param outputEnd
   * @param idx idx of the path
   * @param leafCount the number of leaves in the tree
   * @param cacheLevel the number of top levels to cache
   * @return int the number of nodes in the path
   */
  template <typename Iterator>
  static int GetPathIdx(Iterator outputBegin, Iterator outputEnd,
                        PositionType idx, PositionType leafCount,
                        int cacheLevel) {
    int totalLevel = (int)GetLogBaseTwo(leafCount - 1) + 2;
    // when all the levels are cached
    if ((packLevel == 1 && cacheLevel < totalLevel) ||
        cacheLevel == totalLevel - 1) {
      return GetPathIdx(outputBegin, outputEnd, idx, leafCount, totalLevel);
    }
    // the number of nodes in the path
    int pathLen = (idx | (PositionType)(1UL << (totalLevel - 2))) < leafCount
                      ? totalLevel
                      : totalLevel - 1;
    auto it = outputBegin;
    int i;
    // set top levels
    for (i = 0; i < std::min(pathLen, cacheLevel); ++i, ++it) {
      PositionType prevNodes = (PositionType)((1UL << i) - 1);
      PositionType subTreeIdx = idx & prevNodes;
      if ((1UL << i) <= leafCount) {
        *it = prevNodes + subTreeIdx;
      } else {
        PositionType mask = (PositionType)(1UL << (totalLevel - 2));
        // number of deleted nodes on the last level left to path
        PositionType deletedNode = (PositionType)std::max(
            (int64_t)(std::min(idx, mask) + mask - leafCount), 0L);
        *it = prevNodes + subTreeIdx - deletedNode;
      }
    }
    // return if all the nodes are cached in top levels
    if (i == pathLen) {
      return pathLen;
    }

    // set fully packed levels in the middle
    for (; i < totalLevel - packLevel - 1; i += packLevel) {
      PositionType prevNodes = (PositionType)((1UL << i) - 1);
      // begin offset of the pack
      PositionType beginOffset = prevNodes + (idx & prevNodes) * packed_size;
      // within each pack, nodes follow a standard heap layout
      for (int j = 0; j < packLevel; ++j, ++it) {
        PositionType innerPrevNodes = (PositionType)((1UL << j) - 1);
        *it = beginOffset + innerPrevNodes + ((idx >> i) & innerPrevNodes);
      }
    }

    // set levels with remaining nodes, which forms a packed tree
    PositionType prevNodes = (PositionType)((1UL << i) - 1);
    int remainLevel = totalLevel - i;
    PositionType subTreeIdx = idx & prevNodes;

    // Number of leaves to the left of the path, which equals
    // 2 * (number of (totalLevel - 1)-bit elements less than size, whose last i
    // bits < subtree idx)
    PositionType prevLeafCount = subTreeIdx * (leafCount >> i) +
                                 std::min(subTreeIdx, leafCount & prevNodes);
    // begin offset of the last pack
    PositionType beginOffset = prevNodes + prevLeafCount * 2 - subTreeIdx;
    // number of (totalLevel - 1)-bit elements < size, whose last i bits =
    // subtree idx
    PositionType subtreeLeafCount = ((leafCount - subTreeIdx - 1) >> i) + 1;
    GetPathIdx(it, outputEnd, idx >> i, subtreeLeafCount, remainLevel);
    for (; i < pathLen; ++i, ++it) {
      *it += beginOffset;
    }

    return pathLen;
  }

  /**
   * @brief Retrieve the node on the path
   *
   * @tparam Iterator The type of the iterator to store the path
   * @param pos The index of the path
   * @param pathBegin The begin iterator to store the path
   * @return int The number of nodes in the path
   */
  template <typename Iterator>
  int ReadPath(PositionType pos, Iterator pathBegin) {
    PositionType pathIdx[64];

    int actualLevel = GetPathIdx(&pathIdx[0], &pathIdx[totalLevel], pos,
                                 leafCount, cacheLevel);

    for (int i = 0; i < actualLevel; ++i) {
      PositionType idx = pathIdx[i];
      *(pathBegin + i) = arr.Get(idx);
    }
    return actualLevel;
  }

  /**
   * @brief Write the path to the tree
   *
   * @tparam Iterator The type of the iterator to store the path
   * @param pos The index of the path
   * @param pathBegin The begin iterator to store the path
   * @return int The number of nodes in the path
   */
  template <typename Iterator>
  int WritePath(PositionType pos, const Iterator pathBegin) {
    PositionType pathIdx[64];
    int actualLevel = GetPathIdx(&pathIdx[0], &pathIdx[totalLevel], pos,
                                 leafCount, cacheLevel);
    // we have checked that total level <= 64
    Assert(actualLevel <= totalLevel);
    for (int i = 0; i < actualLevel; ++i) {
      PositionType idx = pathIdx[i];
      arr[idx] = *(pathBegin + i);
    }
    return actualLevel;
  }
};