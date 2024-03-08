#pragma once
#include "common/utils.hpp"
#include "external_memory/cachefrontvector.hpp"

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
struct HeapTree {
  static const int evict_freq_ = evict_freq;
  static constexpr size_t max_node_per_page = divRoundUp(page_size, sizeof(T));
  static constexpr size_t node_per_page_log2 =
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

  HeapTree() {}

  HeapTree(PositionType _size, int _cacheLevel = 62,
           PositionType realCacheSize = -1) {
    Init(_size, _cacheLevel, realCacheSize);
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
  void InitWithDefault(PositionType _size, const T& defaultVal,
                       int _cacheLevel = 62, PositionType realCacheSize = -1) {
    if (totalSize != 0) {
      throw std::runtime_error("Init called on non-empty tree");
    }
    cacheLevel = _cacheLevel;
    leafCount = _size;
    totalLevel = GetLogBaseTwo(_size - 1) + 2;
    Assert(totalLevel <= 64);
    totalSize = 2 * _size - 1;
    cacheSize = std::min((uint64_t)totalSize, (2UL << _cacheLevel) - 1);
    arr.SetSize(totalSize, realCacheSize != -1 ? realCacheSize : cacheSize,
                defaultVal);
    extSize = totalSize - cacheSize;
  }

  void Init(PositionType _size, int _cacheLevel = 62,
            PositionType realCacheSize = -1) {
    InitWithDefault(_size, T(), _cacheLevel, realCacheSize);
  }

  int GetCacheLevel() const { return cacheLevel; }

  size_t GetLeafCount() const { return leafCount; }

  size_t GetNodeCount() const { return totalSize; }

  static uint64_t GetMemoryUsage(PositionType _size, int _cacheLevel = 62,
                                 PositionType realCacheSize = -1) {
    size_t totalSize = 2 * _size - 1;
    size_t cacheSize = std::min((uint64_t)totalSize, (2UL << _cacheLevel) - 1);
    return Vec::GetMemoryUsage(totalSize,
                               realCacheSize != -1 ? realCacheSize : cacheSize);
  }

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
    int totalLevel = GetLogBaseTwo(leafCount - 1) + 2;
    if ((packLevel == 1 && cacheLevel < totalLevel) ||
        cacheLevel == totalLevel - 1) {
      return GetPathIdx(outputBegin, outputEnd, idx, leafCount, totalLevel);
    }
    int pathLen = (idx | (1UL << totalLevel - 2)) < leafCount ? totalLevel
                                                              : totalLevel - 1;
    auto it = outputBegin;
    int i;
    // set top levels
    for (i = 0; i < std::min(pathLen, cacheLevel); ++i, ++it) {
      PositionType prevNodes = (1UL << i) - 1;
      PositionType subTreeIdx = idx & prevNodes;
      if ((1UL << i) <= leafCount) {
        *it = prevNodes + subTreeIdx;
      } else {
        PositionType mask = 1UL << (totalLevel - 2);
        PositionType deletedNode =
            std::max((int64_t)(std::min(idx, mask) + mask - leafCount), 0L);
        *it = prevNodes + subTreeIdx - deletedNode;
      }
    }
    if (i == pathLen) {
      return pathLen;
    }

    // set fully packed levels in the middle
    for (; i < totalLevel - packLevel - 1; i += packLevel) {
      PositionType prevNodes = (1UL << i) - 1;
      PositionType beginOffset = prevNodes + (idx & prevNodes) * packed_size;
      for (int j = 0; j < packLevel; ++j, ++it) {
        PositionType innerPrevNodes = (1UL << j) - 1;
        *it = beginOffset + innerPrevNodes + ((idx >> i) & innerPrevNodes);
      }
    }

    // set levels with remaining nodes
    PositionType prevNodes = (1UL << i) - 1;
    int remainLevel = totalLevel - i;
    PositionType subTreeIdx = idx & prevNodes;

    // 2 * (number of (totalLevel - 1)-bit elements < size, whose last i bits <
    // subtree idx)
    PositionType prevLeafCount = subTreeIdx * (leafCount >> i) +
                                 std::min(subTreeIdx, leafCount & prevNodes);

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
   * @brief Get the index of the node on the path idx at level level
   *
   * @param idx path idx
   * @param leafCount
   * @param level the level of the node we want to get
   * @param totalLevel
   * @param cacheLevel
   * @return int
   */
  static int GetIdx(PositionType idx, PositionType leafCount, int level,
                    int totalLevel, int cacheLevel) {
    if ((packLevel == 1 && cacheLevel < totalLevel) ||
        cacheLevel == totalLevel - 1) {
      return GetIdx(idx, leafCount, level, totalLevel, totalLevel);
    }

    if (level < cacheLevel) {
      PositionType prevNodes = (1UL << level) - 1;
      PositionType subTreeIdx = idx & prevNodes;
      if ((1UL << level) <= leafCount) {
        return prevNodes + subTreeIdx;
      } else {
        // x < idx and x | (1 << totalLevel - 2) >= size &
        PositionType mask = 1UL << (totalLevel - 2);
        PositionType deletedNode =
            std::max((int64_t)(std::min(idx, mask) + mask - leafCount), 0L);
        return prevNodes + subTreeIdx - deletedNode;
      }
    }

    int i = (level - cacheLevel) / packLevel * packLevel + cacheLevel;

    if (i < totalLevel - packLevel - 1) {
      PositionType prevNodes = (1UL << i) - 1;
      PositionType beginOffset = prevNodes + (idx & prevNodes) * packed_size;
      int j = level - i;
      PositionType innerPrevNodes = (1UL << j) - 1;
      return beginOffset + innerPrevNodes + ((idx >> i) & innerPrevNodes);
    }
    if (i == totalLevel - 1) {
      i -= packLevel;
    }

    PositionType prevNodes = (1UL << i) - 1;
    int remainLevel = totalLevel - i;
    PositionType subTreeIdx = idx & prevNodes;

    // 2 * (number of (totalLevel - 1)-bit elements < size, whose last i bits <
    // subtree idx)
    PositionType prevLeafCount = subTreeIdx * (leafCount >> i) +
                                 std::min(subTreeIdx, leafCount & prevNodes);

    PositionType beginOffset = prevNodes + prevLeafCount * 2 - subTreeIdx;
    // number of (totalLevel - 1)-bit elements < size, whose last i bits =
    // subtree idx
    PositionType subtreeLeafCount = ((leafCount - subTreeIdx - 1) >> i) + 1;
    int subtreeTotalLevel = GetLogBaseTwo(subtreeLeafCount - 1) + 2;
    return beginOffset + GetIdx(idx >> i, subtreeLeafCount, level - i,
                                subtreeTotalLevel, subtreeTotalLevel);
  }

  const T& GetByInternalIdx(PositionType idx) { return arr.Get(idx); }

  const T& GetByPathAndLevel(PositionType path, int level) {
    return GetByInternalIdx(
        GetIdx(path, leafCount, level, totalLevel, cacheLevel));
  }

  T& GetMutableByInternalIdx(PositionType idx) { return arr.GetMutable(idx); }

  T& GetMutableByPathAndLevel(PositionType path, int level) {
    return GetMutableByInternalIdx(
        GetIdx(path, leafCount, level, totalLevel, cacheLevel));
  }

  void SetByInternalIdx(PositionType idx, const T& val) { arr[idx] = val; }

  void SetByPathAndLevel(PositionType path, int level, const T& val) {
    SetByInternalIdx(GetIdx(path, leafCount, level, totalLevel, cacheLevel),
                     val);
  }

  typename Vec::Iterator beginInternal() { return arr.begin(); }

  typename Vec::Iterator endInternal() { return arr.end(); }

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

  template <typename Iterator>
  int ReadSubPath(PositionType pos, Iterator pathBegin, int k) {
    PositionType pathIdx[64];

    int actualLevel = GetPathIdx(&pathIdx[0], &pathIdx[totalLevel], pos,
                                 leafCount, cacheLevel);

    for (int i = k; i < actualLevel; ++i) {
      PositionType idx = pathIdx[i];
      *(pathBegin + i - k) = arr.Get(idx);
    }
    return actualLevel;
  }

  template <typename Iterator>
  int WritePath(PositionType pos, const Iterator pathBegin) {
    PositionType pathIdx[64];
    int actualLevel = GetPathIdx(&pathIdx[0], &pathIdx[totalLevel], pos,
                                 leafCount, cacheLevel);

    for (int i = 0; i < actualLevel; ++i) {
      PositionType idx = pathIdx[i];
      arr[idx] = *(pathBegin + i);
    }
    return actualLevel;
  }

  template <typename Iterator>
  int WriteSubPath(PositionType pos, const Iterator pathBegin, int k) {
    PositionType pathIdx[64];
    int actualLevel = GetPathIdx(&pathIdx[0], &pathIdx[totalLevel], pos,
                                 leafCount, cacheLevel);

    for (int i = k; i < actualLevel; ++i) {
      PositionType idx = pathIdx[i];
      arr[idx] = *(pathBegin + i - k);
    }
    return actualLevel;
  }

  template <typename AggT>
  AggT BuildBottomUpHelper(
      PositionType path, int level,
      const std::function<AggT(T&, const AggT&, const AggT&)>& reduceFunc,
      const std::function<AggT(T&, PositionType)>& leafFunc) {
    if (level == totalLevel - 1 ||
        (level == totalLevel - 2 && (path | (1UL << level)) >= leafCount)) {
      return leafFunc(GetMutableByPathAndLevel(path, level), path);
    }
    const AggT& left =
        BuildBottomUpHelper(path, level + 1, reduceFunc, leafFunc);
    const AggT& right = BuildBottomUpHelper(path | (1UL << level), level + 1,
                                            reduceFunc, leafFunc);
    return reduceFunc(GetMutableByPathAndLevel(path, level), left, right);
  }

  template <typename AggT>
  AggT BuildBottomUp(
      const std::function<AggT(T&, const AggT&, const AggT&)>& reduceFunc,
      const std::function<AggT(T&, PositionType)>& leafFunc) {
    return BuildBottomUpHelper<AggT>(0, 0, reduceFunc, leafFunc);
  }

  void getTopKLevel(int k, T* output) {
    if (k <= cacheLevel) {
      PositionType outputSize = (1UL << k) - 1;
      for (int i = 0; i < outputSize; ++i) {
        *output++ = GetByInternalIdx(i);
      }
      return;
    }
    for (int level = 0; level < k; ++level) {
      for (int i = 0; i < (1UL << level); ++i) {
        *output++ = GetByPathAndLevel(i, level);
      }
    }
  }
};

// struct InverseIncrementer {
//   PositionType num = 0;
//   PositionType highestBitMask = 0;
//   InverseIncrementer(int _length) : highestBitMask(1UL << (_length - 1)) {}
//   InverseIncrementer& operator++() {
//     PositionType mask = highestBitMask;
//     while (num & mask) {
//       num ^= mask;
//       mask >>= 1;
//     }
//     num ^= mask;
//     return *this;
//   }
//   PositionType getIndex() const { return num; }
// };
