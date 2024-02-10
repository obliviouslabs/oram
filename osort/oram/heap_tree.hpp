#pragma once
#include "common/utils.hpp"
#include "external_memory/cachefrontvector.hpp"

/**
 * @brief A tree stored in a heap structure, but pack several levels together,
 * and optimize for reverse lexico order eviction
 *
 * @tparam T
 * @tparam PositionType
 * @tparam page_size The size of the page, default to 4096 bytes
 * @tparam evict_freq Number of reverse lexico evictions per access
 */
template <typename T, typename PositionType = uint64_t,
          const size_t page_size = 4096, const int evict_freq = 2>
struct HeapTree {
  static const int evict_freq_ = evict_freq;
  static constexpr size_t max_node_per_page = divRoundUp(page_size, sizeof(T));
  static constexpr size_t node_per_page_log2 =
      GetLogBaseTwoConstexpr(max_node_per_page) + 1;
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
  static constexpr size_t node_per_page = packed_size << bestLexicoGroupLevel;
  static constexpr size_t actual_page_size = node_per_page * sizeof(T);
  using Vec = EM::CacheFrontVector::Vector<T, actual_page_size, true, true>;
  Vec arr;
  int topLevel = 0;
  int totalLevel = 0;
  PositionType leafCount = 0;
  PositionType cacheSize = 0;
  PositionType extSize = 0;
  PositionType totalSize = 0;

  HeapTree() {}

  HeapTree(PositionType _size, int _topLevel = 62,
           PositionType realCacheSize = -1) {
    Init(_size, _topLevel, realCacheSize);
  }

  void InitWithDefault(PositionType _size, const T& defaultVal,
                       int _topLevel = 62, PositionType realCacheSize = -1) {
    if (totalSize != 0) {
      throw std::runtime_error("Init called on non-empty tree");
    }
    topLevel = _topLevel;
    leafCount = _size;
    totalLevel = GetLogBaseTwo(_size - 1) + 2;
    totalSize = 2 * _size - 1;
    cacheSize = std::min((uint64_t)totalSize, (2UL << _topLevel) - 1);
    arr.SetSize(totalSize, realCacheSize != -1 ? realCacheSize : cacheSize,
                defaultVal);
    extSize = totalSize - cacheSize;
  }

  void Init(PositionType _size, int _topLevel = 62,
            PositionType realCacheSize = -1) {
    InitWithDefault(_size, T(), _topLevel, realCacheSize);
  }

  int GetCacheLevel() const { return topLevel; }

  size_t GetLeafCount() const { return leafCount; }

  size_t GetNodeCount() const { return totalSize; }

  static uint64_t GetMemoryUsage(PositionType _size, int _topLevel = 62,
                                 PositionType realCacheSize = -1) {
    size_t totalSize = 2 * _size - 1;
    size_t cacheSize = std::min((uint64_t)totalSize, (2UL << _topLevel) - 1);
    return Vec::GetMemoryUsage(totalSize,
                               realCacheSize != -1 ? realCacheSize : cacheSize);
  }

  uint64_t GetMemoryUsage() const { return arr.GetMemoryUsage(); }

  template <typename Iterator>
  static int GetPathIdx(Iterator outputBegin, Iterator outputEnd,
                        PositionType idx, PositionType leafCount,
                        int topLevel) {
    int totalLevel = GetLogBaseTwo(leafCount - 1) + 2;
    if ((packLevel == 1 && topLevel < totalLevel) ||
        topLevel == totalLevel - 1) {
      return GetPathIdx(outputBegin, outputEnd, idx, leafCount, totalLevel);
    }
    int pathLen = (idx | (1UL << totalLevel - 2)) < leafCount ? totalLevel
                                                              : totalLevel - 1;
    auto it = outputBegin;
    int i;
    for (i = 0; i < std::min(pathLen, topLevel); ++i, ++it) {
      PositionType prevNodes = (1UL << i) - 1;
      PositionType subTreeIdx = idx & prevNodes;
      if ((1UL << i) <= leafCount) {
        *it = prevNodes + subTreeIdx;
      } else {
        // x < idx and x | (1 << totalLevel - 2) >= size &
        PositionType mask = 1UL << (totalLevel - 2);
        PositionType deletedNode =
            std::max((int64_t)(std::min(idx, mask) + mask - leafCount), 0L);
        *it = prevNodes + subTreeIdx - deletedNode;
      }
    }
    if (i == pathLen) {
      return pathLen;
    }

    for (; i < totalLevel - packLevel - 1; i += packLevel) {
      PositionType prevNodes = (1UL << i) - 1;
      PositionType beginOffset = prevNodes + (idx & prevNodes) * packed_size;
      for (int j = 0; j < packLevel; ++j, ++it) {
        PositionType innerPrevNodes = (1UL << j) - 1;
        *it = beginOffset + innerPrevNodes + ((idx >> i) & innerPrevNodes);
      }
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
    GetPathIdx(it, outputEnd, idx >> i, subtreeLeafCount, remainLevel);
    for (; i < pathLen; ++i, ++it) {
      *it += beginOffset;
    }

    return pathLen;
  }

  static int GetIdx(PositionType idx, PositionType leafCount, int level,
                    int totalLevel, int topLevel) {
    if ((packLevel == 1 && topLevel < totalLevel) ||
        topLevel == totalLevel - 1) {
      return GetIdx(idx, leafCount, level, totalLevel, totalLevel);
    }

    if (level < topLevel) {
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

    int i = (level - topLevel) / packLevel * packLevel + topLevel;

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
        GetIdx(path, leafCount, level, totalLevel, topLevel));
  }

  T& GetMutableByInternalIdx(PositionType idx) { return arr.GetMutable(idx); }

  T& GetMutableByPathAndLevel(PositionType path, int level) {
    return GetMutableByInternalIdx(
        GetIdx(path, leafCount, level, totalLevel, topLevel));
  }

  void SetByInternalIdx(PositionType idx, const T& val) { arr[idx] = val; }

  void SetByPathAndLevel(PositionType path, int level, const T& val) {
    SetByInternalIdx(GetIdx(path, leafCount, level, totalLevel, topLevel), val);
  }

  typename Vec::Iterator beginInternal() { return arr.begin(); }

  typename Vec::Iterator endInternal() { return arr.end(); }

  template <typename Iterator>
  int ReadPath(PositionType pos, Iterator pathBegin) {
    std::vector<PositionType> pathIdx(totalLevel);

    int actualLevel =
        GetPathIdx(pathIdx.begin(), pathIdx.end(), pos, leafCount, topLevel);

    for (int i = 0; i < actualLevel; ++i) {
      PositionType idx = pathIdx[i];
      *(pathBegin + i) = arr.Get(idx);
    }
    return actualLevel;
  }

  template <typename Iterator>
  int WritePath(PositionType pos, const Iterator pathBegin) {
    std::vector<PositionType> pathIdx(totalLevel);
    int actualLevel =
        GetPathIdx(pathIdx.begin(), pathIdx.end(), pos, leafCount, topLevel);

    for (int i = 0; i < actualLevel; ++i) {
      PositionType idx = pathIdx[i];
      arr[idx] = *(pathBegin + i);
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