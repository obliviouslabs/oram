#pragma once
#include "common/utils.hpp"
#include "external_memory/cachefrontvector.hpp"
template <typename T, typename PositionType = uint64_t>
struct HeapTree {
  using Vec = EM::CacheFrontVector::Vector<T, 4096, true, true>;
  Vec arr;
  int cacheLevel = 0;
  int totalLevel = 0;
  int packLevel = 1;
  PositionType leafCount = 0;
  PositionType cacheSize = 0;
  PositionType extSize = 0;
  PositionType totalSize = 0;
  bool isPow2 = false;

  HeapTree() {}

  HeapTree(PositionType _size, int _cacheLevel = 62, int _packLevel = 1,
           PositionType realCacheSize = -1) {
    Init(_size, _cacheLevel, _packLevel, realCacheSize);
  }

  void InitWithDefault(PositionType _size, const T& defaultVal,
                       int _cacheLevel = 62, int _packLevel = 1,
                       PositionType realCacheSize = -1) {
    if (totalSize != 0) {
      throw std::runtime_error("Init called on non-empty tree");
    }
    cacheLevel = _cacheLevel;
    packLevel = _packLevel;
    leafCount = _size;
    totalLevel = GetLogBaseTwo(_size - 1) + 2;
    totalSize = 2 * _size - 1;
    cacheSize = std::min((uint64_t)totalSize, (2UL << _cacheLevel) - 1);
    arr.SetSize(totalSize, realCacheSize != -1 ? realCacheSize : cacheSize,
                defaultVal);
    isPow2 = !(_size & (_size - 1));
    extSize = totalSize - cacheSize;
  }

  void Init(PositionType _size, int _cacheLevel = 62, int _packLevel = 1,
            PositionType realCacheSize = -1) {
    InitWithDefault(_size, T(), _cacheLevel, _packLevel, realCacheSize);
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

  static PositionType GetCAIdxPow2(PositionType idx, int level, int totalLevel,
                                   int packLevel = 1, int topLevel = 0) {
    if (totalLevel == 1) {
      return 0;
    }
    if (topLevel == 0) {
      topLevel = totalLevel >> 1;
      if (topLevel > packLevel && packLevel > 1) {
        topLevel = topLevel / packLevel * packLevel;
      }
    }
    int botLevel = totalLevel - topLevel;

    if (level < topLevel) {
      return GetCAIdxPow2(idx, level, topLevel, packLevel);
    }
    // the idx of the bottom subtree
    PositionType bottomIdx = reverseBits(idx, topLevel);
    PositionType offset =
        (1UL << topLevel) - 1 + bottomIdx * ((1UL << botLevel) - 1);

    return offset +
           GetCAIdxPow2(idx >> topLevel, level - topLevel, botLevel, packLevel);
  }

  static PositionType getCABottomOffset(PositionType idx,
                                        PositionType leafCount, int topLevel,
                                        PositionType& bottomLeafCountMutable) {
    PositionType offset = 0;
    PositionType numRight = 1;
    PositionType nodeCount = 2 * leafCount - 1;
    for (int i = 0; i < topLevel; ++i) {
      PositionType leftCount = ((nodeCount >> 2) << 1) + 1;
      offset += numRight;
      if ((idx >> i) & 1) {
        offset += leftCount;
        nodeCount = nodeCount - leftCount - 1;
        numRight = 2 * numRight - 1;
      } else {
        nodeCount = leftCount;
        numRight = 2 * numRight;
      }
    }

    bottomLeafCountMutable = (nodeCount + 1) >> 1;
    return offset;
  }

  static PositionType GetCAIdx(PositionType idx, PositionType leafCount,
                               int level, int totalLevel, int packLevel = 1,
                               int topLevel = 0) {
    if (totalLevel == 1) {
      return 0;
    }
    if (topLevel == 0) {
      topLevel = totalLevel >> 1;
      if (topLevel > packLevel && packLevel > 1) {
        topLevel = topLevel / packLevel * packLevel;
      }
    }
    int botLevel = totalLevel - topLevel;

    if (level < topLevel) {
      return GetCAIdxPow2(idx, level, topLevel, packLevel);
    }
    // the idx of the bottom subtree
    PositionType bottomLeafCount;
    PositionType offset =
        getCABottomOffset(idx, leafCount, topLevel, bottomLeafCount);

    return offset + GetCAIdx(idx >> topLevel, bottomLeafCount, level - topLevel,
                             botLevel, packLevel);
  }

  template <typename Iterator>
  static int GetCAPathIdxPow2(Iterator outputBegin, Iterator outputEnd,
                              PositionType idx, int packLevel = 1,
                              int topLevel = 0) {
    int totalLevel = outputEnd - outputBegin;
    if (totalLevel == 1) {
      *outputBegin = 0;
      return 1;
    }
    if (topLevel == 0) {
      topLevel = totalLevel >> 1;
      if (topLevel > packLevel && packLevel > 1) {
        topLevel = topLevel / packLevel * packLevel;
      }
    } else if (topLevel >= totalLevel) {
      return GetCAPathIdxPow2(outputBegin, outputEnd, idx, packLevel);
    }
    int botLevel = totalLevel - topLevel;

    GetCAPathIdxPow2(outputBegin, outputBegin + topLevel, idx, packLevel);

    // the idx of the bottom subtree
    PositionType bottomIdx = reverseBits(idx, topLevel);
    PositionType offset =
        (1UL << topLevel) - 1 + bottomIdx * ((1UL << botLevel) - 1);
    GetCAPathIdxPow2(outputBegin + topLevel, outputEnd, idx >> topLevel,
                     packLevel);
    for (auto it = outputBegin + topLevel; it != outputEnd; ++it) {
      *it += offset;
    }
    return totalLevel;
  }

  template <typename Iterator>
  static int GetCAPathIdx(Iterator outputBegin, Iterator outputEnd,
                          PositionType idx, PositionType leafCount,
                          int packLevel = 1, int topLevel = 0) {
    Assert(leafCount > 0);
    if (leafCount <= 1) {
      *outputBegin = 0;
      return 1;
    }
    if (leafCount == 2) {
      *outputBegin = 0;
      *(outputBegin + 1) = 1 + (idx & 1);
      return 2;
    }
    int totalLevel = GetLogBaseTwo(leafCount - 1) + 2;
    if (topLevel == 0) {
      topLevel = totalLevel >> 1;
      if (topLevel > packLevel && packLevel > 1) {
        topLevel = topLevel / packLevel * packLevel;
      }
    } else if (topLevel >= totalLevel) {
      return GetCAPathIdx(outputBegin, outputEnd, idx, leafCount, packLevel);
    }
    int botLevel = totalLevel - topLevel;

    GetCAPathIdxPow2(outputBegin, outputBegin + topLevel, idx, packLevel);
    if (botLevel == 1 && (idx | (1UL << (topLevel - 1))) >= leafCount) {
      return topLevel;
    }
    // the idx of the bottom subtree
    PositionType bottomLeafCount;

    PositionType offset =
        getCABottomOffset(idx, leafCount, topLevel, bottomLeafCount);

    int actualBottomLevel =
        GetCAPathIdx(outputBegin + topLevel, outputEnd, idx >> topLevel,
                     bottomLeafCount, packLevel);
    int actualLevel = topLevel + actualBottomLevel;
    for (auto it = outputBegin + topLevel; it != outputBegin + actualLevel;
         ++it) {
      *it += offset;
    }
    return actualLevel;
  }

  const T& Get(PositionType pos, int level) const {
    PositionType idx =
        GetCAIdx(pos, cacheSize, level, totalLevel, packLevel, cacheLevel);
    return GetByInternalIdx(idx);
  }

  void Set(PositionType pos, int level, const T& val) {
    PositionType idx =
        GetCAIdx(pos, cacheSize, level, totalLevel, packLevel, cacheLevel);
    SetByInternalIdx(idx, val);
  }

  const T& GetByInternalIdx(PositionType idx) { return arr.Get(idx); }

  T& GetMutableByInternalIdx(PositionType idx) { return arr.GetMutable(idx); }

  void SetByInternalIdx(PositionType idx, const T& val) { arr[idx] = val; }

  typename Vec::Iterator beginInteranl() { return arr.begin(); }

  typename Vec::Iterator endInteranl() { return arr.end(); }

  template <typename Iterator>
  int ReadPath(PositionType pos, Iterator pathBegin) {
    std::vector<PositionType> pathIdx(totalLevel);

    int actualLevel = isPow2 ? GetCAPathIdxPow2(pathIdx.begin(), pathIdx.end(),
                                                pos, packLevel, cacheLevel)
                             : GetCAPathIdx(pathIdx.begin(), pathIdx.end(), pos,
                                            leafCount, packLevel, cacheLevel);

    for (int i = 0; i < actualLevel; ++i) {
      PositionType idx = pathIdx[i];
      *(pathBegin + i) = arr.Get(idx);
    }
    return actualLevel;
  }

  template <typename Iterator>
  int WritePath(PositionType pos, const Iterator pathBegin) {
    std::vector<PositionType> pathIdx(totalLevel);
    int actualLevel = isPow2 ? GetCAPathIdxPow2(pathIdx.begin(), pathIdx.end(),
                                                pos, packLevel, cacheLevel)
                             : GetCAPathIdx(pathIdx.begin(), pathIdx.end(), pos,
                                            leafCount, packLevel, cacheLevel);

    for (int i = 0; i < actualLevel; ++i) {
      PositionType idx = pathIdx[i];
      arr[idx] = *(pathBegin + i);
    }
    return actualLevel;
  }

  template <typename AggT>
  AggT BuildBottomUp(
      const std::function<AggT(T&, const AggT&, const AggT&)>& reduceFunc,
      const std::function<AggT(T&, PositionType)>& leafFunc) {
    return BuildBottomUp<AggT>(arr.begin(), arr.end(), reduceFunc, leafFunc,
                               cacheLevel, packLevel);
  }

  template <typename AggT, typename Iterator>
  static AggT BuildBottomUp(
      Iterator begin, Iterator end,
      const std::function<AggT(T&, const AggT&, const AggT&)>& reduceFunc,
      const std::function<AggT(T&, PositionType)>& leafFunc,
      int _cacheLevel = 62, int _packLevel = 1) {
    PositionType size = end - begin;
    Assert(size & 1);
    int totalLevel = GetLogBaseTwo(size) + 1;

    int topLevel = _cacheLevel >= totalLevel ? 0 : _cacheLevel;
    return BuildBottomUpHelper<AggT>(begin, end, 0, 0, reduceFunc, leafFunc,
                                     topLevel, _packLevel);
  };

  template <typename AggT, typename Iterator>
  static AggT BuildBottomUpHelper(
      Iterator begin, Iterator end, PositionType path, int level,
      const std::function<AggT(T&, const AggT&, const AggT&)>& reduceFunc,
      const std::function<AggT(T&, PositionType)>& leafFunc, int topLevel,
      int packLevel) {
    PositionType size = end - begin;
    Assert(size & 1);
    PositionType leafCount = (size + 1) / 2;

    if (size == 1) {
      return leafFunc(*begin, path);
    } else if (size == 3) {
      const AggT& left = leafFunc(*(begin + 1), path);
      const AggT& right = leafFunc(*(begin + 2), path | (1UL << level));
      return reduceFunc(*begin, left, right);
    }
    int totalLevel = GetLogBaseTwo(size) + 1;
    if (topLevel == 0) {
      topLevel = totalLevel >> 1;
      if (topLevel > packLevel && packLevel > 1) {
        topLevel = topLevel / packLevel * packLevel;
      }
    }

    int botLevel = totalLevel - topLevel;
    PositionType topSize = (1UL << topLevel) - 1;
    PositionType botOffset = topSize;
    PositionType topLeafCount = 1UL << (topLevel - 1);
    int leafLevel = level + topLevel;

    auto topLeafFunc = [=, &reduceFunc, &leafFunc, &botOffset](
                           T& leaf, PositionType path) {
      PositionType topLeafIdx = path >> level;
      PositionType left1 = leafCount - topLeafIdx - 1;
      if (left1 < topLeafCount) {
        // a path of depth-1
        return leafFunc(leaf, path);
      }
      PositionType right1 = left1 - topLeafCount;
      PositionType leftLeafCount =
          ((leafCount - topLeafIdx - 1) >> topLevel) + 1;
      PositionType rightLeafCount =
          ((leafCount - topLeafIdx - topLeafCount - 1) >> topLevel) + 1;
      PositionType leftSize = 2 * leftLeafCount - 1;
      PositionType rightSize = 2 * rightLeafCount - 1;

      AggT leftResult = BuildBottomUpHelper<AggT, Iterator>(
          begin + botOffset, begin + botOffset + leftSize, path, leafLevel,
          reduceFunc, leafFunc, 0, packLevel);
      botOffset += leftSize;
      AggT rightResult = BuildBottomUpHelper<AggT, Iterator>(
          begin + botOffset, begin + botOffset + rightSize,
          path | (1UL << leafLevel) / 2, leafLevel, reduceFunc, leafFunc, 0,
          packLevel);
      botOffset += rightSize;
      return reduceFunc(leaf, leftResult, rightResult);
    };
    return BuildBottomUpHelper<AggT, Iterator>(begin, begin + topSize, path,
                                               level, reduceFunc, topLeafFunc,
                                               0, packLevel);
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
