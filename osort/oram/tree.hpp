#pragma once
#include "common/utils.hpp"
#include "external_memory/cachefrontvector.hpp"
template <typename T, typename PositionType = uint64_t>
struct HeapTree {
  using Vec = EM::CacheFrontVector::Vector<T>;
  Vec arr;
  uint cacheLevel = 0;
  uint totalLevel = 0;
  uint packLevel = 1;
  size_t leafCount;
  size_t cacheSize;
  size_t extSize;
  size_t totalSize;
  bool isPow2 = false;

  HeapTree(size_t _size, uint _cacheLevel = 62, uint _packLevel = 1,
           size_t realCacheSize = -1)
      : cacheLevel(_cacheLevel),
        packLevel(_packLevel),
        leafCount(_size),
        totalLevel(GetLogBaseTwo(_size - 1) + 2),
        totalSize(2 * _size - 1),
        cacheSize(std::min(totalSize, (2UL << _cacheLevel) - 1)),
        arr(2 * _size - 1,
            realCacheSize != -1 ? realCacheSize : (2UL << _cacheLevel) - 1),
        isPow2(!(_size & (_size - 1))) {
    extSize = totalSize - cacheSize;
  }

  static size_t GetCAIdxPow2(size_t idx, int level, int totalLevel,
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
    size_t bottomIdx = reverseBits(idx, topLevel);
    size_t offset = (1UL << topLevel) - 1 + bottomIdx * ((1UL << botLevel) - 1);

    return offset +
           GetCAIdxPow2(idx >> topLevel, level - topLevel, botLevel, packLevel);
  }

  static size_t getCABottomOffset(size_t idx, size_t leafCount, int topLevel,
                                  size_t& bottomLeafCountMutable) {
    size_t offset = 0;
    size_t numRight = 1;
    size_t nodeCount = 2 * leafCount - 1;
    for (int i = 0; i < topLevel; ++i) {
      size_t leftCount = ((nodeCount >> 2) << 1) + 1;
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

  static size_t GetCAIdx(size_t idx, size_t leafCount, int level,
                         int totalLevel, int packLevel = 1, int topLevel = 0) {
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
    size_t bottomLeafCount;
    size_t offset =
        getCABottomOffset(idx, leafCount, topLevel, bottomLeafCount);

    return offset + GetCAIdx(idx >> topLevel, bottomLeafCount, level - topLevel,
                             botLevel, packLevel);
  }

  template <typename Iterator>
  static int GetCAPathIdxPow2(Iterator outputBegin, Iterator outputEnd,
                              size_t idx, int packLevel = 1, int topLevel = 0) {
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
    size_t bottomIdx = reverseBits(idx, topLevel);
    size_t offset = (1UL << topLevel) - 1 + bottomIdx * ((1UL << botLevel) - 1);
    GetCAPathIdxPow2(outputBegin + topLevel, outputEnd, idx >> topLevel,
                     packLevel);
    for (auto it = outputBegin + topLevel; it != outputEnd; ++it) {
      *it += offset;
    }
    return totalLevel;
  }

  template <typename Iterator>
  static int GetCAPathIdx(Iterator outputBegin, Iterator outputEnd, size_t idx,
                          size_t leafCount, int packLevel = 1,
                          int topLevel = 0) {
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
    // printf("idx: %lu, leafCount: %lu, topLevel: %d, packLevel: %d\n", idx,
    //        leafCount, topLevel, packLevel);

    int botLevel = totalLevel - topLevel;

    GetCAPathIdxPow2(outputBegin, outputBegin + topLevel, idx, packLevel);
    if (botLevel == 1 && (idx | (1UL << (topLevel - 1))) >= leafCount) {
      return topLevel;
    }
    // the idx of the bottom subtree
    size_t bottomLeafCount;

    size_t offset =
        getCABottomOffset(idx, leafCount, topLevel, bottomLeafCount);
    // printf(
    //     "idx: %lu, leafCount: %lu, topLevel: %d, offset: %lu,
    //     bottomLeafCount: "
    //     "%lu\n",
    //     idx, leafCount, topLevel, offset, bottomLeafCount);

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

  const T& Get(size_t pos, int level) const {
    size_t idx =
        GetCAIdx(pos, cacheSize, level, totalLevel, packLevel, cacheLevel);
    return GetByInternalIdx(idx);
  }

  void Set(size_t pos, int level, const T& val) {
    size_t idx =
        GetCAIdx(pos, cacheSize, level, totalLevel, packLevel, cacheLevel);
    SetByInternalIdx(idx, val);
  }

  const T& GetByInternalIdx(size_t idx) { return arr.Get(idx); }

  void SetByInternalIdx(size_t idx, const T& val) { arr[idx] = val; }

  typename Vec::Iterator beginInteranl() { return arr.begin(); }

  typename Vec::Iterator endInteranl() { return arr.end(); }

  template <typename Iterator>
  int ReadPath(size_t pos, Iterator pathBegin) {
    std::vector<size_t> pathIdx(totalLevel);
    // printf("pos: %lu, totalLevel = %d\n", pos, totalLevel);
    // printf("cacheLevel = %d\n", cacheLevel);

    int actualLevel = isPow2 ? GetCAPathIdxPow2(pathIdx.begin(), pathIdx.end(),
                                                pos, packLevel, cacheLevel)
                             : GetCAPathIdx(pathIdx.begin(), pathIdx.end(), pos,
                                            leafCount, packLevel, cacheLevel);
    // printf("Read path with internal index\n");
    for (int i = 0; i < actualLevel; ++i) {
      size_t idx = pathIdx[i];
      // printf("%lu ", idx);
      *(pathBegin + i) = arr.Get(idx);
    }
    // printf("\n");
    return actualLevel;
  }

  template <typename Iterator>
  int WritePath(size_t pos, const Iterator pathBegin) {
    std::vector<size_t> pathIdx(totalLevel);
    // printf("pos: %lu, totalLevel = %d\n", pos, totalLevel);
    int actualLevel = isPow2 ? GetCAPathIdxPow2(pathIdx.begin(), pathIdx.end(),
                                                pos, packLevel, cacheLevel)
                             : GetCAPathIdx(pathIdx.begin(), pathIdx.end(), pos,
                                            leafCount, packLevel, cacheLevel);

    // for (int i = 0; i < actualLevel; ++i) {
    //   printf("%lu ", pathIdx[i]);
    // }
    // printf("\n");
    for (int i = 0; i < actualLevel; ++i) {
      size_t idx = pathIdx[i];
      arr[idx] = *(pathBegin + i);
    }
    return actualLevel;
  }

  template <typename AggT>
  AggT BuildBottomUp(
      const std::function<AggT(T&, const AggT&, const AggT&)>& reduceFunc,
      const std::function<AggT(T&, size_t)>& leafFunc) {
    return BuildBottomUp<AggT>(arr.begin(), arr.end(), reduceFunc, leafFunc,
                               cacheLevel, packLevel);
  }

  template <typename AggT, typename Iterator>
  static AggT BuildBottomUp(
      Iterator begin, Iterator end,
      const std::function<AggT(T&, const AggT&, const AggT&)>& reduceFunc,
      const std::function<AggT(T&, size_t)>& leafFunc, uint _cacheLevel = 62,
      uint _packLevel = 1) {
    size_t size = end - begin;
    Assert(size & 1);
    int totalLevel = GetLogBaseTwo(size) + 1;
    int topLevel = _cacheLevel >= totalLevel ? 0 : _cacheLevel;
    return BuildBottomUpHelper<AggT>(begin, end, 0, 0, reduceFunc, leafFunc,
                                     topLevel, _packLevel);
  };

  template <typename AggT, typename Iterator>
  static AggT BuildBottomUpHelper(
      Iterator begin, Iterator end, size_t path, int level,
      const std::function<AggT(T&, const AggT&, const AggT&)>& reduceFunc,
      const std::function<AggT(T&, size_t)>& leafFunc, uint topLevel,
      uint packLevel) {
    size_t size = end - begin;
    Assert(size & 1);
    size_t leafCount = (size + 1) / 2;
    // printf("size = %lu, leafCount = %lu\n", size, leafCount);
    if (size == 1) {
      // printf("leafFunc due to size = 1\n");
      return leafFunc(*begin, path);
    } else if (size == 3) {
      // printf("leafFunc due to size = 3\n");
      return reduceFunc(*begin, leafFunc(*(begin + 1), path),
                        leafFunc(*(begin + 2), path | (1UL << level)));
    }
    int totalLevel = GetLogBaseTwo(size) + 1;
    if (topLevel == 0) {
      topLevel = totalLevel >> 1;
      if (topLevel > packLevel && packLevel > 1) {
        topLevel = topLevel / packLevel * packLevel;
      }
    }
    int botLevel = totalLevel - topLevel;
    size_t topSize = (1UL << topLevel) - 1;
    // InverseIncrementer topLeafIndexer(topLevel - 1);
    size_t botOffset = topSize;
    size_t topLeafCount = 1UL << (topLevel - 1);
    int leafLevel = level + topLevel;
    std::function<AggT(T&, size_t)> topLeafFunc = [&](T& leaf, size_t path) {
      // size_t topLeafIdx = topLeafIndexer.getIndex();
      // ++topLeafIndexer;  // increment topLeafIdx after every leaf call
      size_t topLeafIdx = path >> level;
      // printf("topLeafIdx = %lu\n", topLeafIdx);
      size_t left1 = leafCount - topLeafIdx - 1;
      if (left1 < topLeafCount) {
        // printf("leafFunc due to lingering node\n");
        return leafFunc(leaf, path);
      }
      size_t right1 = left1 - topLeafCount;
      size_t leftLeafCount = ((leafCount - topLeafIdx - 1) >> topLevel) + 1;
      size_t rightLeafCount =
          ((leafCount - topLeafIdx - topLeafCount - 1) >> topLevel) + 1;
      size_t leftSize = 2 * leftLeafCount - 1;
      size_t rightSize = 2 * rightLeafCount - 1;

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

struct InverseIncrementer {
  size_t num = 0;
  size_t highestBitMask = 0;
  InverseIncrementer(int _length) : highestBitMask(1UL << (_length - 1)) {}
  InverseIncrementer& operator++() {
    size_t mask = highestBitMask;
    while (num & mask) {
      num ^= mask;
      mask >>= 1;
    }
    num ^= mask;
    return *this;
  }
  size_t getIndex() const { return num; }
};

// template <typename Iterator,
//           typename T = typename std::iterator_traits<Iterator>::value_type>
// T& BuildBottomUp(Iterator begin, Iterator end,
//                  const std::function<T(T&, T&)>& reduceFunc,
//                  const std::function<void(T&)>& leafFunc, uint _cacheLevel
//                  = 62, uint _packLevel = 1) {
//   size_t size = end - begin;
//   Assert(size & 1);
//   int totalLevel = GetLogBaseTwo(size) + 1;
//   int topLevel = _cacheLevel >= totalLevel ? 0 : _cacheLevel;
//   return BuildBottomUpHelper<T>(begin, end, reduceFunc, leafFunc, topLevel,
//                                 _packLevel);
// };
