#pragma once
#include "common/utils.hpp"
template <typename T>
struct HeapTree {
  T* cacheArr = NULL;
  T* extArr = NULL;  // TODO: change it to external
  uint cacheLevel = 0;
  uint totalLevel = 0;
  uint packLevel = 1;
  size_t leafCount;
  size_t cacheSize;
  size_t extSize;
  size_t totalSize;
  bool isPow2 = false;

  HeapTree() = default;
  HeapTree(size_t _size, uint _cacheLevel, uint _packLevel = 1) {
    Init(_size, _cacheLevel, _packLevel);
  }

  ~HeapTree() { Destroy(); }

  void Init(size_t _size, uint _cacheLevel = 62, uint _packLevel = 1) {
    isPow2 = !(_size & (_size - 1));
    leafCount = _size;
    cacheLevel = _cacheLevel;
    packLevel = _packLevel;
    totalSize = 2 * _size - 1;
    totalLevel = GetLogBaseTwo(_size - 1) + 2;
    cacheLevel = std::min(cacheLevel, totalLevel);
    cacheSize = std::min(totalSize, (2UL << cacheLevel) - 1);
    extSize = totalSize - cacheSize;
    // printf("totalSize: %lu, totalLevel: %u\n", totalSize, totalLevel);
    // printf("cacheSize: %lu, extSize: %lu\n", cacheSize, extSize);
    cacheArr = new T[cacheSize];
    if (extSize > 0) {
      extArr = new T[extSize];
    }
  }

  void Destroy() {
    if (cacheArr) {
      delete[] cacheArr;
      cacheArr = NULL;
    }
    if (extArr) {
      delete[] extArr;
      extArr = NULL;
    }
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
    int totalLevel = outputEnd - outputBegin;
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

  T& Get(size_t pos, int level) {
    size_t idx =
        GetCAIdx(pos, cacheSize, level, totalLevel, packLevel, cacheLevel);
    if (idx < cacheSize) {
      return cacheArr[idx];
    } else {
      return extArr[idx - cacheSize];
    }
  }

  void Set(size_t pos, int level, const T& val) {
    size_t idx =
        GetCAIdx(pos, cacheSize, level, totalLevel, packLevel, cacheLevel);
    if (idx < cacheSize) {
      cacheArr[idx] = val;
    } else {
      extArr[idx - cacheSize] = val;
    }
  }

  template <typename Iterator>
  int ReadPath(size_t pos, Iterator pathBegin) const {
    std::vector<size_t> pathIdx(totalLevel);
    // printf("pos: %lu, totalLevel = %d\n", pos, totalLevel);
    // printf("cacheLevel = %d\n", cacheLevel);

    int actualLevel = isPow2 ? GetCAPathIdxPow2(pathIdx.begin(), pathIdx.end(),
                                                pos, packLevel, cacheLevel)
                             : GetCAPathIdx(pathIdx.begin(), pathIdx.end(), pos,
                                            leafCount, packLevel, cacheLevel);

    for (int i = 0; i < actualLevel; ++i) {
      size_t idx = pathIdx[i];
      if (idx < cacheSize) {
        *(pathBegin + i) = cacheArr[idx];
      } else {
        *(pathBegin + i) = extArr[idx - cacheSize];
      }
    }
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
      if (idx < cacheSize) {
        cacheArr[idx] = *(pathBegin + i);
      } else {
        extArr[idx - cacheSize] = *(pathBegin + i);
      }
    }
    return actualLevel;
  }

  // struct ReverseLexLeafIterator {
  //   std::vector<size_t> path;
  //   size_t idx;
  //   size_t leafCount;

  // };

  struct ReverseLexLeafPow2Indexer {
    size_t idx;
    size_t pos;
    std::vector<size_t> rzero2inc;

    ReverseLexLeafPow2Indexer(int totalLevel, int cacheLevel, int packLevel) {
      std::vector<int> topLevels;
      int botLevel = totalLevel - cacheLevel;
      if (botLevel <= 0) {
        botLevel = totalLevel;
      } else {
        topLevels.push_back(cacheLevel);
      }
      while (botLevel > 1) {
        size_t topLevel = botLevel >> 1;
        if (topLevel > packLevel && packLevel > 1) {
          topLevel = topLevel / packLevel * packLevel;
        }
        topLevels.push_back(topLevel);
        botLevel -= topLevel;
      }
      size_t accInc = 1;

      auto topIter = topLevels.rbegin();
      size_t accRZero = *topIter;
      for (int i = 0; i < totalLevel; ++i) {
        if (i == accRZero) {
          accInc += (1UL << (*topIter)) - 1;
          accRZero += *(++topIter);
        }
        rzero2inc.push_back(accInc);
      }
      // for (int i = 0; i < totalLevel; ++i) {
      //   printf("%lu ", rzero2inc[i]);
      // }
      idx = rzero2inc.back() - 1;
      pos = 0;
    }

    ReverseLexLeafPow2Indexer(const HeapTree& tree)
        : ReverseLexLeafPow2Indexer(tree.totalLevel, tree.cacheLevel,
                                    tree.packLevel) {}

    ReverseLexLeafPow2Indexer& operator++() {
      int rZero = 0;
      int tmpPos = ++pos;
      while (!(tmpPos & 1)) {
        tmpPos >>= 1;
        ++rZero;
      }
      idx += rzero2inc[rZero];
      return *this;
    }

    size_t getIndex() const { return idx; }
  };

  struct ReverseLexLeafIndexer : public ReverseLexLeafPow2Indexer {
    typedef ReverseLexLeafPow2Indexer Base;
    size_t highestBitMask;
    size_t size;
    bool separateLastLevel = false;
    ReverseLexLeafPow2Indexer prevLevelIndexer;
    size_t lastLevelIdx = 0;
    ReverseLexLeafIndexer(const HeapTree& tree)
        : ReverseLexLeafIndexer(tree.totalLevel, tree.cacheLevel,
                                tree.packLevel, tree.leafCount) {}

    ReverseLexLeafIndexer(int totalLevel, int cacheLevel, int packLevel,
                          size_t _size)
        : Base(totalLevel, cacheLevel, packLevel),
          prevLevelIndexer(std::max(2, totalLevel - 1), cacheLevel, packLevel) {
      separateLastLevel =
          cacheLevel >= 2 && totalLevel == cacheLevel + 1;  // special case
      highestBitMask = 1UL << (totalLevel - 2);
      size = _size;
      lastLevelIdx = Base::idx;
    }

    ReverseLexLeafIndexer& operator++() {
      size_t mask = highestBitMask;
      int consecutiveOnes = 0;
      while (Base::pos & mask) {
        consecutiveOnes++;
        Base::pos ^= mask;
        mask >>= 1;
      }
      Base::pos ^= mask;
      size_t posSibling = Base::pos | highestBitMask;
      if (separateLastLevel) {
        // special case
        // printf("posSibling: %lu, size: %lu\n", posSibling, size);
        if (posSibling >= size) {
          Base::idx = prevLevelIndexer.getIndex();
          // printf("idx becomes %lu\n", Base::idx);
          Base::pos = posSibling;
        } else {
          Base::idx = ++lastLevelIdx;
        }
        if (Base::pos & highestBitMask) {
          // printf("increasing prevLevelIndexer\n");
          ++prevLevelIndexer;
        }
        return *this;
      }
      if (posSibling >= size) {
        Base::idx += Base::rzero2inc[consecutiveOnes] - 1;
        Base::pos = posSibling;
      } else {
        Base::idx += Base::rzero2inc[consecutiveOnes];
      }
      return *this;
    }

    size_t getIndex() const { return Base::idx; }
  };
};