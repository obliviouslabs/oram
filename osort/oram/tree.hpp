#pragma once
#include "common/utils.hpp"
#include "external_memory/cachefrontvector.hpp"
template <typename T>
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

  // struct ReverseLexLeafIterator {
  //   std::vector<size_t> path;
  //   size_t idx;
  //   size_t leafCount;

  // };

  struct ReverseLexLeafPow2Indexer {
    size_t idx;
    size_t pos;
    std::vector<size_t> rzero2inc;
    int totalLevel;

    ReverseLexLeafPow2Indexer(int totalLevel, int cacheLevel, int packLevel)
        : totalLevel(totalLevel) {
      if (totalLevel <= 2) {
        idx = totalLevel - 1;
        return;
      }
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
      if (totalLevel <= 2) {
        ++idx;
        return *this;
      }
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

  template <const bool isPow2 = false>
  struct LeafReader {
    std::conditional_t<isPow2, HeapTree::ReverseLexLeafPow2Indexer,
                       HeapTree::ReverseLexLeafIndexer>
        iter;
    const HeapTree& tree;
    size_t idx = 0;
    LeafReader(const HeapTree& _tree) : iter(_tree), tree(_tree) {}

    const T& get() const { return tree.GetByInternalIdx(iter.getIndex()); }

    const T& read() {
      const T& res = get();
      ++iter;
      ++idx;
      return res;
    }

    bool eof() { return idx >= tree.leafCount; }

    size_t size() { return tree.leafCount; }
  };

  template <const bool isPow2 = false>
  struct LeafWriter {
    std::conditional_t<isPow2, HeapTree::ReverseLexLeafPow2Indexer,
                       HeapTree::ReverseLexLeafIndexer>
        iter;
    HeapTree& tree;
    size_t idx = 0;
    LeafWriter(HeapTree& _tree) : iter(_tree), tree(_tree) {}

    void write(const T& val) {
      tree.SetByInternalIdx(iter.getIndex(), val);
      ++iter;
      ++idx;
    }

    bool eof() { return idx >= tree.leafCount; }

    void flush() {}

    size_t size() { return tree.leafCount; }
  };
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
