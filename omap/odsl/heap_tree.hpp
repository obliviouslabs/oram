#pragma once
#include "common/utils.hpp"
#include "external_memory/cachefrontvector.hpp"
#define MAX_CACHE_LEVEL 62
#define MAX_CACHE_SIZE (1UL << MAX_CACHE_LEVEL)
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
  // The actual data array, let freshness be checked at the application level
  using Vec = EM::CacheFrontVector::Vector<
      T, actual_page_size, EM::CacheFrontVector::EncryptType::ENCRYPT_AND_AUTH>;
  Vec arr;  // the actual data
  int cacheLevel = 0;
  int totalLevel = 0;
  PositionType leafCount = 0;
  PositionType cacheSize = 0;
  PositionType extSize = 0;
  PositionType totalSize = 0;  // total number of nodes

 public:
  HeapTree() {}

  explicit HeapTree(PositionType _size, int _cacheLevel = MAX_CACHE_LEVEL) {
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
                       int _cacheLevel = MAX_CACHE_LEVEL) {
    if (totalSize != 0) {
      throw std::runtime_error("Init called on non-empty tree");
    }
    cacheLevel = _cacheLevel;
    leafCount = _size;
    totalLevel = GetLogBaseTwo(_size - 1) + 2;
    if (totalLevel > 64) {
      throw std::runtime_error("ORAM size too large.");
    }
    totalSize = 2 * _size - 1;
    cacheSize =
        (PositionType)std::min((uint64_t)totalSize, (2UL << _cacheLevel) - 1);
    arr.SetSize(totalSize, cacheSize, defaultVal);
    extSize = totalSize - cacheSize;
  }

  void Init(PositionType _size, int _cacheLevel = MAX_CACHE_LEVEL) {
    InitWithDefault(_size, T(), _cacheLevel);
  }

  /**
   * @brief Get the memory usage of the tree
   *
   * @param _size Number of leaves in the tree
   * @param _cacheLevel The number of top levels to cache.
   * @return uint64_t The memory usage in bytes
   */
  static uint64_t GetMemoryUsage(PositionType _size,
                                 int _cacheLevel = MAX_CACHE_LEVEL) {
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
   * @brief Return whether the tree is fully cached
   *
   * @return true if the tree is fully cached
   */
  INLINE bool IsFullyCached() const { return cacheSize == totalSize; }

  /**
   * @brief Retrieve the node indices on a path idx
   *
   * @tparam Iterator
   * @param outputBegin
   * @param idx idx of the path
   * @param leafCount the number of leaves in the tree
   * @param totalLevel the total number of levels in the tree
   * @param cacheLevel the number of top levels to cache
   * @return int the number of nodes in the path
   */
  template <class Iterator>
  static int GetNodeIdxArr(Iterator outputBegin, PositionType idx,
                           PositionType leafCount, int totalLevel,
                           int cacheLevel) {
    // when all the levels are cached
    if ((packLevel == 1 && cacheLevel < totalLevel) ||
        cacheLevel == totalLevel - 1) {
      return GetNodeIdxArr(outputBegin, idx, leafCount, totalLevel, totalLevel);
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
        PositionType leafCountPlusDeleted = std::min(idx, mask) + mask;
        PositionType deletedNode = leafCountPlusDeleted >= leafCount
                                       ? leafCountPlusDeleted - leafCount
                                       : 0;
        *it = prevNodes + subTreeIdx - (PositionType)deletedNode;
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
    GetNodeIdxArr(it, idx >> i, subtreeLeafCount, remainLevel, remainLevel);
    for (; i < pathLen; ++i, ++it) {
      *it += beginOffset;
    }

    return pathLen;
  }

  /**
   * @brief Get indices of the nodes on the path
   *
   * @tparam Iterator The type of the iterator to store the path indices
   * @param outputBegin The begin iterator to store the path indices
   * @param idx  The index of the path
   * @return int The number of nodes in the path
   */
  template <class Iterator>
  int GetNodeIdxArr(Iterator outputBegin, PositionType idx) const {
    return GetNodeIdxArr(outputBegin, idx, leafCount, totalLevel, cacheLevel);
  }

  /**
   * @brief Get individual node By index
   *
   * @param idx The index of the node
   * @return T& The reference to the node
   */
  T& GetNodeByIdx(PositionType idx) { return arr[idx]; }

  struct BatchAccessor {
    HeapTree& tree;
    typename Vec::BatchAccessor arrAccessor;

    BatchAccessor(HeapTree& _tree) : tree(_tree), arrAccessor(_tree.arr) {}

    template <class Iterator>
    int GetNodeIdxArrAndPrefetch(Iterator outputBegin, PositionType idx,
                                 uint32_t& pathPrefetchReceipt) {
      int pathLen = tree.GetNodeIdxArr(outputBegin, idx);
      uint32_t levelReceipt = -1;
      // printf("read path %d with %d nodes\n", (int)idx, pathLen);
      for (int i = 0; i < pathLen; ++i) {
        levelReceipt = arrAccessor.Prefetch(outputBegin[i]);
      }
      // the receipt of cached levels can be arbitrary
      pathPrefetchReceipt = levelReceipt - (pathLen - 1);
      return pathLen;
    };

    void FlushRead() { arrAccessor.FlushRead(); }

    T& GetPrefetchedNode(PositionType idx, int level,
                         uint32_t pathPrefetchReceipt) {
      // printf("get node %d at level %d with receipt %d\n", (int)idx, level,
      //  (int)pathPrefetchReceipt + level);
      return arrAccessor.At(pathPrefetchReceipt + level, idx);
    }

    void FlushWrite() { arrAccessor.FlushWrite(); }
  };
};