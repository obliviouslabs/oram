#pragma once

#include "oram_common.hpp"

/// @brief This file implements Circuit ORAM (https://eprint.iacr.org/2014/672),
/// which turns out to be faster than Path ORAM when full obliviousness is
/// desired.

namespace ODSL::CircuitORAM {
/// @brief Circuit ORAM implementation.
/// @tparam T   Type of data stored in the ORAM
/// @tparam PositionType   Type of the position, default to uint64_t. Each
/// position corresponds to a path in the ORAM tree.
/// @tparam UidType   Type of the unique id, default to uint64_t. Each block in
/// the ORAM has a unique id. The unique id is used to identify the block in the
/// ORAM. Dummy blocks have a unique id of DUMMY<UidType>().
/// @tparam Z   The number of blocks in each bucket
/// @tparam stashSize   The number of blocks in the stash, not counting the root
/// bucket.
/// @tparam page_size   The size of the page if part of the data is stored in
/// external memory, default to 4096 bytes
/// @tparam check_freshness Whether to check freshness of the buckets swapped
/// from external memory
/// @tparam evict_freq  The number of reverse lexico evictions per random
/// access
/// @tparam evict_group The number of evictions to perform for each position.
/// Experiment shows that performing two evictions on one path is not much worse
/// than performing two evictions on two paths.
template <typename T, const int Z = 2, const int stashSize = 20,
          typename PositionType = uint64_t, typename UidType = uint64_t,
          const uint64_t page_size = 4096, const bool check_freshness = true,
          int evict_freq = 2, int evict_group = 2>
struct ORAM {
  using Stash = Bucket<T, stashSize, PositionType, UidType>;
  using Block_ = Block<T, PositionType, UidType>;
  using Bucket_ = Bucket<T, Z, PositionType, UidType>;

  using FreshORAMNode_ = FreshORAMNode<Bucket_>;
  using NonFreshORAMNode_ = NonFreshORAMNode<Bucket_>;

  using TreeNode_ =
      std::conditional_t<check_freshness, FreshORAMNode_, NonFreshORAMNode_>;

  using HeapTree_ = HeapTree<TreeNode_, PositionType, page_size,
                             divRoundUp(evict_freq, evict_group)>;

 private:
  HeapTree_ tree;  // underlying tree structure

  int depth = 0;           // number of levels in the tree
  PositionType _size = 0;  // size of the ORAM (number of leaves in the tree)
  PositionType capacity = 0;

  PositionType evictCounter = 0;  // counter for reverse lexicographic eviction

  std::vector<Block_> path;  // a buffer for reading and writing paths, the
  // first stashSize blocks are the stash

#ifdef DISK_IO
  uint32_t latestPrefetchReceipt = 0;
  typename HeapTree_::BatchAccessor* globalTreeAccessor = NULL;
#endif

  uint64_t rootNonce = 0;  // nonce for the root bucket

  /**
   * @brief Evict the buffered path
   *
   * @param pos The position of the path
   * @param depth The depth of the path, i.e., the number of buckets excluding
   * stash
   */
  void evictPath(PositionType pos, int depth) {
    // allocate metadata on stack to avoid contention of heap allocation
    int deepest[64];     // the deepest level an block in this level can go
    int deepestIdx[64];  // the index of the block in this level that goes
                         // deepest
    int target[64];      // the target level an block in this level should go
    bool hasEmpty[64];   // whether each level has empty slot
    std::fill(&deepest[0], &deepest[64], -1);
    std::fill(&deepestIdx[0], &deepestIdx[64], 0);
    std::fill(&target[0], &target[64], -1);
    std::fill(&hasEmpty[0], &hasEmpty[64], false);

    // first pass, root to leaf

    // pack the two integers so that we can set them using one obliMove
    struct SrcGoal {
      int src;
      int goal;
    };
    SrcGoal sg = {-1, -1};
    // first level including the stash
    for (int idx = 0; idx < stashSize + Z; ++idx) {
      int deepestLevel = CommonSuffixLength(path[idx].position, pos);
      bool deeperFlag = (!path[idx].IsDummy()) & (deepestLevel > sg.goal);
      obliMove(deeperFlag, sg.goal, deepestLevel);
      obliMove(deeperFlag, deepestIdx[0], idx);
    }
    obliMove(sg.goal >= 0, sg.src, 0);
    // the remaining levels
    for (int i = 1, idx = stashSize + Z; i < depth; ++i) {
      obliMove(sg.goal >= i, deepest[i], sg.src);
      int bucketDeepestLevel = -1;
      for (int j = 0; j < Z; ++j, ++idx) {
        int deepestLevel = CommonSuffixLength(path[idx].position, pos);
        bool isEmpty = path[idx].IsDummy();
        hasEmpty[i] |= isEmpty;
        // cache if each empty has empty slot, so in the second scan, we don't
        // need to access the buckets
        bool deeperFlag = (!isEmpty) & (deepestLevel > bucketDeepestLevel);
        obliMove(deeperFlag, bucketDeepestLevel, deepestLevel);
        obliMove(deeperFlag, deepestIdx[i], idx);
      }
      bool deeperFlag = bucketDeepestLevel > sg.goal;
      obliMove(deeperFlag, sg, {i, bucketDeepestLevel});
    }

    // second pass, leaf to root

    struct SrcDest {
      int src;
      int dest;
    };
    SrcDest sd = {-1, -1};
    for (int i = depth - 1; i > 0; --i) {
      bool isSrc = i == sd.src;
      obliMove(isSrc, target[i], sd.dest);
      obliMove(isSrc, sd, {-1, -1});
      bool changeFlag = (((sd.dest == -1) & hasEmpty[i]) | (target[i] != -1)) &
                        (deepest[i] != -1);
      obliMove(changeFlag, sd, {deepest[i], i});
    }
    obliMove(0 == sd.src, target[0], sd.dest);

    // actually moving the data

    Block_ hold;
    for (int idx = 0; idx < stashSize + Z; ++idx) {
      bool isDeepest = idx == deepestIdx[0];
      bool readAndRemoveFlag = isDeepest & (target[0] != -1);
      obliMove(readAndRemoveFlag, hold, path[idx]);
      obliMove(readAndRemoveFlag, path[idx].uid, DUMMY<UidType>());
    }
    int dest = target[0];
    for (int i = 1, idx = stashSize + Z; i < depth - 1; ++i) {
      bool hasTargetFlag = target[i] != -1;
      bool placeDummyFlag = (i == dest) & (!hasTargetFlag);
      for (int j = 0; j < Z; ++j, ++idx) {
        // case 0: level i is neither a dest and not a src
        //         hasTargetFlag = false, placeDummyFlag = false
        //         nothing will change
        // case 1: level i is a dest but not a src
        //         hasTargetFlag = false, placeDummyFlag = true
        //         hold will be swapped with each dummy slot
        //         after the first swap, hold will become dummy, and the
        //         subsequent swaps have no effect.
        // case 2: level i is a src but not a dest
        //         hasTargetFlag = true, placeDummyFlag = false
        //         hold must be dummy originally (eviction cannot carry two
        //         blocks). hold will be swapped with the slot that evicts to
        //         deepest.
        // case 3: level i is both a src and a dest
        //         hasTargetFlag = true, placeDummyFlag = false
        //         hold will be swapped with the slot that evicts to deepest,
        //         which fulfills both src and dest requirements.
        bool isDeepest = (idx == deepestIdx[i]);
        bool readAndRemoveFlag = isDeepest & hasTargetFlag;
        bool writeFlag = path[idx].IsDummy() & placeDummyFlag;
        bool swapFlag = readAndRemoveFlag | writeFlag;
        // for large element, obliSwap is faster than performing two obliMove
        obliSwap(swapFlag, hold, path[idx]);
      }
      obliMove(hasTargetFlag | placeDummyFlag, dest, target[i]);
    }

    bool placeDummyFlag = (depth - 1 == dest);
    int offset = stashSize + (depth - 1) * Z;
    bool written = false;
    for (int j = 0; j < Z; ++j) {
      int idx = offset + j;
      bool writeFlag = (path[idx].IsDummy() & placeDummyFlag) & (!written);
      written |= writeFlag;
      obliMove(writeFlag, path[idx], hold);
    }
  }

  /**
   * @brief Read the path indexed by pos, evict it, and writeback.
   *
   * @tparam _evict_group The number of evictions to perform for the path
   * @param pos The position of the path to evict
   */
  template <const int _evict_group = evict_group>
  void evict(PositionType pos) {
    PositionType nodeIdxArr[64];
    int pathDepth = readPathAndGetNodeIdxArr(pos, nodeIdxArr);
    for (int i = 0; i < evict_group; ++i) {
      evictPath(pos, pathDepth);
    }
    writeBackPath(pos, pathDepth, nodeIdxArr);
  }

  /**
   * @brief Perform multiple evictions on paths given by evict counter.
   *
   * @tparam _evict_freq The number of evictions to perform
   * @tparam _evict_group The number of evictions to perform for each position.
   * Experiment shows that performing two evictions on one path is not much
   * worse than performing two evictions on two paths.
   */
  template <const int _evict_freq = evict_freq,
            const int _evict_group = evict_group>
  void evict() {
    static constexpr int numGroup = _evict_freq / _evict_group;
    static constexpr int remaining = _evict_freq % _evict_group;
    for (int i = 0; i < numGroup; ++i) {
      evict<_evict_group>((evictCounter++) % _size);
    }
    if constexpr (remaining) {
      evict<remaining>((evictCounter++) % _size);
    }
  }

  INLINE void assertNodeFreshness(const TreeNode_& node,

                                  uint64_t expectedNonce) {
    static_assert(check_freshness);
    if (node.leftNonce + node.rightNonce != expectedNonce) {
      throw std::runtime_error("ORAM freshness check failed");
    }
  }

  /**
   * @brief Get the nonce of the child node on the path.
   *
   * @param node The parent node
   * @param level The level of the parent node
   * @param pos The position of the path
   * @return uint64_t The nonce of the child node
   */
  INLINE uint64_t getChildNonce(const TreeNode_& node, int level,
                                PositionType pos) const {
    static_assert(check_freshness);
    if ((pos >> level) & 1) {
      return node.rightNonce;
    } else {
      return node.leftNonce;
    }
  }

  /**
   * @brief Increment the nonce for the child node on the path.
   *
   * @param node The parent node
   * @param level The level of the parent node
   * @param pos The position of the path
   * @return uint64_t The nonce of the child node
   */
  INLINE void incrementNonce(TreeNode_& node, int level, PositionType pos) {
    static_assert(check_freshness);
    if ((pos >> level) & 1) {
      ++node.rightNonce;
    } else {
      ++node.leftNonce;
    }
  }

  /**
   * @brief Increment the nonce for the child node on the path and return the
   * original child nonce.
   *
   * @param node The parent node
   * @param level The level of the parent node
   * @param pos The position of the path
   * @return uint64_t The original nonce of the child node
   */
  INLINE uint64_t incrementNonceAndGetChildNonce(TreeNode_& node, int level,
                                                 PositionType pos) {
    static_assert(check_freshness);
    if ((pos >> level) & 1) {
      return node.rightNonce++;
    } else {
      return node.leftNonce++;
    }
  }

  /**
   * @brief Write a new block to the top of the tree, evict the currently
   * buffered path, and write back. If the write fails, perform some more
   * evictions and retry. The retrying leaks some information, but it happens
   * only with negligible probability.
   *
   * @tparam _evict_freq Number of evictions to perform
   * @param newBlock The new block to write
   * @param pos The position of the path
   * @param pathDepth The depth of the path
   * @param nodeIdxArr The index of the nodes in the path
   * @param retry Maximum number of retries
   */
  template <const int _evict_freq = evict_freq>
  INLINE void writeBlockWithRetry(const Block_& newBlock, PositionType pos,
                                  int pathDepth, PositionType nodeIdxArr[64],
                                  int retry = 10) {
    while (true) {
      bool success = WriteNewBlockToPath(
          path.begin(), path.begin() + stashSize + Z, newBlock);
      evictPath(pos, pathDepth);
      writeBackPath(pos, pathDepth, nodeIdxArr);
      evict<_evict_freq>();
      if (success) {
        break;
      }
      PERFCTR_INCREMENT(CIRCUITORAM_OVERFLOW);
      if (!retry) {
        // printState();
        throw std::runtime_error("ORAM overflows");
      }
      --retry;
      pos = (evictCounter++) % _size;
      pathDepth = readPathAndGetNodeIdxArr(pos, nodeIdxArr);
    }
  }

  /**
   * @brief Read a path in the ORAM tree and output the path index for later
   * write back.
   *
   * @param pos The position of the path
   * @param nodeIdxArr The buffer to output the path index
   * @return int The depth of the path
   */
  int readPathAndGetNodeIdxArr(PositionType pos, PositionType nodeIdxArr[64]) {
    uint64_t expectedNonce;

    int pathDepth;

#ifdef DISK_IO
    if (!tree.IsFullyCached()) {
      pathDepth = globalTreeAccessor->GetNodeIdxArrAndPrefetch(
          nodeIdxArr, pos, latestPrefetchReceipt);
      globalTreeAccessor->FlushRead();
      readPathFromAccessor(*globalTreeAccessor, pos, pathDepth, nodeIdxArr,
                           latestPrefetchReceipt);
      return pathDepth;
    }
#endif
    if constexpr (check_freshness) {
      expectedNonce = rootNonce;
    }
    pathDepth = tree.GetNodeIdxArr(nodeIdxArr, pos);
    for (int i = 0; i < pathDepth; ++i) {
      const TreeNode_& node = tree.GetNodeByIdx(nodeIdxArr[i]);
      if constexpr (check_freshness) {
        assertNodeFreshness(node, expectedNonce);
        expectedNonce = getChildNonce(node, i, pos);
      }
      memcpy(&path[stashSize + i * Z], node.bucket.blocks, Z * sizeof(Block_));
    }

    return pathDepth;
  }

  /**
   * @brief Write back the buffered path to the ORAM tree
   *
   * @param pos The position of the path
   */
  void writeBackPath(PositionType pos, int pathDepth,
                     const PositionType nodeIdxArr[64]) {
#ifdef DISK_IO
    if (!tree.IsFullyCached()) {
      writePathToAccessor(*globalTreeAccessor, pos, pathDepth, nodeIdxArr,
                          latestPrefetchReceipt);
      globalTreeAccessor->FlushWrite();
      // uint32_t prefetchReceipt = 0;
      // pathDepth = treeAccessor.GetNodeIdxArrAndPrefetch(nodeIdxArr,
      // pos,
      // prefetchReceipt);
      // treeAccessor.FlushRead();
      // readPathFromAccessor(treeAccessor, pos, pathDepth, nodeIdxArr,
      //                      prefetchReceipt);
      // return pathDepth;
      // return;
    }
#endif
    if constexpr (check_freshness) {
      ++rootNonce;
    }
    for (int i = 0; i < pathDepth; ++i) {
      TreeNode_& node = tree.GetNodeByIdx(nodeIdxArr[i]);
      if constexpr (check_freshness) {
        incrementNonce(node, i, pos);
      }
      memcpy(node.bucket.blocks, &path[stashSize + i * Z], Z * sizeof(Block_));
    }
  }

  void readPathFromAccessor(typename HeapTree_::BatchAccessor& treeAccessor,
                            PositionType pos, int pathDepth,
                            PositionType* nodeIdxArr,
                            uint32_t pathPrefetchReceipt) {
    uint64_t expectedNonce;
    if constexpr (check_freshness) {
      expectedNonce = rootNonce;
    }
    for (int i = 0; i < pathDepth; ++i) {
      const TreeNode_& node =
          treeAccessor.GetPrefetchedNode(nodeIdxArr[i], i, pathPrefetchReceipt);
      if constexpr (check_freshness) {
        assertNodeFreshness(node, expectedNonce);
        expectedNonce = getChildNonce(node, i, pos);
      }
      memcpy(&path[stashSize + i * Z], node.bucket.blocks, Z * sizeof(Block_));
    }
  }

  void writePathToAccessor(typename HeapTree_::BatchAccessor& treeAccessor,
                           PositionType pos, int pathDepth,
                           const PositionType* nodeIdxArr,
                           uint32_t pathPrefetchReceipt) {
    if constexpr (check_freshness) {
      ++rootNonce;
    }
    for (int i = 0; i < pathDepth; ++i) {
      TreeNode_& node =
          treeAccessor.GetPrefetchedNode(nodeIdxArr[i], i, pathPrefetchReceipt);
      if constexpr (check_freshness) {
        incrementNonce(node, i, pos);
      }
      memcpy(node.bucket.blocks, &path[stashSize + i * Z], Z * sizeof(Block_));
    }
  }

  /**
   * @brief Duplicate the positions of blocks with the same uid. Keep the
   * position of the first block with the uid.
   *
   * @param batchSize The number of blocks
   * @param newPos The new positions of the blocks. The positions of duplicate
   * blocks may be arbitrary.
   * @param uid The uids of the blocks. Requires uid to be sorted.
   * @return const std::vector<PositionType> A vector with duplicated
   * positions
   */
  const std::vector<PositionType> duplicateNewPoses(uint64_t batchSize,
                                                    const PositionType* newPos,
                                                    const UidType* uid) {
    std::vector<PositionType> dupNewPos(newPos, newPos + batchSize);
    for (int i = 1; i < batchSize; ++i) {
      obliMove(uid[i] == uid[i - 1], dupNewPos[i], dupNewPos[i - 1]);
    }
    return dupNewPos;
  }

  /**
   * @brief Change positions of duplicate blocks to random positions, to avoid
   * leaking information.
   *
   * @param batchSize The number of blocks
   * @param pos The positions of the blocks
   * @param uid The uids of the blocks. Requires uid to be sorted.
   */
  void deDuplicatePoses(uint64_t batchSize, PositionType* pos,
                        const UidType* uid) {
    for (uint64_t i = 1; i < batchSize; ++i) {
      PositionType randPos = GetRandPos();
      obliMove(uid[i] == uid[i - 1], pos[i], randPos);
    }
  }

  /**
   * @brief Duplicate the values of blocks with the same uid. Keep the value
   * of the first block with the uid.
   *
   * @param batchSize The number of blocks
   * @param out The output blocks
   * @param uid The uids of the blocks. Requires uid to be sorted.
   */
  void duplicateVal(uint64_t batchSize, T* out, const UidType* uid) {
    for (uint64_t i = 1; i < batchSize; ++i) {
      obliMove(uid[i] == uid[i - 1], out[i], out[i - 1]);
    }
  }

 public:
  ORAM() {}

  /**
   * @brief Construct a new ORAM object
   *
   * @param size Capacity of the ORAM
   */
  explicit ORAM(PositionType size) { SetSize(size); }

  /**
   * @brief Construct a new ORAM object
   *
   * @param size Capacity of the ORAM
   * @param cacheBytes The maximum size of the tree-top cache in bytes.
   */
  ORAM(PositionType size, uint64_t cacheBytes) { SetSize(size, cacheBytes); }

  /**
   * @brief Get the Stash object
   *
   * @return Stash& Returns the stash
   */
  const Stash& GetStash() { return *reinterpret_cast<Stash*>(&path[0]); }

  /**
   * @brief Set the size of the ORAM, and allocate the cache.
   *
   * @param size Capacity of the ORAM
   * @param cacheBytes The maximum size of the tree-top cache in bytes.
   */
  void SetSize(PositionType size, uint64_t cacheBytes = MAX_CACHE_SIZE) {
    if (_size) {
      throw std::runtime_error("Circuit ORAM double initialization");
    }
    capacity = size;
    _size = divRoundUp(size, Z);
    int cacheLevel = GetMaxCacheLevel<T, Z, stashSize, PositionType, UidType,
                                      check_freshness>(_size, cacheBytes);
    if (cacheLevel < 0) {
      throw std::runtime_error("Circuit ORAM cache size too small");
    }

    tree.InitWithDefault(_size, TreeNode_::DUMMY(), cacheLevel);
    depth = GetLogBaseTwo(_size - 1) + 2;
    if (depth > 64) {
      throw std::runtime_error("Circuit ORAM too large");
    }
    path.resize(stashSize + Z * depth);
#ifdef DISK_IO
    globalTreeAccessor = new typename HeapTree_::BatchAccessor(tree);
#endif
  }

  ~ORAM() {
#ifdef DISK_IO
    if (globalTreeAccessor) {
      delete globalTreeAccessor;
    }
#endif
  }

  /**
   * @brief Return the memory usage of the ORAM
   *
   * @param size Capacity of the ORAM
   * @param cacheBytes The maximum size of the tree-top cache in bytes.
   * @return uint64_t The memory usage in bytes
   */
  static uint64_t GetMemoryUsage(PositionType size,
                                 uint64_t cacheBytes = MAX_CACHE_SIZE) {
    size_t numPaths = divRoundUp(size, Z);
    return sizeof(Stash) +
           HeapTree_::GetMemoryUsage(
               numPaths,
               GetMaxCacheLevel<T, Z, stashSize, PositionType, UidType,
                                check_freshness>(numPaths, cacheBytes));
  }

  /**
   * @brief Get a random position in the ORAM
   *
   * @return PositionType the random position
   */
  INLINE PositionType GetRandPos() { return UniformRandom(_size - 1); }

  /**
   * @brief Get the an array of random positions.
   *
   * @param newPos The array to store the random positions
   * @param batchSize The number of random positions to generate
   */
  void GetRandNewPoses(PositionType* newPos, uint64_t batchSize) {
    for (int i = 0; i < batchSize; ++i) {
      newPos[i] = GetRandPos();
    }
  }

  /**
   * @brief Get the memory usage of this ORAM
   *
   * @return uint64_t The memory usage in bytes
   */
  uint64_t GetMemoryUsage() const {
    return sizeof(Stash) + tree.GetMemoryUsage();
  }

  /**
   * @brief Initialize the ORAM from a reader and output the position map.
   *
   * @tparam Reader Type of the reader
   * @tparam PosMapWriter Type of the position map writer
   * @param reader A sequential reader of the RAM data, starting from index 0.
   * @param posMapWriter A writer of the position map, the data it writes has
   * type UidBlock<PositionType, UidType>.
   */
  template <typename Reader, class PosMapWriter>
    requires Readable<Reader, T> &&
             Writable<PosMapWriter, UidBlock<PositionType, UidType>>
  void InitFromReader(Reader& reader, PosMapWriter& posMapWriter) {
    uint64_t initSize = reader.size();
    for (UidType uid = 0; uid != (UidType)initSize; ++uid) {
      PositionType newPos = Write(uid, reader.read());
      posMapWriter.write(UidBlock<PositionType, UidType>(newPos, uid));
    }
  }

  /**
   * @brief Return the capcity of the ORAM
   *
   * @return PositionType
   */
  PositionType size() const { return capacity; }

  /**
   * @brief Read a block from the ORAM and assign it to a new position
   *
   * @param pos The current position of the block
   * @param uid The unique id of the block
   * @param out The output block
   * @param newPos The new position of the block
   * @return PositionType The new position of the block
   */
  PositionType Read(PositionType pos, const UidType& uid, T& out,
                    PositionType newPos) {
    PositionType nodeIdxArr[64];
    int pathDepth = readPathAndGetNodeIdxArr(pos, nodeIdxArr);

    bool findFlag = ReadElementAndRemoveFromPath(
        path.begin(), path.begin() + (pathDepth * Z + stashSize), uid, out);

    Block_ newBlock(out, newPos, uid);
    obliMove(!findFlag, newBlock.uid, DUMMY<UidType>());
    writeBlockWithRetry(newBlock, pos, pathDepth, nodeIdxArr);
    return newPos;
  }

  /**
   * @brief Read a block from the ORAM and assign it to a random new position
   *
   * @param pos The current position of the block
   * @param uid The unique id of the block
   * @param out The output block
   * @return PositionType The random new position of the block
   */
  PositionType Read(PositionType pos, const UidType& uid, T& out) {
    PositionType newPos = GetRandPos();
    return Read(pos, uid, out, newPos);
  }

  /**
   * @brief Write a new block to the ORAM
   *
   * @tparam _evict_freq The number of evictions to perform after write
   * @param uid The unique id of the block
   * @param in The new block to write
   * @param newPos The new position of the block
   * @return PositionType The new position of the block
   */
  template <const int _evict_freq = evict_freq>
  PositionType Write(const UidType& uid, const T& in, PositionType newPos) {
    PositionType pos = (evictCounter++) % _size;
    PositionType nodeIdxArr[64];
    int pathDepth = readPathAndGetNodeIdxArr(pos, nodeIdxArr);
    Block_ newBlock(in, newPos, uid);
    writeBlockWithRetry<_evict_freq>(newBlock, pos, pathDepth, nodeIdxArr);
    return newPos;
  }

  /**
   * @brief Write a new block to the ORAM and assign it to a random new
   * position
   *
   * @tparam _evict_freq The number of evictions to perform after write
   * @param uid The unique id of the block
   * @param in The new block to write
   * @return PositionType The random new position of the block
   */
  template <const int _evict_freq = evict_freq>
  PositionType Write(const UidType& uid, const T& in) {
    PositionType newPos = GetRandPos();
    return Write<evict_freq>(uid, in, newPos);
  }

  /**
   * @brief Update a block in the ORAM and assign it to a random new position
   *
   * @tparam Func The type of the update function. The update function should
   * take a reference to the block and return a bool. If the return value is
   * true, the block is kept, otherwise it is deleted.
   * @param pos The current position of the block
   * @param uid The unique id of the block
   * @param updateFunc The update function
   * @return PositionType The random new position of the block
   */
  template <class Func>
    requires UpdateOrRemoveFunction<Func, T>
  PositionType Update(PositionType pos, const UidType& uid,
                      const Func& updateFunc) {
    T out;
    return Update(pos, uid, updateFunc, out);
  }

  /**
   * @brief Update a block in the ORAM and assign it to a new position
   *
   * @tparam Func The type of the update function. The update function should
   * take a reference to the block and return a bool. If the return value is
   * true, the block is kept, otherwise it is deleted.
   * @param pos The current position of the block
   * @param uid The unique id of the block
   * @param newPos The new position of the block
   * @param updateFunc The update function
   * @return PositionType
   */
  template <class Func>
    requires UpdateOrRemoveFunction<Func, T>
  PositionType Update(PositionType pos, const UidType& uid, PositionType newPos,
                      const Func& updateFunc) {
    T out = T();
    return Update(pos, uid, newPos, updateFunc, out);
  }

  /**
   * @brief Update a block in the ORAM and assign it to a random new position
   *
   * @tparam Func The type of the update function. The update function should
   * take a reference to the block and return a bool. If the return value is
   * true, the block is kept, otherwise it is deleted.
   * @param pos The current position of the block
   * @param uid The unique id of the block
   * @param updateFunc The update function
   * @param out Output the updated block
   * @return PositionType The random new position of the block
   */
  template <class Func>
    requires UpdateOrRemoveFunction<Func, T>
  PositionType Update(PositionType pos, const UidType& uid,
                      const Func& updateFunc, T& out) {
    return Update(pos, uid, updateFunc, out, uid);  // does not change uid
  }

  /**
   * @brief Update a block in the ORAM and assign it to a new position
   *
   * @tparam Func The type of the update function. The update function should
   * take a reference to the block and return a bool. If the return value is
   * true, the block is kept, otherwise it is deleted.
   * @param pos The current position of the block
   * @param uid The unique id of the block
   * @param newPos The new position of the block
   * @param updateFunc The update function
   * @param out Output the updated block
   * @return PositionType The new position of the block
   */
  template <class Func>
    requires UpdateOrRemoveFunction<Func, T>
  PositionType Update(PositionType pos, const UidType& uid, PositionType newPos,
                      const Func& updateFunc, T& out) {
    return Update(pos, uid, newPos, updateFunc, out,
                  uid);  // does not change uid
  }

  /**
   * @brief Update a block in the ORAM and assign it to a random new position.
   * Possibly change the uid.
   *
   * @tparam Func The type of the update function. The update function should
   * take a reference to the block and return a bool. If the return value is
   * true, the block is kept, otherwise it is deleted.
   * @param pos The current position of the block
   * @param uid The unique id of the block
   * @param updateFunc The update function
   * @param out Output the updated block
   * @param updatedUid The new uid of the block
   * @return PositionType The random new position of the block
   */
  template <class Func>
    requires UpdateOrRemoveFunction<Func, T>
  PositionType Update(PositionType pos, const UidType& uid,
                      const Func& updateFunc, T& out,
                      const UidType& updatedUid) {
    PositionType newPos = GetRandPos();
    return Update(pos, uid, newPos, updateFunc, out, updatedUid);
  }

  /**
   * @brief Update a block in the ORAM and assign it to a new position.
   * Possibly change the uid.
   *
   * @tparam Func The type of the update function. The update function should
   * take a reference to the block and return a bool. If the return value is
   * true, the block is kept, otherwise it is deleted.
   * @param pos The current position of the block
   * @param uid The unique id of the block
   * @param newPos The new position of the block
   * @param updateFunc The update function
   * @param out Output the updated block
   * @param updatedUid The new uid of the block
   * @return PositionType The random new position of the block
   */
  template <class Func>
    requires UpdateOrRemoveFunction<Func, T>
  PositionType Update(PositionType pos, const UidType& uid, PositionType newPos,
                      const Func& updateFunc, T& out,
                      const UidType& updatedUid) {
    PositionType nodeIdxArr[64];
    int pathDepth = readPathAndGetNodeIdxArr(pos, nodeIdxArr);
    ReadElementAndRemoveFromPath(
        path.begin(), path.begin() + (pathDepth * Z + stashSize), uid, out);
    bool keepFlag = updateFunc(out);
    UidType newUid = DUMMY<UidType>();
    obliMove(keepFlag, newUid, updatedUid);
    Block_ newBlock(out, newPos, newUid);
    writeBlockWithRetry(newBlock, pos, pathDepth, nodeIdxArr);
    return newPos;
  }

  /**
   * @brief Read a batch of blocks and remove them from the ORAM.
   *
   * @param batchSize The number of blocks
   * @param pos The positions of the blocks
   * @param uid The uids of the blocks
   * @param out The output blocks
   *
   */
  void BatchReadAndRemove(uint64_t batchSize, PositionType* pos,
                          const UidType* uid, T* out) {
    // mask duplicate positions
    deDuplicatePoses(batchSize, pos, uid);

// read and remove
#ifdef DISK_IO
    if (!tree.IsFullyCached()) {
      // if data is swapped to disk, minimize the number of ocalls by
      // prefetching and grouping all the writebacks together
      // freshness should only be checked when some data is out of the EPC
      typename HeapTree_::BatchAccessor treeAccessor(tree);
      const uint32_t maxPrefetchBatchSize = 16;
      PositionType nodeIdxArrs[maxPrefetchBatchSize][64];
      int pathDepths[maxPrefetchBatchSize];
      uint32_t prefetchReceipts[maxPrefetchBatchSize] = {0};
      for (uint64_t i = 0; i < batchSize;) {
        uint32_t prefetchBatchSize =
            (uint32_t)std::min((uint64_t)maxPrefetchBatchSize, batchSize - i);
        for (uint32_t j = 0; j < prefetchBatchSize; ++j) {
          pathDepths[j] = treeAccessor.GetNodeIdxArrAndPrefetch(
              &nodeIdxArrs[j][0], pos[i + j], prefetchReceipts[j]);
        }
        treeAccessor.FlushRead();

        for (uint32_t j = 0; j < prefetchBatchSize; ++j, ++i) {
          ReadElementAndRemoveFromPath(path.begin(), path.begin() + stashSize,
                                       uid[i], out[i]);
          uint64_t expectedNonce;
          if constexpr (check_freshness) {
            expectedNonce = rootNonce++;
          }
          for (int k = 0; k < pathDepths[j]; ++k) {
            TreeNode_& node = treeAccessor.GetPrefetchedNode(
                nodeIdxArrs[j][k], k, prefetchReceipts[j]);
            if constexpr (check_freshness) {
              assertNodeFreshness(node, expectedNonce);
              expectedNonce = incrementNonceAndGetChildNonce(node, k, pos[i]);
            }
            ReadElementAndRemoveFromPath(
                node.bucket.blocks, node.bucket.blocks + Z, uid[i], out[i]);
          }
        }
        treeAccessor.FlushWrite();
      }
      duplicateVal(batchSize, out, uid);
      return;
    }
#endif
    PositionType nodeIdxArr[64];
    // read and remove
    for (uint64_t i = 0; i < batchSize; ++i) {
      int pathDepth = tree.GetNodeIdxArr(&nodeIdxArr[0], pos[i]);
      ReadElementAndRemoveFromPath(path.begin(), path.begin() + stashSize,
                                   uid[i], out[i]);

      uint64_t expectedNonce;
      if constexpr (check_freshness) {
        expectedNonce = rootNonce++;
      }
      for (int j = 0; j < pathDepth; ++j) {
        TreeNode_& node = tree.GetNodeByIdx(nodeIdxArr[j]);
        if constexpr (check_freshness) {
          assertNodeFreshness(node, expectedNonce);
          expectedNonce = incrementNonceAndGetChildNonce(node, j, pos[i]);
        }
        ReadElementAndRemoveFromPath(node.bucket.blocks, node.bucket.blocks + Z,
                                     uid[i], out[i]);
      }
    }

    // propagate duplicate values
    duplicateVal(batchSize, out, uid);
  }

  /**
   * @brief Write back a batch of blocks to the ORAM. For duplicate blocks,
   * ignore all but the first block.
   *
   * @param batchSize The number of blocks
   * @param uid The uids of the blocks. Requires uid to be sorted.
   * @param newPos The new positions of the blocks
   * @param in The new blocks
   * @param writeBackFlags The write back flags. If writeBackFlags[i] is true,
   * the block is written back, otherwise it is not written back.
   */
  void BatchWriteBack(uint64_t batchSize, const UidType* uid,
                      const PositionType* newPos, const T* in,
                      const std::vector<bool>& writeBackFlags) {
#ifdef DISK_IO
    // if data is swapped to disk, minimize the number of ocalls by
    // prefetching and grouping all the writebacks together
    if (!tree.IsFullyCached()) {
      // freshness should only be checked when some data is out of the EPC
      typename HeapTree_::BatchAccessor treeAccessor(tree);
      constexpr uint32_t maxPrefetchBatchSize = 8;
      constexpr int numPathPerAccess = 1 + divRoundUp(evict_freq, evict_group);
      constexpr uint32_t maxPathCount = maxPrefetchBatchSize * numPathPerAccess;
      PositionType nodeIdxArrs[maxPathCount][64];
      int pathDepths[maxPathCount];
      uint32_t prefetchReceipts[maxPathCount] = {0};
      for (uint64_t i = 0; i < batchSize;) {
        uint64_t prefetchEvictCounter = evictCounter;
        uint32_t prefetchBatchSize =
            (uint32_t)std::min((uint64_t)maxPrefetchBatchSize, batchSize - i);

        for (uint32_t pathOffset = 0;
             pathOffset < prefetchBatchSize * numPathPerAccess; ++pathOffset) {
          PositionType p = prefetchEvictCounter++ % _size;
          pathDepths[pathOffset] = treeAccessor.GetNodeIdxArrAndPrefetch(
              nodeIdxArrs[pathOffset], p, prefetchReceipts[pathOffset]);
        }
        treeAccessor.FlushRead();
        for (uint32_t j = 0, pathOffset = 0; j < prefetchBatchSize; ++j, ++i) {
          Block_ toWrite = {in[i], newPos[i], DUMMY<UidType>()};
          bool writeBack = writeBackFlags[i];
          if (i > 0) {
            // don't write back duplicate blocks
            writeBack &= (uid[i] != uid[i - 1]);
          }
          obliMove(writeBack, toWrite.uid, uid[i]);
          bool success = false;
          for (int k = 0; k < numPathPerAccess; ++k, ++pathOffset) {
            int pathDepth = pathDepths[pathOffset];
            PositionType p = evictCounter++ % _size;
            readPathFromAccessor(treeAccessor, p, pathDepth,
                                 nodeIdxArrs[pathOffset],
                                 prefetchReceipts[pathOffset]);
            if (k == 0) {
              success = WriteNewBlockToPath(
                  path.begin(), path.begin() + stashSize + Z, toWrite);
              evictPath(p, pathDepth);
            } else {
              for (int h = 0; h < evict_group; ++h) {
                evictPath(p, pathDepth);
              }
            }

            writePathToAccessor(treeAccessor, p, pathDepth,
                                nodeIdxArrs[pathOffset],
                                prefetchReceipts[pathOffset]);
          }
          if (!success) {
            PERFCTR_INCREMENT(CIRCUITORAM_OVERFLOW);
            --i;  // retry the current element
          }
        }
        treeAccessor.FlushWrite();
      }
      return;
    }
#endif
    PositionType nodeIdxArr[64];
    for (uint64_t i = 0; i < batchSize; ++i) {
      PositionType p = evictCounter++ % _size;
      int pathDepth = readPathAndGetNodeIdxArr(p, nodeIdxArr);
      Block_ toWrite = {in[i], newPos[i], DUMMY<UidType>()};
      bool writeBack = writeBackFlags[i];
      if (i > 0) {
        // don't write back duplicate blocks
        writeBack &= (uid[i] != uid[i - 1]);
      }
      obliMove(writeBack, toWrite.uid, uid[i]);
      writeBlockWithRetry(toWrite, p, pathDepth, nodeIdxArr);
    }
  }

  /**
   * @brief Update a batch of blocks in the ORAM and assign them to new
   * positions
   *
   * @tparam Func The type of the update function
   * @param batchSize The number of blocks
   * @param pos The positions of the blocks. May be modified.
   * @param uid The uids of the blocks. Requires uid to be sorted.
   * @param newPos The new positions of the blocks
   * @param updateFunc The update function. The update function should take a
   * integer batchSize and a pointer to an array of batchSize blocks, and
   * return a std::vector<bool> of batchSize, indicating whether each block
   * should be written back.
   * @param out Pointer to an array of output blocks
   */
  template <class Func>
    requires BatchUpdateOrRemoveFunction<Func, T>
  void BatchUpdate(uint64_t batchSize, PositionType* pos, const UidType* uid,

                   const PositionType* newPos, const Func& updateFunc, T* out) {
    BatchReadAndRemove(batchSize, pos, uid, out);
    const std::vector<bool>& writeBackFlags = updateFunc(batchSize, out);
    BatchWriteBack(batchSize, uid, newPos, out, writeBackFlags);
  }

  /**
   * @brief Update a batch of blocks in the ORAM and assign them to new
   * positions
   *
   * @tparam Func The type of the update function
   * @param batchSize The number of blocks
   * @param pos The positions of the blocks. May be modified.
   * @param uid The uids of the blocks. Requires uid to be sorted.
   * @param newPos The new positions of the blocks
   * @param updateFunc The update function. The update function should take a
   * integer batchSize and a pointer to an array of batchSize blocks, and
   * return a std::vector<bool> of batchSize, indicating whether each block
   * should be written back.
   */
  template <class Func>
    requires BatchUpdateOrRemoveFunction<Func, T>
  void BatchUpdate(uint64_t batchSize, PositionType* pos, const UidType* uid,
                   const PositionType* newPos, const Func& updateFunc) {
    std::vector<T> out(batchSize);
    BatchUpdate(batchSize, pos, uid, newPos, updateFunc, &out[0]);
  }

  /**
   * @brief Update a batch of blocks in the ORAM and assign them to random new
   * positions
   *
   * @tparam Func The type of the update function
   * @param batchSize The number of blocks
   * @param pos The positions of the blocks. May be modified.
   * @param uid The uids of the blocks. Requires uid to be sorted.
   * @param updateFunc The update function. The update function should take a
   * integer batchSize and a pointer to an array of batchSize blocks, and
   * return a std::vector<bool> of batchSize, indicating whether each block
   * should be written back.
   * @return const std::vector<PositionType> The random new positions of the
   * blocks.
   */
  template <class Func>
    requires BatchUpdateOrRemoveFunction<Func, T>
  std::vector<PositionType> BatchUpdate(uint64_t batchSize, PositionType* pos,
                                        const UidType* uid,
                                        const Func& updateFunc) {
    std::vector<PositionType> newPoses(batchSize);
    GetRandNewPoses(&newPoses[0], batchSize);
    BatchUpdate(batchSize, pos, uid, &newPoses[0], updateFunc);
    return duplicateNewPoses(batchSize, &newPoses[0], uid);
  }

  void printState() {
#ifndef ENCLAVE_MODE
    std::cout << "ORAM state: " << std::endl;
    std::cout << "stash: " << std::endl;
    for (int i = 0; i < stashSize; ++i) {
      std::cout << path[i] << " ";
    }
    std::cout << std::endl;
    for (int i = 0; i < _size; ++i) {
      PositionType nodeIdxArr[64];
      int pathDepth = tree.GetNodeIdxArr(nodeIdxArr, i);
      std::cout << "\npath " << i << ": ";
      for (int j = 0; j < pathDepth; ++j) {
        const TreeNode_& node = tree.GetNodeByIdx(nodeIdxArr[j]);
        std::cout << "[";
        for (int k = 0; k < Z; ++k) {
          std::cout << node.bucket.blocks[k] << " ";
        }
        std::cout << "] ";
      }
    }
#endif
  }
};

}  // namespace ODSL::CircuitORAM