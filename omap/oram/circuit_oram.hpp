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
/// @tparam evict_freq  The number of reverse lexico evictions per random access
template <typename T, const int Z = 2, const int stashSize = 20,
          typename PositionType = uint64_t, typename UidType = uint64_t,
          const uint64_t page_size = 4096, int evict_freq = 2>
struct ORAM {
  using Stash = Bucket<T, stashSize, PositionType, UidType>;
  using Block_ = Block<T, PositionType, UidType>;
  using Bucket_ = Bucket<T, Z, PositionType, UidType>;
  using HeapTree_ = HeapTree<Bucket_, PositionType, page_size, evict_freq>;

 private:
  HeapTree_ tree;  // underlying tree structure

  Stash* stash = nullptr;  // stash for Circuit ORAM
  int depth = 0;           // number of levels in the tree
  PositionType _size = 0;  // size of the ORAM (number of leaves in the tree)

  PositionType evictCounter = 0;  // counter for reverse lexicographic eviction

  std::vector<Block_> path;  // a temporary buffer for reading and writing paths

  /**
   * @brief Evict a path of pathLen blocks
   *
   * @param path The path to evict
   * @param pos The position of the path
   * @param pathLen Only evict the first pathLen blocks in the path
   */
  static void evictPath(std::vector<Block_>& path, PositionType pos,
                        int pathLen) {
    int depth = (pathLen - stashSize) / Z;
    return evictPath(path, pos, depth, stashSize);
  }

  /**
   * @brief Evict an entire path
   *
   * @param path The path to evict
   * @param pos The position of the path
   */
  static void evictPath(std::vector<Block_>& path, PositionType pos) {
    int depth = (path.size() - stashSize) / Z;
    return evictPath(path, pos, depth, stashSize);
  }

  /**
   * @brief Evict path
   *
   * @param path The path to evict. The first actualStashSize blocks are stash.
   * @param pos The position of the path
   * @param depth The depth of the path, i.e., the number of buckets excluding
   * stash
   * @param actualStashSize The actual size of the stash in the path
   */
  static void evictPath(std::vector<Block_>& path, PositionType pos, int depth,
                        int actualStashSize) {
    // allocate metadata on stack to avoid contention of heap allocation
    int deepest[64];
    int deepestIdx[64];
    int target[64];
    std::fill(&deepest[0], &deepest[depth], -1);
    std::fill(&deepestIdx[0], &deepestIdx[depth], 0);
    std::fill(&target[0], &target[depth], -1);
    int src = -1;
    int goal = -1;

    for (int i = 0; i < depth; ++i) {
      obliMove(goal >= i, deepest[i], src);
      int bucketDeepestLevel = -1;
      for (int j = (i == 0 ? -actualStashSize : 0); j < Z; ++j) {
        int idx = actualStashSize + i * Z + j;
        int deepestLevel = commonSuffixLength(path[idx].position, pos);
        bool deeperFlag =
            !path[idx].IsDummy() & (deepestLevel > bucketDeepestLevel);
        obliMove(deeperFlag, bucketDeepestLevel, deepestLevel);
        obliMove(deeperFlag, deepestIdx[i], idx);
      }
      bool deeperFlag = bucketDeepestLevel > goal;
      obliMove(deeperFlag, goal, bucketDeepestLevel);
      obliMove(deeperFlag, src, i);
    }

    int dest = -1;
    src = -1;
    for (int i = depth - 1; i > 0; --i) {
      obliMove(i == src, target[i], dest);
      obliMove(i == src, dest, -1);
      obliMove(i == src, src, -1);
      bool hasEmpty = false;
      for (int j = 0; j < Z; ++j) {
        int idx = actualStashSize + i * Z + j;
        hasEmpty |= path[idx].IsDummy();
      }
      bool changeFlag =
          (((dest == -1) & hasEmpty) | (target[i] != -1)) & (deepest[i] != -1);
      obliMove(changeFlag, src, deepest[i]);
      obliMove(changeFlag, dest, i);
    }
    obliMove(0 == src, target[0], dest);

    Block_ hold;
    hold.uid = DUMMY<UidType>();
    dest = -1;

    for (int i = 0; i < depth; ++i) {
      Block_ toWrite = hold;
      bool placeFlag = !hold.IsDummy() & (i == dest);
      obliMove(!placeFlag, toWrite.uid, DUMMY<UidType>());
      obliMove(placeFlag, hold.uid, DUMMY<UidType>());
      obliMove(placeFlag, dest, -1);
      bool hasTargetFlag = target[i] != -1;
      for (int j = (i == 0 ? -actualStashSize : 0); j < Z; ++j) {
        int idx = actualStashSize + i * Z + j;
        bool isDeepest = (idx == deepestIdx[i]);
        bool readAndRemoveFlag = isDeepest & hasTargetFlag;
        obliMove(readAndRemoveFlag, hold, path[idx]);
        obliMove(readAndRemoveFlag, path[idx].uid, DUMMY<UidType>());
      }
      obliMove(hasTargetFlag, dest, target[i]);
      if (i == 0) {
        continue;  // no need to write to stash
      }
      bool needWriteFlag = !toWrite.IsDummy();
      for (int j = 0; j < Z; ++j) {
        int idx = actualStashSize + i * Z + j;
        bool writeFlag = needWriteFlag & path[idx].IsDummy();
        obliMove(writeFlag, path[idx], toWrite);
        needWriteFlag &= !writeFlag;
      }
    }
  }

  /**
   * @brief Read the path indexed by pos, evict it, and writeback.
   *
   * @param pos The position of the path to evict
   */
  void evict(PositionType pos) {
    int len = readPath(pos, path);
    evictPath(path, pos, len);
    writeBackPath(path, pos);
  }

  /**
   * @brief Perform multiple evictions on paths given by evict counter.
   *
   * @tparam _evict_freq The number of evictions to perform
   */
  template <const int _evict_freq = evict_freq>
  void evict() {
    for (int i = 0; i < _evict_freq; ++i) {
      evict((evictCounter++) % size());
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
   * @param pathLen The length of the path
   * @param retry Maximum number of retries
   */
  template <const int _evict_freq = evict_freq>
  INLINE void writeBlockWithRetry(const Block_& newBlock, PositionType pos,
                                  int pathLen, int retry = 10) {
    while (true) {
      bool success = WriteNewBlockToTreeTop(path, newBlock, stashSize + Z);
      evictPath(path, pos, pathLen);
      writeBackPath(path, pos);
      evict<_evict_freq>();
      if (success) {
        break;
      }
      if (!retry) {
        throw std::runtime_error("ORAM read failed");
      }
      --retry;
      pos = (evictCounter++) % size();
      pathLen = readPath(pos, path);
    }
  }

  /**
   * @brief Read a path in the ORAM
   *
   * @param pos The position of the path
   * @param path The buffer to store the path
   * @return int The length of the path
   */
  int readPath(PositionType pos, std::vector<Block_>& path) {
    Assert(path.size() == stashSize + Z * depth);

    memcpy(&path[0], stash->blocks, stashSize * sizeof(Block_));

    int level = tree.ReadPath(pos, (Bucket_*)(&path[stashSize]));
    return stashSize + Z * level;
  }

  /**
   * @brief Write back a path to the ORAM
   *
   * @param path The path to write back
   * @param pos The position of the path
   */
  void writeBackPath(const std::vector<Block_>& path, PositionType pos) {
    memcpy(stash->blocks, &path[0], stashSize * sizeof(Block_));
    tree.WritePath(pos, (const Bucket_*)(&path[stashSize]));
  }

  /**
   * @brief Get a random position in the ORAM
   *
   * @return PositionType the random position
   */
  INLINE PositionType getRandPos() {
    return (PositionType)UniformRandom(size() - 1);
  }

  /**
   * @brief Get a vector of random positions.
   *
   * @param batchSize The number of random positions to generate
   * @return const std::vector<PositionType> The vector of random positions
   */
  const std::vector<PositionType> getRandNewPoses(uint64_t batchSize) {
    std::vector<PositionType> newPos(batchSize);
    for (int i = 0; i < batchSize; ++i) {
      newPos[i] = getRandPos();
    }
    return newPos;
  }

  /**
   * @brief Duplicate the positions of blocks with the same uid. Keep the
   * position of the first block with the uid.
   *
   * @param batchSize The number of blocks
   * @param newPos The new positions of the blocks. The positions of duplicate
   * blocks may be arbitrary.
   * @param uid The uids of the blocks. Requires uid to be sorted.
   * @return const std::vector<PositionType> A vector with duplicated positions
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
      PositionType randPos = getRandPos();
      obliMove(uid[i] == uid[i - 1], pos[i], randPos);
    }
  }

  /**
   * @brief Duplicate the values of blocks with the same uid. Keep the value of
   * the first block with the uid.
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
  ORAM(PositionType size) { SetSize(size); }

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
  const Stash& GetStash() { return *stash; }

  /**
   * @brief Set the size of the ORAM, and allocate the cache.
   *
   * @param size Capacity of the ORAM
   * @param cacheBytes The maximum size of the tree-top cache in bytes.
   */
  void SetSize(PositionType size, uint64_t cacheBytes = 1UL << 62) {
    if (_size) {
      throw std::runtime_error("Circuit ORAM double initialization");
    }

    _size = size;
    int cacheLevel = GetMaxCacheLevel<T, Z, stashSize, PositionType, UidType>(
        size, cacheBytes);
    if (cacheLevel < 0) {
      throw std::runtime_error("Circuit ORAM cache size too small");
    }
    Bucket_ dummyBucket;
    for (int i = 0; i < Z; ++i) {
      dummyBucket.blocks[i].uid = DUMMY<UidType>();
    }
    tree.InitWithDefault(size, dummyBucket, cacheLevel);
    depth = (int)GetLogBaseTwo(size - 1) + 2;
    if (depth > 64) {
      throw std::runtime_error("Circuit ORAM too large");
    }
    stash = new Stash();
    path.resize(stashSize + Z * depth);
  }

  /**
   * @brief Return the memory usage of the ORAM
   *
   * @param size Capacity of the ORAM
   * @param cacheBytes The maximum size of the tree-top cache in bytes.
   * @return uint64_t The memory usage in bytes
   */
  static uint64_t GetMemoryUsage(PositionType size,
                                 uint64_t cacheBytes = 1UL << 62) {
    return sizeof(Stash) +
           HeapTree_::getMemoryUsage(
               size, GetMaxCacheLevel<T, Z, stashSize, PositionType, UidType>(
                         size, cacheBytes));
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
   * @tparam Reader
   * @tparam PosMapWriter
   * @param reader A sequential reader of the RAM data, starting from index 0.
   * @param posMapWriter A writer of the position map, the data it writes has
   * type UidBlock<PositionType, UidType>.
   */
  template <typename Reader, class PosMapWriter>
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
  PositionType size() const { return _size; }

  ~ORAM() {
    if (stash) {
      delete stash;
      stash = nullptr;
    }
  }

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
    int len = readPath(pos, path);

    ReadElementAndRemoveFromPath(path.begin(), path.begin() + len, uid, out);

    Block_ newBlock(out, newPos, uid);
    writeBlockWithRetry(newBlock, pos, len);
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
    PositionType newPos = getRandPos();
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
    PositionType pos = (evictCounter++) % size();
    int len = readPath(pos, path);
    Block_ newBlock(in, newPos, uid);
    writeBlockWithRetry<_evict_freq>(newBlock, pos, len);
    return newPos;
  }

  /**
   * @brief Write a new block to the ORAM and assign it to a random new position
   *
   * @tparam _evict_freq The number of evictions to perform after write
   * @param uid The unique id of the block
   * @param in The new block to write
   * @return PositionType The random new position of the block
   */
  template <const int _evict_freq = evict_freq>
  PositionType Write(const UidType& uid, const T& in) {
    PositionType newPos = getRandPos();
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
  PositionType Update(PositionType pos, const UidType& uid,
                      const Func& updateFunc, T& out,
                      const UidType& updatedUid) {
    PositionType newPos = getRandPos();
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
  PositionType Update(PositionType pos, const UidType& uid, PositionType newPos,
                      const Func& updateFunc, T& out,
                      const UidType& updatedUid) {
    int len = readPath(pos, path);
    ReadElementAndRemoveFromPath(path.begin(), path.begin() + len, uid, out);
    bool keepFlag = updateFunc(out);
    UidType newUid = DUMMY<UidType>();
    obliMove(keepFlag, newUid, updatedUid);
    Block_ newBlock(out, newPos, updatedUid);
    writeBlockWithRetry(newBlock, pos, len);
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

    // only copy stash once
    memcpy(&path[0], stash->blocks, stashSize * sizeof(Block_));

    // read and remove
    for (uint64_t i = 0; i < batchSize; ++i) {
      int actualLevel = tree.ReadPath(pos[i], (Bucket_*)&(path[stashSize]));
      ReadElementAndRemoveFromPath(path.begin(),
                                   path.begin() + stashSize + Z * actualLevel,
                                   uid[i], out[i]);
      tree.WritePath(pos[i], (Bucket_*)&(path[stashSize]));
    }

    // copy stash back
    memcpy(stash->blocks, &path[0], stashSize * sizeof(Block_));

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
    // only copy stash once
    memcpy(&path[0], stash->blocks, stashSize * sizeof(Block_));
    for (uint64_t i = 0; i < batchSize; ++i) {
      PositionType p = evictCounter++ % size();
      int actualLevel = tree.ReadPath(p, (Bucket_*)&(path[stashSize]));
      Block_ toWrite = {in[i], newPos[i], DUMMY<UidType>()};
      bool writeBack = writeBackFlags[i];
      if (i > 0) {
        // don't write back duplicate blocks
        writeBack &= (uid[i] != uid[i - 1]);
      }
      obliMove(writeBack, toWrite.uid, uid[i]);
      writeBlockWithRetry(toWrite, p, stashSize + Z * actualLevel);
    }

    memcpy(stash->blocks, &path[0], stashSize * sizeof(Block_));
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
   * integer batchSize and a pointer to an array of batchSize blocks, and return
   * a std::vector<bool> of batchSize, indicating whether each block should be
   * written back.
   * @param out Pointer to an array of output blocks
   */
  template <class Func>
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
   * integer batchSize and a pointer to an array of batchSize blocks, and return
   * a std::vector<bool> of batchSize, indicating whether each block should be
   * written back.
   */
  template <class Func>
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
   * integer batchSize and a pointer to an array of batchSize blocks, and return
   * a std::vector<bool> of batchSize, indicating whether each block should be
   * written back.
   * @return const std::vector<PositionType> The random new positions of the
   * blocks.
   */
  template <class Func>
  std::vector<PositionType> BatchUpdate(uint64_t batchSize, PositionType* pos,
                                        const UidType* uid,
                                        const Func& updateFunc) {
    const std::vector<PositionType>& newPoses = getRandNewPoses(batchSize);
    BatchUpdate(batchSize, pos, uid, &newPoses[0], updateFunc);
    return duplicateNewPoses(batchSize, &newPoses[0], uid);
  }
};

}  // namespace ODSL::CircuitORAM