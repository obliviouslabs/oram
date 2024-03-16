#pragma once

#include "oram_common.hpp"

/// @brief This file implements Circuit ORAM

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
  HeapTree_ tree;

  Stash* stash = nullptr;
  int depth = 0;
  PositionType _size = 0;

  PositionType evictCounter = 0;

  std::vector<Block_> path;

  static void EvictPath(std::vector<Block_>& path, PositionType pos,
                        int pathLen) {
    int depth = (pathLen - stashSize) / Z;
    return EvictPath(path, pos, depth, stashSize, 0);
  }

  static void EvictPath(std::vector<Block_>& path, PositionType pos) {
    int depth = (path.size() - stashSize) / Z;
    return EvictPath(path, pos, depth, stashSize, 0);
  }

  static void EvictPath(std::vector<Block_>& path, PositionType pos, int depth,
                        int actualStashSize, int k) {
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
        int deepestLevel = commonSuffixLength(path[idx].position, pos) - k;
        bool deeperFlag =
            !path[idx].isDummy() & (deepestLevel > bucketDeepestLevel);
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
        hasEmpty |= path[idx].isDummy();
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
      bool placeFlag = !hold.isDummy() & (i == dest);
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
      bool needWriteFlag = !toWrite.isDummy();
      for (int j = 0; j < Z; ++j) {
        int idx = actualStashSize + i * Z + j;
        bool writeFlag = needWriteFlag & path[idx].isDummy();
        obliMove(writeFlag, path[idx], toWrite);
        needWriteFlag &= !writeFlag;
      }
    }
  }

  void Evict(PositionType pos) {
    int len = ReadPath(pos, path);
    EvictPath(path, pos, len);
    WriteBackPath(path, pos);
  }

  template <const int _evict_freq = evict_freq>
  void Evict() {
    for (int i = 0; i < _evict_freq; ++i) {
      Evict((evictCounter++) % size());
    }
  }

  template <const int _evict_freq = evict_freq>
  INLINE void writeBlockWithRetry(const Block_& newBlock, PositionType pos,
                                  int pathLen, int retry = 10) {
    while (true) {
      bool success = WriteNewBlockToTreeTop(path, newBlock, stashSize + Z);
      EvictPath(path, pos, pathLen);
      WriteBackPath(path, pos);
      Evict<_evict_freq>();
      if (success) {
        break;
      }
      if (!retry) {
        throw std::runtime_error("ORAM read failed");
      }
      --retry;
      pos = (evictCounter++) % size();
      pathLen = ReadPath(pos, path);
    }
  }

 public:
  ORAM() {}

  ORAM(PositionType size) { SetSize(size); }

  ORAM(PositionType size, size_t cacheBytes) { SetSize(size, cacheBytes); }

  Stash* getStash() { return stash; }

  void SetSize(PositionType size, size_t cacheBytes = 1UL << 62) {
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
    depth = GetLogBaseTwo(size - 1) + 2;
    if (depth > 64) {
      throw std::runtime_error("Circuit ORAM too large");
    }
    stash = new Stash();
    path.resize(stashSize + Z * depth);
  }

  static size_t GetMemoryUsage(PositionType size,
                               size_t cacheBytes = 1UL << 62) {
    return sizeof(Stash) +
           HeapTree_::getMemoryUsage(
               size, GetMaxCacheLevel<T, Z, stashSize, PositionType, UidType>(
                         size, cacheBytes));
  }

  size_t GetMemoryUsage() const {
    return sizeof(Stash) + tree.GetMemoryUsage();
  }

  template <typename Reader, class PosMapWriter>
  void InitFromReader(Reader& reader, PosMapWriter& posMapWriter) {
    size_t initSize = reader.size();
    for (UidType uid = 0; uid != (UidType)initSize; ++uid) {
      PositionType newPos = Write(uid, reader.read());
      posMapWriter.write(UidBlock<PositionType, UidType>(newPos, uid));
    }
  }

  PositionType size() const { return _size; }

  ~ORAM() {
    if (stash) {
      delete stash;
      stash = nullptr;
    }
  }

  ORAM& operator=(const ORAM& other) = default;
  ORAM& operator=(ORAM&& other) = default;

  PositionType Read(PositionType pos, const UidType& uid, T& out,
                    PositionType newPos) {
    int len = ReadPath(pos, path);

    ReadElementAndRemoveFromPath(path.begin(), path.begin() + len, uid, out);

    Block_ newBlock(out, newPos, uid);
    writeBlockWithRetry(newBlock, pos, len);
    return newPos;
  }

  PositionType Read(PositionType pos, const UidType& uid, T& out) {
    PositionType newPos = UniformRandom(size() - 1);
    return Read(pos, uid, out, newPos);
  }

  template <const int _evict_freq = evict_freq>
  PositionType Write(const UidType& uid, const T& in, PositionType newPos) {
    PositionType pos = (evictCounter++) % size();
    int len = ReadPath(pos, path);
    Block_ newBlock(in, newPos, uid);
    writeBlockWithRetry<_evict_freq>(newBlock, pos, len);
    return newPos;
  }

  template <const int _evict_freq = evict_freq>
  PositionType Write(const UidType& uid, const T& in) {
    PositionType newPos = UniformRandom(size() - 1);
    return Write<evict_freq>(uid, in, newPos);
  }

  template <class Func>
  PositionType Update(PositionType pos, const UidType& uid,
                      const Func& updateFunc) {
    T out;
    return Update(pos, uid, updateFunc, out);
  }

  template <class Func>
  PositionType Update(PositionType pos, const UidType& uid, PositionType newPos,
                      const Func& updateFunc) {
    T out;
    return Update(pos, uid, newPos, updateFunc, out);
  }

  template <class Func>
  PositionType Update(PositionType pos, const UidType& uid,
                      const Func& updateFunc, T& out) {
    return Update(pos, uid, updateFunc, out, uid);  // does not change uid
  }

  template <class Func>
  PositionType Update(PositionType pos, const UidType& uid, PositionType newPos,
                      const Func& updateFunc, T& out) {
    return Update(pos, uid, newPos, updateFunc, out,
                  uid);  // does not change uid
  }

  template <class Func>
  PositionType Update(PositionType pos, const UidType& uid,
                      const Func& updateFunc, T& out,
                      const UidType& updatedUid) {
    PositionType newPos = UniformRandom(size() - 1);
    return Update(pos, uid, newPos, updateFunc, out, updatedUid);
  }

  template <class Func>
  PositionType Update(PositionType pos, const UidType& uid, PositionType newPos,
                      const Func& updateFunc, T& out,
                      const UidType& updatedUid) {
    int len = ReadPath(pos, path);
    ReadElementAndRemoveFromPath(path.begin(), path.begin() + len, uid, out);
    bool keepFlag = updateFunc(out);
    UidType newUid = DUMMY<UidType>();
    obliMove(keepFlag, newUid, updatedUid);
    Block_ newBlock(out, newPos, updatedUid);
    writeBlockWithRetry(newBlock, pos, len);
    return newPos;
  }

  int ReadPath(PositionType pos, std::vector<Block_>& path) {
    Assert(path.size() == stashSize + Z * depth);

    memcpy(&path[0], stash->blocks, stashSize * sizeof(Block_));

    int level = tree.ReadPath(pos, (Bucket_*)(&path[stashSize]));
    return stashSize + Z * level;
  }

  void WriteBackPath(const std::vector<Block_>& path, PositionType pos) {
    memcpy(stash->blocks, &path[0], stashSize * sizeof(Block_));
    tree.WritePath(pos, (Bucket_*)(&path[stashSize]));
  }

  // requires uid to be sorted
  const std::vector<PositionType> GetRandNewPoses(uint64_t batchSize) {
    std::vector<PositionType> newPos(batchSize);
    for (int i = 0; i < batchSize; ++i) {
      newPos[i] = UniformRandom(size() - 1);
    }
    return newPos;
  }

  const std::vector<PositionType> duplicateNewPoses(uint64_t batchSize,
                                                    const PositionType* newPos,
                                                    const UidType* uid) {
    std::vector<PositionType> dupNewPos(newPos, newPos + batchSize);
    for (int i = 1; i < batchSize; ++i) {
      obliMove(uid[i] == uid[i - 1], dupNewPos[i], dupNewPos[i - 1]);
    }
    return dupNewPos;
  }

  // requires uid to be sorted
  void deDuplicatePoses(uint64_t batchSize, PositionType* pos,
                        const UidType* uid) {
    for (int i = 1; i < batchSize; ++i) {
      PositionType randPos = UniformRandom(size() - 1);
      obliMove(uid[i] == uid[i - 1], pos[i], randPos);
    }
  }

  void duplicateVal(uint64_t batchSize, T* out, const UidType* uid) {
    for (int i = 1; i < batchSize; ++i) {
      obliMove(uid[i] == uid[i - 1], out[i], out[i - 1]);
    }
  }

  // requries uid to be sorted
  void BatchReadAndRemove(uint64_t batchSize, PositionType* pos,
                          const UidType* uid, T* out) {
    deDuplicatePoses(batchSize, pos, uid);

    // reads elements in the subtree
    std::vector<Block_> localPath(stashSize + Z * depth);
    memcpy(&localPath[0], stash->blocks, stashSize * sizeof(Block_));
    for (int i = 0; i < batchSize; ++i) {
      int actualLevel =
          tree.ReadPath(pos[i], (Bucket_*)&(localPath[stashSize]));
      ReadElementAndRemoveFromPath(
          localPath.begin(), localPath.begin() + stashSize + Z * actualLevel,
          uid[i], out[i]);
      tree.WritePath(pos[i], (Bucket_*)&(localPath[stashSize]));
    }
    memcpy(stash->blocks, &localPath[0], stashSize * sizeof(Block_));
    duplicateVal(batchSize, out, uid);
  }

  void BatchWriteBack(uint64_t batchSize, const UidType* uid,
                      const PositionType* newPos, const T* in,
                      const std::vector<bool>& writeBackFlags) {
    std::vector<Block_> localPath(stashSize + Z * depth);
    std::vector<Block_> toWrite(batchSize);
    for (uint64_t i = 0; i < batchSize; ++i) {
      toWrite[i].uid = DUMMY<UidType>();
      bool writeBack = writeBackFlags[i];
      if (i > 0) {
        writeBack &= (uid[i] != uid[i - 1]);
      }
      obliMove(writeBack, toWrite[i].uid, uid[i]);
      // for duplicate uids the positions should be random
      toWrite[i].position = newPos[i];
      toWrite[i].data = in[i];
    }

    // std::vector<Block_> localPath(stashSize + Z * depth);
    memcpy(&localPath[0], stash->blocks, stashSize * sizeof(Block_));
    for (uint64_t i = 0; i < batchSize; ++i) {
      const Block_& toWriteBlock = toWrite[i];

      // this part may be slow (5% run time for batch size = 1000)
      bool success = WriteNewBlockToTreeTop(localPath, toWriteBlock, stashSize);
      Assert(success || toWriteBlock.isDummy());
      for (int i = 0; i < evict_freq + 1; ++i) {
        PositionType p = evictCounter++ % size();
        int actualLevel = tree.ReadPath(p, (Bucket_*)&(localPath[stashSize]));
        EvictPath(localPath, p, actualLevel, stashSize, 0);
        tree.WritePath(p, (Bucket_*)&(localPath[stashSize]));
      }
    }

    memcpy(stash->blocks, &localPath[0], stashSize * sizeof(Block_));
  }

  template <class Func>
  void BatchUpdate(uint64_t batchSize, PositionType* pos, const UidType* uid,

                   const PositionType* newPos, const Func& updateFunc, T* out) {
    BatchReadAndRemove(batchSize, pos, uid, out);
    const std::vector<bool>& writeBackFlags = updateFunc(batchSize, out);
    BatchWriteBack(batchSize, uid, newPos, out, writeBackFlags);
  }

  template <class Func>
  void BatchUpdate(uint64_t batchSize, PositionType* pos, const UidType* uid,
                   const PositionType* newPos, const Func& updateFunc) {
    std::vector<T> out(batchSize);
    BatchUpdate(batchSize, pos, uid, newPos, updateFunc, &out[0]);
  }

  template <class Func>
  std::vector<PositionType> BatchUpdate(uint64_t batchSize, PositionType* pos,
                                        const UidType* uid,
                                        const Func& updateFunc) {
    const std::vector<PositionType>& newPoses = GetRandNewPoses(batchSize);
    BatchUpdate(batchSize, pos, uid, &newPoses[0], updateFunc);
    return duplicateNewPoses(batchSize, &newPoses[0], uid);
  }
};

}  // namespace ODSL::CircuitORAM