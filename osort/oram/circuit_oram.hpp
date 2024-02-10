#pragma once

#include "oram_common.hpp"

namespace ODSL::CircuitORAM {
template <typename T, const int Z = 2, const int stashSize = 50,
          typename PositionType = uint64_t, typename UidType = uint64_t,
          const uint64_t page_size = 4096>
struct ORAM {
  // using Node_ = Node<T, Z, PositionType, UidType>;
  using Stash = Bucket<T, stashSize, PositionType, UidType>;
  using Block_ = Block<T, PositionType, UidType>;
  using Bucket_ = Bucket<T, Z, PositionType, UidType>;
  using UidBlock_ = UidBlock<T, UidType>;
  using HeapTree_ = HeapTree<Bucket_, PositionType, page_size>;
  // Node_* root = nullptr;
  HeapTree_ tree;
  Stash* stash = nullptr;
  int depth = 0;
  PositionType _size = 0;

  PositionType evictCounter = 0;

  ORAM() {}

  // ORAM(const ORAM& other) = default;
  ORAM(ORAM&& other) {
    _size = other._size;
    tree = std::move(other.tree);
    stash = other.stash;
    other.stash = nullptr;
    depth = other.depth;
  }

  ORAM(PositionType size, size_t cacheBytes = 1UL << 62) {
    SetSize(size, cacheBytes);
  }

  void SetSize(PositionType size, size_t cacheBytes = 1UL << 62) {
    if (_size) {
      throw std::runtime_error("Path ORAM double initialization");
    }

    _size = size;
    int cacheLevel = GetMaxCacheLevel<T, Z, stashSize, PositionType, UidType>(
        size, cacheBytes);
    if (cacheLevel < 0) {
      throw std::runtime_error("Path ORAM cache size too small");
    }
    // printf(
    //     "tree size = %lu, cacheBytes = %lu, cacheLevel = %d, element size = "
    //     "%d\n",
    //     size, cacheBytes, cacheLevel, sizeof(T));
    Bucket_ dummyBucket;
    for (int i = 0; i < Z; ++i) {
      dummyBucket.blocks[i].uid = DUMMY<UidType>();
    }
    tree.InitWithDefault(size, dummyBucket, cacheLevel);
    depth = GetLogBaseTwo(size - 1) + 2;
    stash = new Stash();
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
    if (initSize * 5 < size()) {
      for (UidType uid = 0; uid != (UidType)initSize; ++uid) {
        PositionType newPos = Write(uid, reader.read());
        posMapWriter.write(UidBlock<PositionType, UidType>(newPos, uid));
      }
      return;
    }
    FastInitFromReader<T, Z, stashSize, PositionType, UidType>(
        reader, posMapWriter, tree);
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
    std::vector<Block_> path = ReadPath(pos);

    ReadElementAndRemoveFromPath(path, uid, out);

    Block_ newBlock(out, newPos, uid);
    WriteNewBlockToTreeTop(path, newBlock, stashSize + Z);
    EvictPath(path, pos);
    WriteBackPath(path, pos);
    Evict();
    return newPos;
  }

  void Evict(PositionType pos) {
    std::vector<Block_> path = ReadPath(pos);
    EvictPath(path, pos);
    WriteBackPath(path, pos);
  }

  void Evict() {
    Evict((evictCounter++) % size());
    Evict((evictCounter++) % size());
  }

  PositionType Read(PositionType pos, const UidType& uid, T& out) {
    PositionType newPos = UniformRandom(size() - 1);
    return Read(pos, uid, out, newPos);
  }

  PositionType Write(const UidType& uid, const T& in, PositionType newPos) {
    PositionType pos = UniformRandom(size() - 1);
    std::vector<Block_> path = ReadPath(pos);
    Block_ newBlock(in, newPos, uid);
    WriteNewBlockToTreeTop(path, newBlock, stashSize + Z);
    EvictPath(path, pos);
    WriteBackPath(path, pos);
    Evict();
    return newPos;
  }

  PositionType Write(const UidType& uid, const T& in) {
    PositionType newPos = UniformRandom(size() - 1);
    return Write(uid, in, newPos);
  }

  PositionType Update(PositionType pos, const UidType& uid,
                      std::function<void(T&)> updateFunc) {
    T out;
    return Update(pos, uid, updateFunc, out);
  }

  PositionType Update(PositionType pos, const UidType& uid, PositionType newPos,
                      std::function<void(T&)> updateFunc) {
    T out;
    return Update(pos, uid, newPos, updateFunc, out);
  }

  PositionType Update(PositionType pos, const UidType& uid,
                      std::function<void(T&)> updateFunc, T& out) {
    return Update(pos, uid, updateFunc, out, uid);  // does not change uid
  }

  PositionType Update(PositionType pos, const UidType& uid, PositionType newPos,
                      std::function<void(T&)> updateFunc, T& out) {
    return Update(pos, uid, newPos, updateFunc, out,
                  uid);  // does not change uid
  }

  PositionType Update(PositionType pos, const UidType& uid,
                      std::function<void(T&)> updateFunc, T& out,
                      const UidType& updatedUid) {
    PositionType newPos = UniformRandom(size() - 1);
    return Update(pos, uid, newPos, updateFunc, out, updatedUid);
  }

  PositionType Update(PositionType pos, const UidType& uid, PositionType newPos,
                      std::function<void(T&)> updateFunc, T& out,
                      const UidType& updatedUid) {
    std::vector<Block_> path = ReadPath(pos);

    ReadElementAndRemoveFromPath(path, uid, out);
    updateFunc(out);
    Block_ newBlock(out, newPos, updatedUid);
    WriteNewBlockToTreeTop(path, newBlock, stashSize + Z);
    EvictPath(path, pos);

    WriteBackPath(path, pos);
    Evict();
    return newPos;
  }

  std::vector<Block_> ReadPath(PositionType pos) {
    std::vector<Block_> path(stashSize + Z * depth);

    memcpy(&path[0], stash->blocks, stashSize * sizeof(Block_));

    int level = tree.ReadPath(pos, (Bucket_*)(&path[stashSize]));
    path.resize(stashSize + Z * level);
#ifndef NDEBUG
    for (int i = stashSize; i < path.size(); ++i) {
      int lv = (i - stashSize) / Z;
      size_t mask = (1UL << lv) - 1;
      if (!(path[i].uid == DUMMY<UidType>() ||
            (path[i].position & mask) == (pos & mask))) {
        printf(
            "position mismatch pos = %lu at level %d, path[%d].position = %lu, "
            "mask = "
            "%lu\n",
            pos, lv, i, path[i].position, mask);
      }
    }
#endif
    return path;
  }

  void WriteBackPath(const std::vector<Block_>& path, PositionType pos) {
    memcpy(stash->blocks, &path[0], stashSize * sizeof(Block_));
    tree.WritePath(pos, (Bucket_*)(&path[stashSize]));
  }

  static void EvictPath(std::vector<Block_>& path, PositionType pos) {
    int depth = (path.size() - stashSize) / Z;
    std::vector<int> deepest(depth, -1);
    std::vector<int> deepestIdx(depth, 0);
    int src = -1;
    int goal = -1;

    for (int i = 0; i < depth; ++i) {
      obliMove(goal >= i, deepest[i], src);
      int bucketDeepestLevel = -1;
      for (int j = (i == 0 ? -stashSize : 0); j < Z; ++j) {
        int idx = stashSize + i * Z + j;
        int deepestLevel = commonSuffixLength(path[idx].position, pos);
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
    std::vector<int> target(depth, -1);
    for (int i = depth - 1; i > 0; --i) {
      obliMove(i == src, target[i], dest);
      obliMove(i == src, dest, -1);
      obliMove(i == src, src, -1);
      bool hasEmpty = false;
      for (int j = 0; j < Z; ++j) {
        int idx = stashSize + i * Z + j;
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
      for (int j = (i == 0 ? -stashSize : 0); j < Z; ++j) {
        int idx = stashSize + i * Z + j;
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
        int idx = stashSize + i * Z + j;
        bool writeFlag = needWriteFlag & path[idx].isDummy();
        obliMove(writeFlag, path[idx], toWrite);
        needWriteFlag &= !writeFlag;
      }
    }
    // for (int i = 0; i < path.size(); ++i) {
    //   printf("%ld ", path[i].isDummy() ? -1UL : (int64_t)path[i].position);
    //   if (i >= stashSize + Z && (path.size() - i) % Z == 1) {
    //     printf("\n");
    //   }
    // }
  }
};
}  // namespace ODSL::CircuitORAM