#pragma once

#include "omp.h"
#include "oram_common.hpp"
namespace ODSL::CircuitORAM {
template <typename T, const int Z = 2, const int stashSize = 20,
          typename PositionType = uint64_t, typename UidType = uint64_t,
          const uint64_t page_size = 4096, int evict_freq = 2>
struct ORAM {
  // using Node_ = Node<T, Z, PositionType, UidType>;
  using Stash = Bucket<T, stashSize, PositionType, UidType>;
  using Block_ = Block<T, PositionType, UidType>;
  using Bucket_ = Bucket<T, Z, PositionType, UidType>;
  using UidBlock_ = UidBlock<T, UidType>;
  using HeapTree_ = HeapTree<Bucket_, PositionType, page_size, evict_freq>;
  // Node_* root = nullptr;
  HeapTree_ tree;
  Stash* stash = nullptr;
  int depth = 0;
  PositionType _size = 0;

  PositionType evictCounter = 0;

  std::vector<Block_> path;

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
      throw std::runtime_error("Circuit ORAM double initialization");
    }

    _size = size;
    int cacheLevel = GetMaxCacheLevel<T, Z, stashSize, PositionType, UidType>(
        size, cacheBytes);
    if (cacheLevel < 0) {
      throw std::runtime_error("Circuit ORAM cache size too small");
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
    int len = ReadPath(pos, path);

    ReadElementAndRemoveFromPath(path.begin(), path.begin() + len, uid, out);

    Block_ newBlock(out, newPos, uid);
    int retry = 10;
    while (true) {
      bool success = WriteNewBlockToTreeTop(path, newBlock, stashSize + Z);
      EvictPath(path, pos, len);
      WriteBackPath(path, pos);
      Evict();
      if (success) {
        break;
      }
      if (!retry) {
        throw std::runtime_error("ORAM read failed");
      }
      --retry;
      pos = (evictCounter++) % size();
      len = ReadPath(pos, path);
    }
    return newPos;
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

  PositionType Read(PositionType pos, const UidType& uid, T& out) {
    PositionType newPos = UniformRandom(size() - 1);
    return Read(pos, uid, out, newPos);
  }

  template <const int _evict_freq = evict_freq>
  PositionType Write(const UidType& uid, const T& in, PositionType newPos) {
    PositionType pos = (evictCounter++) % size();
    int len = ReadPath(pos, path);
    Block_ newBlock(in, newPos, uid);
    int retry = 10;
    while (true) {
      bool success = WriteNewBlockToTreeTop(path, newBlock, stashSize + Z);
      EvictPath(path, pos, len);
      WriteBackPath(path, pos);
      Evict<_evict_freq>();
      if (success) {
        break;
      }
      if (!retry) {
        throw std::runtime_error("ORAM write failed");
      }
      --retry;
      pos = (evictCounter++) % size();
      len = ReadPath(pos, path);
    }
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

    ReadElementAndRemoveFromPath(path, uid, out);
    bool keepFlag = updateFunc(out);
    UidType newUid = DUMMY<UidType>();
    obliMove(keepFlag, newUid, updatedUid);
    Block_ newBlock(out, newPos, updatedUid);
    int retry = 10;
    while (true) {
      bool success = WriteNewBlockToTreeTop(path, newBlock, stashSize + Z);
      success |= !keepFlag;
      EvictPath(path, pos, len);
      WriteBackPath(path, pos);
      Evict();
      if (success) {
        break;
      }
      if (!retry) {
        throw std::runtime_error("ORAM update failed");
      }
      --retry;
      pos = (evictCounter++) % size();
      len = ReadPath(pos, path);
    }
    return newPos;
  }

  std::vector<PositionType> BatchUpdate(
      const std::vector<PositionType>& pos, const std::vector<UidType>& uid,
      std::function<std::vector<bool>(std::vector<T>&)> updateFunc) {
    std::vector<PositionType> newPos(pos.size());
    for (int i = 0; i < pos.size(); ++i) {
      newPos[i] = UniformRandom(size() - 1);
    }
    BatchUpdate(pos, uid, newPos, updateFunc);
    return newPos;
  }

  void BatchUpdate(
      const std::vector<PositionType>& pos, const std::vector<UidType>& uid,
      const std::vector<PositionType>& newPos,
      std::function<std::vector<bool>(std::vector<T>&)> updateFunc) {
    std::vector<T> out(pos.size());
    BatchUpdate(pos, uid, newPos, updateFunc, out);
  }

  void BatchUpdate(const std::vector<PositionType>& pos,
                   const std::vector<UidType>& uid,
                   const std::vector<PositionType>& newPos,
                   std::function<std::vector<bool>(std::vector<T>&)> updateFunc,
                   std::vector<T>& out) {
    uint64_t batchSize = pos.size();
    Assert(batchSize == uid.size());
    Assert(batchSize == newPos.size());
    Assert(batchSize == out.size());
    for (uint64_t i = 0; i < batchSize; ++i) {
      int len = ReadPath(pos[i], path);
      ReadElementAndRemoveFromPath(path, uid[i], out[i]);
      EvictPath(path, pos[i], len);
      WriteBackPath(path, pos[i]);
    }
    const std::vector<bool>& writeBackFlags = updateFunc(out);
    for (uint64_t i = 0; i < batchSize; ++i) {
      UidType u = DUMMY<UidType>();
      obliMove(writeBackFlags[i], u, uid[i]);
      Write<1>(u, out[i], newPos[i]);
    }
  }

  int ReadPath(PositionType pos, std::vector<Block_>& path) {
    Assert(path.size() == stashSize + Z * depth);

    memcpy(&path[0], stash->blocks, stashSize * sizeof(Block_));

    int level = tree.ReadPath(pos, (Bucket_*)(&path[stashSize]));
#ifndef NDEBUG
    for (int i = stashSize; stashSize + Z * level; ++i) {
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
    return stashSize + Z * level;
  }

  void WriteBackPath(const std::vector<Block_>& path, PositionType pos) {
    memcpy(stash->blocks, &path[0], stashSize * sizeof(Block_));
    tree.WritePath(pos, (Bucket_*)(&path[stashSize]));
  }

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
    // for (int i = 0; i < path.size(); ++i) {
    //   printf("%ld ", path[i].isDummy() ? -1UL : (int64_t)path[i].position);
    //   if (i >= actualStashSize + Z && (path.size() - i) % Z == 1) {
    //     printf("\n");
    //   }
    // }
  }

  // requires uid to be sorted
  const std::vector<PositionType> GetRandNewPoses(uint64_t batchSize) {
    std::vector<PositionType> newPos(batchSize);
    for (int i = 0; i < batchSize; ++i) {
      newPos[i] = UniformRandom(size() - 1);
    }
    return newPos;
  }

  const std::vector<PositionType> duplicateNewPoses(
      const std::vector<PositionType>& newPos,
      const std::vector<UidType>& uid) {
    std::vector<PositionType> dupNewPos(newPos);
    for (int i = 1; i < uid.size(); ++i) {
      obliMove(uid[i] == uid[i - 1], dupNewPos[i], dupNewPos[i - 1]);
    }
    return dupNewPos;
  }

  // requires uid to be sorted
  void deDuplicatePoses(std::vector<PositionType>& pos,
                        const std::vector<UidType>& uid) {
    for (int i = 1; i < pos.size(); ++i) {
      PositionType randPos = UniformRandom(size() - 1);
      obliMove(uid[i] == uid[i - 1], pos[i], randPos);
    }
  }

  std::vector<Block_> getTopKLevel(int k) {
    std::vector<Block_> topK(stashSize + Z * ((1UL << k) - 1));
    memcpy(&topK[0], stash->blocks, stashSize * sizeof(Block_));
    tree.getTopKLevel(k, (Bucket_*)(&topK[stashSize]));
    return topK;
  }

  std::vector<Block_> packCopiesTopKHelper(
      const std::vector<std::vector<Block_>>& topKCopies, int k, int level,
      PositionType offset) {
    if (k == level) {
      return topKCopies[offset];
    }
    const std::vector<Block_>& left =
        packCopiesTopKHelper(topKCopies, k, level + 1, offset);
    const std::vector<Block_>& right =
        packCopiesTopKHelper(topKCopies, k, level + 1, offset + (1UL << level));
    std::vector<Block_> result(left.size() + right.size());
    std::copy(left.begin(), left.end(), result.begin());
    std::copy(right.begin(), right.end(), result.begin() + left.size());
    EM::Algorithm::OrCompact(result.begin(), result.end(),
                             [](const auto& b) { return !b.isDummy(); });
    Assert(result[stashSize + Z * (level + 1)].isDummy());
    result.resize(stashSize + Z * (level + 1));  // truncate

    Bucket_ combinedBucket;
    for (int i = 0; i < Z; ++i) {
      combinedBucket.blocks[i] = result[i];
    }
    tree.SetByPathAndLevel(offset, level, combinedBucket);
    return std::vector<Block_>(result.begin() + Z, result.end());
  }

  void packCopiesTopK(const std::vector<std::vector<Block_>>& topKCopies,
                      int k) {
    const std::vector<Block_>& newStash =
        packCopiesTopKHelper(topKCopies, k, 0, 0);
    Assert(newStash.size() == stashSize);
    memcpy(stash->blocks, &newStash[0], stashSize * sizeof(Block_));
  }

  // void ParBatchReadAndRemove(int k, std::vector<PositionType>& pos,
  //                            const std::vector<UidType>& uid,
  //                            std::vector<T>& out) {}

  void duplicateVal(std::vector<T>& out, const std::vector<UidType>& uid) {
    for (int i = 1; i < out.size(); ++i) {
      obliMove(uid[i] == uid[i - 1], out[i], out[i - 1]);
    }
  }

  void ParBatchUpdate(
      std::vector<PositionType>& pos, const std::vector<UidType>& uid,
      const std::vector<PositionType>& newPos,
      std::function<std::vector<bool>(std::vector<T>&)> updateFunc,
      std::vector<T>& out, int numThreads) {
    uint64_t batchSize = pos.size();
    Assert(batchSize == uid.size());
    Assert(batchSize == newPos.size());
    Assert(batchSize == out.size());
    // std::cout << "Received batch requests" << std::endl;
    // for (uint64_t i = 0; i < batchSize; ++i) {
    //   std::cout << "pos = " << pos[i] << ", uid = " << uid[i] << std::endl;
    // }
    deDuplicatePoses(pos, uid);
    // std::cout << "\nDe-duplicated batch requests" << std::endl;
    // for (uint64_t i = 0; i < batchSize; ++i) {
    //   std::cout << "pos = " << pos[i] << ", uid = " << uid[i] << std::endl;
    // }
    const int k = GetLogBaseTwo(numThreads);
    const std::vector<Block_>& topK = getTopKLevel(k);
    const PositionType topKSize = topK.size();
    // std::cout << "Top K:" << std::endl;
    // for (const Block_& b : topK) {
    //   if (!b.isDummy()) {
    //     std::cout << b << ", ";
    //   }
    // }
    // std::cout << std::endl;
    // std::vector<uint8_t> found(uid.size(), 0);
    // EM::Algorithm::OrShuffle(topK);
    // std::unordered_map<UidType, PositionType>
    //     topKIndexer;  // TODO: use a secure cuckoo hash table and hide
    //     whether
    //                   // an element is found
    // for (int i = 0; i < topK.size(); ++i) {
    //   topKIndexer[topK[i].uid] = i;
    // }
    // for (int i = 0; i < uid.size(); ++i) {
    //   UidType u = uid[i];
    //   if (i > 0) {
    //     // if duplicate, replace with random
    //     obliMove(u == uid[i - 1], u, UniformRandom(DUMMY<UidType>() - 1));
    //   }
    //   // TODO replace the following to be oblivious
    //   auto it = topKIndexer.find(uid[i]);
    //   if (it != topKIndexer.end()) {
    //     out[i] = topK[it->second].data;
    //     topK[it->second].uid = DUMMY<UidType>();
    //     found[i] = 1;
    //   }
    // }
    // std::cout << "Top K after batched linear scan:" << std::endl;
    // for (const Block_& b : topK) {
    //   if (!b.isDummy()) {
    //     std::cout << b << ", ";
    //   }
    // }
    // std::cout << std::endl;
    // std::cout << std::endl;
    const PositionType numSubtree = 1UL << k;
    std::vector<std::vector<Block_>> topKCopies(numSubtree,
                                                std::vector<Block_>(topK));
    const int subStashSize = stashSize + Z * k;

#pragma omp parallel for num_threads(numSubtree)
    for (int subtreeIdx = 0; subtreeIdx < numSubtree; ++subtreeIdx) {
      auto& topKCopy = topKCopies[subtreeIdx];
      std::vector<PositionType> subPos;
      std::vector<UidType> subUid;
      std::vector<size_t> outputIdx;
      subPos.reserve(uid.size());
      subUid.reserve(uid.size());
      for (int i = 0; i < uid.size(); ++i) {
        // it's fine to reveal pos of each request
        if ((pos[i] ^ subtreeIdx) & (1UL << k) - 1) {
          continue;  // not in the subtree
        }
        subPos.push_back(pos[i]);
        subUid.push_back(uid[i]);
        outputIdx.push_back(i);
      }

      // std::cout << "Subtree " << subtreeIdx << " requests:" << std::endl;
      // for (int i = 0; i < subPos.size(); ++i) {
      //   std::cout << "pos = " << subPos[i] << ", uid = " << subUid[i]
      //             << std::endl;
      // }
      // filter elements in top k levels belonging to the subtree
      std::vector<int> prefixSum(topKSize + 1, 0);
      prefixSum[0] = 0;
      for (int i = 0; i < topKSize; ++i) {
        bool otherSubtreeFlag =
            (topKCopy[i].position ^ subtreeIdx) & (1UL << k) - 1;
        obliMove(otherSubtreeFlag, topKCopy[i].uid, DUMMY<UidType>());
        bool match = !topKCopy[i].isDummy();
        prefixSum[i + 1] = prefixSum[i] + match;
      }

      Assert(prefixSum.back() <= subStashSize);
      EM::Algorithm::OrCompactSeparateMark(topKCopy.begin(), topKCopy.end(),
                                           prefixSum.begin());
      topKCopy.resize(subStashSize);
      // std::cout << "Top K Copy after filtering:" << std::endl;
      // for (const Block_& b : topKCopy) {
      //   if (!b.isDummy()) {
      //     std::cout << b << ", ";
      //   }
      // }
      // std::cout << std::endl << std::endl;
      if (!subPos.empty()) {
        // reads elements in the subtree
        std::vector<Block_> localPath(stashSize + Z * depth);
        std::copy(topKCopy.begin(), topKCopy.end(), localPath.begin());
        for (int i = 0; i < subPos.size(); ++i) {
          int actualLevel = tree.ReadSubPath(
              subPos[i], (Bucket_*)&(localPath[subStashSize]), k);
          T tempOut;  // avoid interprocessor write contention
          ReadElementAndRemoveFromPath(
              localPath.begin(),
              localPath.begin() + stashSize + Z * actualLevel, subUid[i],
              tempOut);
          // change to following if requests are already served by searching the
          // top K levels:
          // ReadElementAndRemoveFromPath(path.begin() +
          // subStashSize,
          //                              path.begin() + stashSize + Z *
          //                              actualLevel, subUid[i], tempOut);
          // obliMove(!found[outputIdx[i]], out[outputIdx[i]], tempOut);
          out[outputIdx[i]] = tempOut;
          // EvictPath(path, subPos[i], actualLevel - k, subStashSize, k);
          tree.WriteSubPath(subPos[i], (Bucket_*)&(localPath[subStashSize]), k);
        }
        std::copy(localPath.begin(), localPath.begin() + subStashSize,
                  topKCopy.begin());
      }
      // std::cout << "remain in top K Copy of subtree " << subtreeIdx << ": ";

      // for (auto b : topKCopy) {
      //   if (!b.isDummy()) {
      //     std::cout << b << ", ";
      //   }
      // }
      // std::cout << std::endl;
    }

    duplicateVal(out, uid);
    // std::cout << "\nRead results" << std::endl;
    // for (uint64_t i = 0; i < batchSize; ++i) {
    //   std::cout << out[i] << std::endl;
    // }
    const std::vector<bool>& writeBackFlags = updateFunc(out);
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
      toWrite[i].data = out[i];
    }

    PositionType evictCounterCopy = evictCounter;
#pragma omp parallel for num_threads(numSubtree)
    for (PositionType subtreeIdx = 0; subtreeIdx < numSubtree; ++subtreeIdx) {
      PositionType localEvictCounter =
          evictCounterCopy;  // avoid race condition
      std::vector<Block_>& topKCopy = topKCopies[subtreeIdx];
      std::vector<Block_> localPath(stashSize + Z * depth);
      std::copy(topKCopy.begin(), topKCopy.end(), localPath.begin());
      int freq = std::max(1UL, numSubtree / evict_freq);
      for (uint64_t i = 0; i < batchSize; ++i) {
        Block_ toWriteBlock = toWrite[i];
        bool otherSubtreeFlag =
            (toWriteBlock.position ^ subtreeIdx) & (numSubtree - 1);
        obliMove(otherSubtreeFlag, toWriteBlock.uid, DUMMY<UidType>());
        // if (!toWriteBlock.isDummy()) {
        //   std::cout << toWriteBlock << std::endl;
        // }

        // this part may be slow
        bool success =
            WriteNewBlockToTreeTop(localPath, toWriteBlock, subStashSize);
        Assert(success || toWriteBlock.isDummy());
        for (int i = 0; i < evict_freq + 1; ++i) {
          ++localEvictCounter;
          PositionType p = localEvictCounter % size();
          if (p % numSubtree == subtreeIdx) {
            // std::cout << "Subtree " << subtreeIdx << " Evicting at pos " << p
            //           << std::endl;
            int actualLevel =
                tree.ReadSubPath(p, (Bucket_*)&(localPath[subStashSize]), k);
            EvictPath(localPath, p, actualLevel - k, subStashSize, k);
            tree.WriteSubPath(p, (Bucket_*)&(localPath[subStashSize]), k);
          }
        }
      }
      std::copy(localPath.begin(), localPath.begin() + subStashSize,
                topKCopy.begin());
    }

    // combine top K Copies recursively, pack blocks closer to leaves
    packCopiesTopK(topKCopies, k);
    evictCounter += batchSize * (evict_freq + 1);
    // below is for debug
    // std::vector<Block_> topKAfterBatch = getTopKLevel(k);
    // std::cout << "Top K after batched read:" << std::endl;
    // for (const Block_& b : topKAfterBatch) {
    //   if (!b.isDummy()) {
    //     std::cout << b << ", ";
    //   }
    // }
    // std::cout << std::endl << std::endl;
  }

  void ParBatchUpdate(
      std::vector<PositionType>& pos, const std::vector<UidType>& uid,
      const std::vector<PositionType>& newPos,
      std::function<std::vector<bool>(std::vector<T>&)> updateFunc,
      int numThreads) {
    std::vector<T> out(pos.size());
    ParBatchUpdate(pos, uid, newPos, updateFunc, out, numThreads);
  }

  std::vector<PositionType> ParBatchUpdate(
      std::vector<PositionType>& pos, const std::vector<UidType>& uid,
      std::function<std::vector<bool>(std::vector<T>&)> updateFunc,
      int numThreads) {
    const std::vector<PositionType>& newPoses = GetRandNewPoses(uid.size());
    ParBatchUpdate(pos, uid, newPoses, updateFunc, numThreads);
    return duplicateNewPoses(newPoses, uid);
  }
};
}  // namespace ODSL::CircuitORAM