#pragma once
#include <functional>
#include <vector>

#include "bucket.hpp"
#include "common/probability.hpp"
#include "external_memory/algorithm/bitonic.hpp"
#include "external_memory/algorithm/kway_butterfly_sort.hpp"
#include "external_memory/algorithm/or_compact_shuffle.hpp"
#include "external_memory/noncachedvector.hpp"
#include "external_memory/stdvector.hpp"
#include "external_memory/virtualvector.hpp"
#include "tree.hpp"

namespace ORAM::PathORAM {
template <typename T, const int Z = 5, const int stashSize = 63,
          typename PositionType = uint64_t, typename UidType = uint64_t>
struct PathORAM {
  // using Node_ = Node<T, Z, PositionType, UidType>;
  using Stash = Bucket<T, stashSize, PositionType, UidType>;
  using Block_ = Block<T, PositionType, UidType>;
  using Bucket_ = Bucket<T, Z, PositionType, UidType>;
  using UidBlock_ = UidBlock<T, UidType>;
  // Node_* root = nullptr;
  HeapTree<Bucket_> tree;
  Stash* stash = nullptr;
  int depth = 0;
  int cacheLevel = 62;
  PositionType _size = 0;

  // PathORAM(const PathORAM& other) = default;
  PathORAM(PathORAM&& other) {
    _size = other._size;
    tree = std::move(other.tree);
    stash = other.stash;
    other.stash = nullptr;
    depth = other.depth;
    cacheLevel = other.cacheLevel;
  }

  PathORAM(PositionType size, int cacheLevel = 62)
      : _size(size), tree(size, cacheLevel), cacheLevel(cacheLevel) {
    stash = new Stash();
    depth = GetLogBaseTwo(size - 1) + 2;
  }

  template <typename Reader, class PosMapWriter>
  void InitFromReader(Reader& reader, PosMapWriter& posMapWriter) {
    size_t initSize = reader.size();
    // printf("Init from reader\n");
    size_t numBucket = 2 * size() - 1;
    size_t numBlock = numBucket * Z;
    // StdVector<char> loadVec(numBucket);
    struct Positions {
      PositionType pos[Z];
      Positions() {
        for (int i = 0; i < Z; ++i) {
          pos[i] = DUMMY<PositionType>();
        }
      }
    };
    struct PositionStash {
      PositionType pos[16 - Z];
      PositionStash() {
        for (int i = 0; i < 16 - Z; ++i) {
          pos[i] = DUMMY<PositionType>();
        }
      }
    };
    // HeapTree<
    NoReplaceSampler sampler(initSize, size());
    std::function<PositionStash(Positions&, const PositionStash&,
                                const PositionStash&)>
        reduceFunc = [](Positions& root, const PositionStash& left,
                        const PositionStash& right) -> PositionStash {
      PositionType combinedStash[16];
      for (int i = 0; i < 16; ++i) {
        combinedStash[i] = DUMMY<PositionType>();
      }
      for (int i = 0; i < 16 - Z; ++i) {
        obliMove(left.pos[i] != DUMMY<PositionType>(), combinedStash[i],
                 left.pos[i]);
      }
      for (int i = 0; i < 16 - Z; ++i) {
        obliMove(right.pos[i] != DUMMY<PositionType>(), combinedStash[15 - i],
                 right.pos[i]);
      }
      // TODO check that no element is lost

      EM::Algorithm::BitonicMergePow2(
          combinedStash, combinedStash + 16,
          [](auto a, auto b) { return a < b; }, true);
      memcpy(root.pos, combinedStash, sizeof(PositionType) * Z);
      return *(PositionStash*)(combinedStash + Z);
    };

    std::function<PositionStash(Positions&, size_t)> leafFunc =
        [&](Positions& leaf, size_t path) {
          for (int i = 0; i < Z; ++i) {
            leaf.pos[i] = DUMMY<PositionType>();
          }
          int num = sampler.Sample();
          PositionStash stash;
          for (int i = 0; i < Z; ++i) {
            obliMove(i < num, leaf.pos[i], path);
          }
          for (int i = Z; i < 16; ++i) {
            obliMove(i < num, stash.pos[i - Z], path);
          }
          return stash;
        };
    using PosVec = StdVector<PositionType>;
    PosVec positionVec(numBlock);
    PosVec prefixSum(numBlock + 1);
    // printf("Generate positions bottom up\n");
    {
      HeapTree<Positions> positions(size(), cacheLevel, 1, 1UL << 62);
      positions.template BuildBottomUp<PositionStash>(reduceFunc, leafFunc);
      // for (size_t i = 0; i < 2 * size - 1; ++i) {
      //   for (int j = 0; j < Z; ++j) {
      //     printf("%ld ", (int64_t)positions.GetByInternalIdx(i).pos[j]);
      //   }
      //   printf("\n");
      // }
      // printf("calc prefix sum\n");
      prefixSum[0] = 0;
      for (PositionType i = 0; i < numBucket; ++i) {
        const Positions& posBucket = positions.GetByInternalIdx(i);
        for (PositionType j = 0; j < Z; ++j) {
          prefixSum[i * Z + j + 1] =
              prefixSum[i * Z + j] +
              (posBucket.pos[j] != DUMMY<PositionType>());
          positionVec[i * Z + j] = posBucket.pos[j];
        }
      }
    }
    // EM::Algorithm::OrDistributeSeparateMark(
    if (prefixSum[numBlock] != initSize) {
      printf("prefixSum %lu initSize %lu\n", (int64_t)prefixSum[numBlock],
             (int64_t)initSize);
      throw std::runtime_error("Stash overflows.");
    }
    using DistributeVec = StdVector<UidBlock_>;
    using DistributeReader = typename DistributeVec::Reader;
    using DistributeWriter = typename DistributeVec::Writer;
    using UidVec = StdVector<UidType>;
    using UidReader = typename UidVec::Reader;
    using UidWriter = typename UidVec::Writer;

    DistributeVec distributeVec(numBlock);
    DistributeWriter distributeInputWriter(distributeVec.begin(),
                                           distributeVec.begin() + initSize);
    UidVec uidVec(initSize);
    UidWriter uidWriter(uidVec.begin(), uidVec.end());
    PositionType inputIdx = 0;
    EM::VirtualVector::WrappedReader<UidBlock_, Reader> shuffleReader(
        reader, [&](const T& val) { return UidBlock_(val, inputIdx++); });
    EM::VirtualVector::WrappedWriter<UidBlock_, UidWriter> shuffleWriter(
        uidWriter, [&](const UidBlock_& block) {
          distributeInputWriter.write(block);
          return block.uid;
        });
    // printf("shuffle elements\n");
    EM::Algorithm::KWayButterflyOShuffle(shuffleReader, shuffleWriter);
    distributeInputWriter.flush();
    // printf("distribute elements\n");
    EM::Algorithm::OrDistributeSeparateMark(
        distributeVec.begin(), distributeVec.end(), prefixSum.begin());
    DistributeReader distributeOutputReader(distributeVec.begin(),
                                            distributeVec.end());
    for (size_t i = 0; i < numBucket; ++i) {
      Bucket_ bucket;
      for (int j = 0; j < Z; ++j) {
        const UidBlock_& uidBlock = distributeOutputReader.read();
        PositionType pos = positionVec[i * Z + j];
        UidType uid = uidBlock.uid;
        // obliMove(pos == DUMMY<PositionType>(), uid, DUMMY<UidType>());
        bucket.blocks[j] = Block_(uidBlock.data, pos, uid);
      }
      tree.SetByInternalIdx(i, bucket);
    }

    // printf("compact positions\n");
    EM::Algorithm::OrCompactSeparateMark(positionVec.begin(), positionVec.end(),
                                         prefixSum.begin());
    for (size_t i = 0; i < positionVec.size(); ++i) {
    }
    size_t posMapIdx = 0;

    UidReader uidReader(uidVec.begin(), uidVec.end());

    EM::VirtualVector::WrappedReader<UidBlock<PositionType>, UidReader>
        posMapReader(uidReader, [&](const UidType& uid) {
          return UidBlock<PositionType>(positionVec[posMapIdx++], uid);
        });
    // printf("Final sort\n");
    EM::Algorithm::KWayButterflySort(posMapReader, posMapWriter);

    /* A quick test for reduceFunc
PositionStash testLeft, testRight;
int leftSize = 6;
int rightSize = 5;
for (int i = 0; i < 16 - Z; ++i) {
  if (i < leftSize) {
    testLeft.pos[i] = i;
  } else {
    testLeft.pos[i] = DUMMY<PositionType>();
  }
  if (i < rightSize) {
    testRight.pos[i] = i;
  } else {
    testRight.pos[i] = DUMMY<PositionType>();
  }
}
Positions testRoot;
PositionStash testRootStash = reduceFunc(testRoot, testLeft, testRight);
for (int i = 0; i < 16 - Z; ++i) {
  printf("%d ", (int)testRootStash.pos[i]);
}
printf("\n");
for (int i = 0; i < Z; ++i) {
  printf("%d ", (int)testRoot.pos[i]);
}

=> prints 2 3 3 4 4 5 -1 -1 -1 -1 -1
          0 0 1 1 2 2
*/
  }

  size_t size() const { return _size; }

  ~PathORAM() {
    if (stash) {
      delete stash;
      stash = nullptr;
    }
  }

  PathORAM& operator=(const PathORAM& other) = default;
  PathORAM& operator=(PathORAM&& other) = default;

  PositionType Read(PositionType pos, const UidType& uid, T& out) {
    PositionType newPos = UniformRandom(size() - 1);
    std::vector<Block_> path = ReadPath(pos);

    // for (int i = 0; i < path.size(); ++i) {
    //   printf("(%ld %ld)", (int64_t)path[i].uid, (int64_t)path[i].position);
    //   if (i >= stashSize - 1 && (path.size() - i - 1) % Z == 0) {
    //     printf(", ");
    //   }
    // }
    // printf("\n");
    ReadElementAndRemoveFromPath(path, uid, out);

    Block_ newBlock(out, newPos, uid);
    WriteNewBlockToPath(path, newBlock);
    EvictPath(path, pos);
    // printf("write %ld to position %ld, evict path %lu\n", (int64_t)uid,
    //        (int64_t)newPos, (int64_t)pos);
    // for (int i = 0; i < path.size(); ++i) {
    //   printf("uid %ld pos %ld\n", (int64_t)path[i].uid,
    //          (int64_t)path[i].position);
    // }
    // printf("\n");
    WriteBackPath(path, pos);
    return newPos;
  }

  PositionType Write(const UidType& uid, const T& in) {
    PositionType pos = UniformRandom(size() - 1);
    PositionType newPos = UniformRandom(size() - 1);
    std::vector<Block_> path = ReadPath(pos);
    Block_ newBlock(in, newPos, uid);
    WriteNewBlockToPath(path, newBlock);
    EvictPath(path, pos);
    WriteBackPath(path, pos);
    return newPos;
  }

  PositionType Update(PositionType pos, const UidType& uid,
                      std::function<void(T&)> updateFunc) {
    T out;
    return Update(pos, uid, updateFunc, out);
  }

  PositionType Update(PositionType pos, const UidType& uid,
                      std::function<void(T&)> updateFunc, T& out) {
    return Update(pos, uid, updateFunc, out, uid);  // does not change uid
  }

  PositionType Update(PositionType pos, const UidType& uid,
                      std::function<void(T&)> updateFunc, T& out,
                      const UidType& updatedUid) {
    PositionType newPos = UniformRandom(size() - 1);
    std::vector<Block_> path = ReadPath(pos);

    ReadElementAndRemoveFromPath(path, uid, out);
    updateFunc(out);
    Block_ newBlock(out, newPos, updatedUid);
    WriteNewBlockToPath(path, newBlock);
    EvictPath(path, pos);

    WriteBackPath(path, pos);
    return newPos;
  }

  std::vector<Block_> ReadPath(PositionType pos) {
    std::vector<Block_> path(stashSize + Z * depth);

    memcpy(&path[0], stash->blocks, stashSize * sizeof(Block_));

    size_t level = tree.ReadPath(pos, (Bucket_*)(&path[stashSize]));
    path.resize(stashSize + Z * level);
    // for (int i = stashSize; i < path.size(); ++i) {
    //   printf("%ld ", (int64_t)path[i].position);
    // }
    // printf("\n");
#ifndef NDEBUG
    for (int i = stashSize; i < path.size(); ++i) {
      int lv = (i - stashSize) / Z;
      size_t mask = (1UL << lv) - 1;
      Assert(path[i].uid == DUMMY<UidType>() ||
                 (path[i].position & mask) == (pos & mask),
             "Path position mismatch");
    }
#endif
    return path;
  }

  void WriteBackPath(const std::vector<Block_>& path, PositionType pos) {
    memcpy(stash->blocks, &path[0], stashSize * sizeof(Block_));
    tree.WritePath(pos, (Bucket_*)(&path[stashSize]));
    // for (int i = 0; i < path.size(); ++i) {
    //   printf("%ld ", (int64_t)path[i].uid);
    //   if (i >= stashSize - 1 && (path.size() - i - 1) % Z == 0) {
    //     printf(", ");
    //   }
    // }
    // printf("\n");
  }

  static void ReadElementAndRemoveFromPath(std::vector<Block_>& path,
                                           const UidType& uid, T& out) {
    for (Block_& b : path) {
      b.invalidateAndCopyDataIfUidMatch(uid, out);
    }
  }

  template <const bool toStashOnly = false>
  static void WriteNewBlockToPath(std::vector<Block_>& path,
                                  const Block_& block) {
    int endIdx = toStashOnly ? stashSize + Z : path.size();
    bool cond = true;
    // fill the first slot that's empty
    for (int i = 0; i < endIdx; i++) {
      cond &= !path[i].conditionalFillIfDummy(cond, block);
    }
    Assert(!cond, "No empty slot in path");
  }

  static void EvictPath(std::vector<Block_>& path, PositionType pos) {
    std::vector<int> deepest(path.size());
    std::vector<int> deeperDummyCount(path.size());
    int depth_ = (path.size() - stashSize) / Z;
    for (int i = 0; i < path.size(); i++) {
      deepest[i] = commonSuffixLength(path[i].position, pos);
      obliMove(path[i].isDummy(), deepest[i], -1);
      obliMove(deepest[i] >= depth_, deepest[i], depth_ - 1);
    }
    EM::Algorithm::BitonicSortSepPayload(deepest.begin(), deepest.end(),
                                         path.begin());
    int rank = 0;
    int numDummy = 0;
    for (int i = path.size() - 1; i >= 0; i--) {
      numDummy += path[i].isDummy();
      // least number of blocks deeper than this block
      int minRank = (depth_ - 1 - deepest[i]) * Z;
      obliMove(rank < minRank, rank, minRank);
      // # dummies that needs to be added deeper than this block
      deeperDummyCount[i] = rank - (path.size() - 1 - i);
      if ((deepest[i] >= 0) & (rank >= path.size())) {
        throw std::runtime_error("Stash overflow");
      }
      obliMove(deepest[i] == -1, deeperDummyCount[i],
               0);  // mark dummy as 0, so that it won't be moved
      ++rank;
    }

    // for (int i = 0; i < path.size(); ++i) {
    //   printf("%d ", deeperDummyCount[i]);
    //   if (i >= stashSize - 1 && (path.size() - i - 1) % Z == 0) {
    //     printf("\n");
    //   }
    // }
    // printf("\n");

    int distriLevel = GetLogBaseTwo(path.size() - 1) + 1;

    for (int level = distriLevel - 1; level >= 0; --level) {
      int pow2Level = 1 << level;
      for (int i = pow2Level; i < path.size(); ++i) {
        bool shiftFlag =
            (deeperDummyCount[i] & (2 * pow2Level - 1)) >= pow2Level;
        obliSwap(shiftFlag, deeperDummyCount[i - pow2Level],
                 deeperDummyCount[i]);
        // could uncomment this instead of OrDistribute
        // causes twice the swap
        // obliSwap(shiftFlag, path[i - pow2Level], path[i]);
      }
    }
    std::vector<int> markArr(path.size() + 1);
    markArr[0] = path.size() - numDummy;
    for (int i = 0; i < path.size(); ++i) {
      --numDummy;
      obliMove(numDummy < deeperDummyCount[i], numDummy, deeperDummyCount[i]);
      markArr[i + 1] = path.size() - 1 - i - numDummy;
    }
    EM::Algorithm::OrDistributeSeparateMark(path.rbegin(), path.rend(),
                                            markArr.rbegin());
  }

  static int commonSuffixLength(PositionType a, PositionType b) {
    return std::countr_zero(a ^ b);
  }

 private:
};
}  // namespace ORAM::PathORAM