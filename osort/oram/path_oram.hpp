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
  PositionType size = 0;

  PathORAM(const PathORAM& other) = default;
  PathORAM(PathORAM&& other) = default;

  PathORAM(PositionType size, int cacheLevel = 62)
      : size(size), tree(size, cacheLevel), cacheLevel(cacheLevel) {
    stash = new Stash();
    depth = GetLogBaseTwo(size - 1) + 2;
  }

  template <typename Reader, class PosMapWriter>
  void InitFromReader(Reader& reader, PosMapWriter& posMapWriter) {
    size_t numBucket = 2 * size - 1;
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
    NoReplaceSampler sampler(size);
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
    {
      HeapTree<Positions> positions(size, cacheLevel, 1, 1UL << 62);
      positions.template BuildBottomUp<PositionStash>(reduceFunc, leafFunc);
      // for (size_t i = 0; i < 2 * size - 1; ++i) {
      //   for (int j = 0; j < Z; ++j) {
      //     printf("%ld ", (int64_t)positions.GetByInternalIdx(i).pos[j]);
      //   }
      //   printf("\n");
      // }

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
    Assert(prefixSum[numBlock] == size);
    using DistributeVec = StdVector<UidBlock_>;
    using DistributeReader = typename DistributeVec::Reader;
    using DistributeWriter = typename DistributeVec::Writer;
    using UidVec = StdVector<UidType>;
    using UidReader = typename UidVec::Reader;
    using UidWriter = typename UidVec::Writer;

    DistributeVec distributeVec(numBlock);
    DistributeWriter distributeInputWriter(distributeVec.begin(),
                                           distributeVec.begin() + size);
    UidVec uidVec(size);
    UidWriter uidWriter(uidVec.begin(), uidVec.end());
    PositionType heapIdx = 0;
    PositionType inputIdx = 0;
    EM::VirtualVector::WrappedReader<UidBlock_, Reader> shuffleReader(
        reader, [&](const T& val) { return UidBlock_(val, inputIdx++); });
    EM::VirtualVector::WrappedWriter<UidBlock_, UidWriter> shuffleWriter(
        uidWriter, [&](const UidBlock_& block) {
          distributeInputWriter.write(block);
          return block.uid;
        });
    EM::Algorithm::KWayButterflyOShuffle(shuffleReader, shuffleWriter);
    distributeInputWriter.flush();

    EM::Algorithm::OrDistributeSeparateMark(
        distributeVec.begin(), distributeVec.end(), prefixSum.begin());
    DistributeReader distributeOutputReader(distributeVec.begin(),
                                            distributeVec.end());
    for (size_t i = 0; i < numBucket; ++i) {
      Bucket_ bucket;
      for (int j = 0; j < Z; ++j) {
        const UidBlock_& uidBlock = distributeOutputReader.read();
        bucket.blocks[j] =
            Block_(uidBlock.data, positionVec[i * Z + j], uidBlock.uid);
      }
      tree.SetByInternalIdx(i, bucket);
    }
    EM::Algorithm::OrCompactSeparateMark(positionVec.begin(), positionVec.end(),
                                         prefixSum.begin());

    size_t posMapIdx = 0;

    UidReader uidReader(uidVec.begin(), uidVec.end());

    EM::VirtualVector::WrappedReader<UidBlock<PositionType>, UidReader>
        posMapReader(uidReader, [&](const UidType& uid) {
          return UidBlock<PositionType>(positionVec[posMapIdx++], uid);
        });

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

  void Destroy() {
    if (stash) {
      delete stash;
      stash = nullptr;
    }
  }

  ~PathORAM() { Destroy(); }

  PathORAM& operator=(const PathORAM& other) = default;
  PathORAM& operator=(PathORAM&& other) = default;

  PositionType Read(PositionType pos, const UidType& uid, T& out) {
    PositionType newPos = UniformRandom(size - 1);
    std::vector<Block_> path = ReadPath(pos);

    // for (int i = 0; i < path.size(); ++i) {
    //   printf("%ld ", (int64_t)path[i].uid);
    //   if (i >= stashSize - 1 && (path.size() - i - 1) % Z == 0) {
    //     printf(", ");
    //   }
    // }
    // printf("\n");
    ReadElementAndRemoveFromPath(path, uid, out);

    Block_ newBlock(out, newPos, uid);
    WriteNewBlockToPath(path, newBlock);
    EvictPath(path, pos);

    WriteBackPath(path, pos);
    return newPos;
  }

  PositionType Write(const UidType& uid, const T& in) {
    PositionType pos = UniformRandom(size - 1);
    PositionType newPos = UniformRandom(size - 1);
    std::vector<Block_> path = ReadPath(pos);
    Block_ newBlock(in, newPos, uid);
    WriteNewBlockToPath(path, newBlock);
    EvictPath(path, pos);
    // printf("write %lu to position %lu, evict path %lu\n", uid, newPos, pos);
    // for (int i = 0; i < path.size(); ++i) {
    //   printf("%lu ", path[i].uid);
    // }
    // printf("\n");
    WriteBackPath(path, pos);
    return newPos;
  }

  PositionType Update(PositionType pos, const UidType& uid,
                      std::function<T(T)> updateFunc) {
    T out;
    return Update(pos, uid, updateFunc, out);
  }

  PositionType Update(PositionType pos, const UidType& uid,
                      std::function<T(T)> updateFunc, T& out) {
    return Update(pos, uid, updateFunc, out, uid);  // does not change uid
  }

  PositionType Update(PositionType pos, const UidType& uid,
                      std::function<T(T)> updateFunc, T& out,
                      const UidType& updatedUid) {
    PositionType newPos = UniformRandom(size - 1);
    std::vector<Block_> path = ReadPath(pos);

    ReadElementAndRemoveFromPath(path, uid, out);
    out = updateFunc(out);
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

    // for (int i = 0; i < path.size(); ++i) {
    //   printf("%d ", deeperDummyCount[i]);
    //   if (i >= stashSize - 1 && (path.size() - i - 1) % Z == 0) {
    //     printf("\n");
    //   }
    // }
    // printf("\n");

    // for (int i = 0; i < markArr.size(); ++i) {
    //   printf("%d ", markArr[i]);
    //   if (i >= stashSize && (markArr.size() - i - 1) % Z == 0) {
    //     printf("\n");
    //   }
    // }
    // printf("\n");

    // for (int i = 0; i < path.size(); ++i) {
    //   printf("%d ",
    //          path[i].isDummy() ? 0 : commonSuffixLength(path[i].position,
    //          pos));
    //   if (i >= stashSize - 1 && (path.size() - i - 1) % Z == 0) {
    //     printf("\n");
    //   }
    // }
    // printf("\n");
  }

  static int commonSuffixLength(PositionType a, PositionType b) {
    return std::countr_zero(a ^ b);
  }

 private:
  // // TODO: make it cache efficient
  // Node_* initTree(size_t size) {
  //   Node_* root = new Node_();
  //   if (size == 1) {
  //     return root;
  //   }
  //   size_t rightSize = size >> 1;
  //   root->left = initTree(size - rightSize);
  //   root->right = initTree(rightSize);
  //   return root;
  // }

  // template <typename Iterator>
  // Node_* initTree(Iterator begin, Iterator end) {
  //   size_t size = end - begin;
  //   Node_* root = new Node_();
  //   if (size == 1) {
  //     root->bucket.blocks[0] = *begin;
  //     printf("init pos %lu uid %lu, key = %lu\n",
  //            root->bucket.blocks[0].position, root->bucket.blocks[0].uid,
  //            root->bucket.blocks[0].data.key);
  //     return root;
  //   }
  //   Iterator mid = begin + ((size + 1) >> 1);
  //   root->left = initTree(begin, mid);
  //   root->right = initTree(mid, end);
  //   return root;
  // }

  // void destroyTree(Node_* node) {
  //   if (node == nullptr) {
  //     return;
  //   }
  //   destroyTree(node->left);
  //   destroyTree(node->right);
  //   delete node;
  // }
};
}  // namespace ORAM::PathORAM