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
  PositionType size = 0;

  PathORAM() = default;
  PathORAM(const PathORAM& other) = default;
  PathORAM(PathORAM&& other) = default;

  void Init(PositionType size) {
    tree.Init(size);
    stash = new Stash();
    depth = GetLogBaseTwo(size - 1) + 2;
    this->size = size;
  }

  template <typename Reader, class PosMapWriter>
  void InitFromReader(Reader& reader, PosMapWriter& posMapWriter) {
    size = reader.size();
    stash = new Stash();
    depth = GetLogBaseTwo(size - 1) + 2;
    tree.Init(size);
    StdVector<char> loadVec(2 * size - 1);
    NoReplaceSampler sampler(size);
    std::function<char(char&, char&)> reduceFunc = [](char& i,
                                                      char& j) -> char {
      char overflowLeft = 0;
      char overflowRight = 0;
      obliMove(i > Z, overflowLeft, (char)(i - Z));
      obliMove(j > Z, overflowRight, (char)(j - Z));
      i -= overflowLeft;
      j -= overflowRight;
      return overflowLeft + overflowRight;
    };
    std::function<void(char&)> leafFunc = [&](char& i) {
      i = sampler.Sample();
      // char overflow = 0;
      // obliMove(i > Z, overflow, (char)(i - Z));
      // i -= overflow;
      // return overflow;
    };
    char overflow =
        BuildBottomUp(loadVec.begin(), loadVec.end(), reduceFunc, leafFunc);
    for (int i = 0; i < loadVec.size(); ++i) {
      printf("%d ", (int)loadVec[i]);
    }
  }

  template <typename Vec, class Writer>
  void InitFromVector(Vec& vec, Writer& posMapWriter) {
    InitFromVector(vec, posMapWriter, 0, vec.size());
  }

  template <typename Vec, class Writer>
  void InitFromVector(Vec& vec, Writer& posMapWriter, UidType beginIdx,
                      UidType endIdx) {
    size = endIdx - beginIdx;

    stash = new Stash();
    StdVector<UidType> uidVec(size);
    depth = GetLogBaseTwo(size - 1) + 2;
    // EM::VirtualVector::Vector<UidBlock_, Vec> v(
    //     vec,
    //     [](size_t i, const T& val) { return UidBlock_(val, i); },  //
    //     virtualize
    //     [&](size_t i, const UidBlock_& block) {
    //       uidVec[i] = block.uid;
    //       return block.data;
    //     }  // devirtualize
    // );
    // EM::Algorithm::KWayButterflyOShuffle(v);

    for (auto val : vec) {
      printf("%lu ", val.key);
    }
    printf("\n");
    // for (const uint64_t& uid : uidVec) {
    //   printf("%lu ", uid);
    // }
    // printf("\n");

    // EM::VirtualVector::Vector<Block_, Vec> vBlock(
    //     vec,
    //     [&](size_t i, const T& val) {
    //       printf("set oram entry pos =  %lu uid = %lu\n", i, uidVec[i]);
    //       return Block_(val, i, uidVec[i]);
    //     },  // virtualize
    //     [](size_t i, const Block_& block) {
    //       return block.data;
    //     }  // devirtualize
    // );

    // root = initTree(vBlock.begin(), vBlock.end());

    // EM::VirtualVector::Vector<UidBlock<PositionType>, StdVector<UidType>>
    // vUid(
    //     uidVec,
    //     [](size_t i, const UidType& uid) { return UidBlock<UidType>(i, uid);
    //     },
    //     [&](size_t i, const UidBlock<UidType>& block) {
    //       posMapWriter.write(block.data);
    //       return block.uid;
    //     });
    // EM::Algorithm::KWayButterflySortInternal(vUid);
    // for (const uint64_t& pos : positionVec) {
    //   printf("%lu ", pos);
    // }
    // printf("\n");
  }

  void Destroy() {
    tree.Destroy();
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

  std::vector<Block_> ReadPath(PositionType pos) const {
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