#pragma once

#include <functional>
#include <vector>

#include "bucket.hpp"
#include "common/probability.hpp"
#include "external_memory/algorithm/kway_butterfly_sort.hpp"
#include "external_memory/algorithm/or_compact_shuffle.hpp"
#include "external_memory/noncachedvector.hpp"
#include "external_memory/stdvector.hpp"
#include "external_memory/virtualvector.hpp"
#include "heap_tree.hpp"
namespace ODSL {
template <typename T, const int Z = 5, const int stashSize = 63,
          typename PositionType = uint64_t, typename UidType = uint64_t,
          typename Reader, class PosMapWriter, class HeapTree_>
void FastInitFromReader(Reader& reader, PosMapWriter& posMapWriter,
                        HeapTree_& tree) {
  using Stash = Bucket<T, stashSize, PositionType, UidType>;
  using Block_ = Block<T, PositionType, UidType>;
  using Bucket_ = Bucket<T, Z, PositionType, UidType>;
  using UidBlock_ = UidBlock<T, UidType>;

  PositionType initSize = reader.size();
  // printf("initSize = %lu\n", initSize);

  PositionType numBucket = 2 * tree.GetLeafCount() - 1;
  PositionType numBlock = numBucket * Z;
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
  NoReplaceSampler sampler(initSize, tree.GetLeafCount());
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
        combinedStash, combinedStash + 16, [](auto a, auto b) { return a < b; },
        true);
    memcpy(root.pos, combinedStash, sizeof(PositionType) * Z);
    return *(PositionStash*)(combinedStash + Z);
  };

  std::function<PositionStash(Positions&, PositionType)> leafFunc =
      [&](Positions& leaf, PositionType path) {
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
    constexpr size_t max_node_per_page = HeapTree_::max_node_per_page;
    static constexpr int page_size = max_node_per_page * sizeof(Positions);
    constexpr int evict_freq = HeapTree_::evict_freq_;
    HeapTree<Positions, PositionType, page_size, evict_freq> positions(
        tree.GetLeafCount(), tree.GetCacheLevel(),
        2 * tree.GetLeafCount());  // ensure the structure is the same
    positions.template BuildBottomUp<PositionStash>(reduceFunc, leafFunc);
    prefixSum[0] = 0;
    for (PositionType i = 0; i < numBucket; ++i) {
      const Positions& posBucket = positions.GetByInternalIdx(i);
      for (PositionType j = 0; j < Z; ++j) {
        prefixSum[i * Z + j + 1] =
            prefixSum[i * Z + j] + (posBucket.pos[j] != DUMMY<PositionType>());
        positionVec[i * Z + j] = posBucket.pos[j];
      }
    }
  }
  // EM::Algorithm::OrDistributeSeparateMark(
  if (prefixSum[numBlock] != initSize) {
    throw std::runtime_error("Stash overflows.");
  }

  EM::Algorithm::OrCompactSeparateMark(positionVec.begin(), positionVec.end(),
                                       prefixSum.begin());

  typename PosVec::Reader positionReader(positionVec.begin(),
                                         positionVec.begin() + initSize);

  // using DistributeVec = StdVector<UidBlock_>;
  using DistributeVec = typename EM::VirtualVector::Vector<Block_>;
  using DistributeReader = typename DistributeVec::Reader;
  using DistributeWriter = typename DistributeVec::Writer;
  using UidVec = StdVector<UidType>;
  using UidReader = typename UidVec::Reader;
  using UidWriter = typename UidVec::Writer;
  UidVec uidVec(initSize);

  std::function<Block_&(size_t)> virtualize = [&](size_t idx) -> Block_& {
    return tree.GetMutableByInternalIdx(idx / Z).blocks[idx % Z];
  };
  // std::function<Block_&(size_t)> virtualize;
  DistributeVec distributeVec(numBlock, virtualize);

  DistributeWriter distributeInputWriter(distributeVec.begin(),
                                         distributeVec.begin() + initSize);

  UidWriter uidWriter(uidVec.begin(), uidVec.end());
  PositionType inputIdx = 0;
  EM::VirtualVector::WrappedReader<UidBlock_, Reader> shuffleReader(
      reader, [&](const T& val) { return UidBlock_(val, inputIdx++); });
  EM::VirtualVector::WrappedWriter<UidBlock_, UidWriter> shuffleWriter(
      uidWriter, [&](const UidBlock_& uidBlock) {
        Block_ block(uidBlock.data, positionReader.read(), uidBlock.uid);
        distributeInputWriter.write(block);
        return uidBlock.uid;
      });
  EM::Algorithm::KWayButterflyOShuffle(shuffleReader, shuffleWriter);
  distributeInputWriter.flush();

  EM::Algorithm::OrDistributeSeparateMark(
      distributeVec.begin(), distributeVec.end(), prefixSum.begin());

  PositionType posMapIdx = 0;

  UidReader uidReader(uidVec.begin(), uidVec.end());

  EM::VirtualVector::WrappedReader<UidBlock<PositionType, UidType>, UidReader>
      posMapReader(uidReader, [&](const UidType& uid) {
        return UidBlock<PositionType, UidType>(positionVec[posMapIdx++], uid);
      });

  EM::Algorithm::KWayButterflySort(posMapReader, posMapWriter);
}

template <typename PositionType>
int commonSuffixLength(PositionType a, PositionType b) {
  return std::countr_zero(a ^ b);
}

template <class PathVec, typename UidType, typename T>
void ReadElementAndRemoveFromPath(PathVec& path, const UidType& uid, T& out) {
  using Block_ = PathVec::value_type;
  for (Block_& b : path) {
    b.invalidateAndCopyDataIfUidMatch(uid, out);
  }
}

// Write a new block to the path, return false if the path is full
template <class PathVec, typename Block_>
bool WriteNewBlockToPath(PathVec& path, const Block_& block) {
  int endIdx = path.size();
  bool cond = true;
  // fill the first slot that's empty
  for (int i = 0; i < endIdx; i++) {
    cond &= !path[i].conditionalFillIfDummy(cond, block);
  }
  return !cond;
}

// Write a new block to the top of the tree, return false if the top is full
template <class PathVec, typename Block_>
bool WriteNewBlockToTreeTop(PathVec& path, const Block_& block, int topSize) {
  bool cond = true;
  // fill the first slot that's empty
  for (int i = 0; i < topSize; i++) {
    cond &= !path[i].conditionalFillIfDummy(cond, block);
  }
  return !cond;
}

template <typename T, const int Z, const int stashSize,
          typename PositionType = uint64_t, typename UidType = uint64_t>
int GetMaxCacheLevel(PositionType size, size_t cacheBytes = 1UL << 62) {
  using Bucket_ = Bucket<T, Z, PositionType, UidType>;
  using Stash = Bucket<T, stashSize, PositionType, UidType>;
  using HeapTree_ = HeapTree<Bucket_, PositionType>;
  int maxCacheLevel = 0;
  if (cacheBytes >= HeapTree_::GetMemoryUsage(size, 62)) {
    // default value
    return 62;
  }
  if (cacheBytes < sizeof(Stash)) {
    return -1;
  }
  cacheBytes -= sizeof(Stash);
  for (;; ++maxCacheLevel) {
    size_t treeUsage = HeapTree_::GetMemoryUsage(size, maxCacheLevel);
    if (cacheBytes < treeUsage) {
      break;
    }
  }
  return maxCacheLevel - 1;
}

}  // namespace ODSL
