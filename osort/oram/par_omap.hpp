#pragma once

#include <omp.h>

#include "external_memory/algorithm/param_select.hpp"
#include "omap.hpp"
namespace ODSL {
template <typename K, typename V>
struct ParOMap {
  std::vector<OMap<K, V>> shards;
  uint8_t randSalt[16];

  static uint64_t maxQueryPerShard(uint64_t batchSize, uint64_t shardCount,
                                   double logFailProb = -40) {
    auto satisfy = [&](uint64_t n) {
      double logSf = EM::Algorithm::binomLogSf(n, batchSize, 1.0 / shardCount);
      return logSf < logFailProb;
    };
    return EM::Algorithm::lowerBound(divRoundUp(batchSize, shardCount),
                                     batchSize, satisfy);
  }

  ParOMap(uint64_t mapSize, uint64_t shardCount) : shards(shardCount) {
    uint64_t shardSize = maxQueryPerShard(mapSize, shardCount, -60);
    for (auto& shard : shards) {
      shard.SetSize(shardSize);
    }
    read_rand(randSalt, sizeof(randSalt));
  }

  void Init() {
    for (auto& shard : shards) {
      shard.Init();
    }
  }

  template <class Reader>
  void InitFromReaderInPlace(Reader& reader) {
    uint64_t shardCount = shards.size();
    std::vector<uint64_t> shardCounter(shardCount);
    uint64_t initSizePerShard =
        maxQueryPerShard(reader.size(), shardCount, -60);
    struct ShuffleElement {
      uint32_t shardIdx;
      std::pair<K, V> kvPair;
    };
    uint64_t shuffleVecSize = initSizePerShard * shardCount;
    EM::VirtualVector::VirtualReader<ShuffleElement> shuffleReader(
        shuffleVecSize, [&](uint64_t i) {
          ShuffleElement shuffleElement;
          if (i < reader.size()) {
            shuffleElement.kvPair = reader.read();
            uint64_t hash = secure_hash_with_salt(
                (uint8_t*)&shuffleElement.kvPair.first, sizeof(K), randSalt);
            uint32_t shardIdx = hash % shardCount;
            for (uint32_t j = 0; j < shardCount; ++j) {
              bool matchFlag = shardIdx == j;
              shardCounter[j] += matchFlag;
            }
          } else {
            shuffleElement.isDummy = true;
            bool selectedShardFlag = false;
            for (uint32_t j = 0; j < shardCount; ++j) {
              bool matchFlag =
                  (shardCounter[j] < initSizePerShard) & !selectedShardFlag;
              shardCounter[j] += matchFlag;
              obliMove(matchFlag, shuffleElement.shardIdx, j);
              selectedShardFlag |= matchFlag;
            }
          }
          return shuffleElement;
        });
    using ShardInitVec = EM::NonCachedVector::Vector<std::pair<K, V>>;
    std::vector<ShardInitVec*> shardInitVecs(shardCount);
    for (uint32_t i = 0; i < shardCount; ++i) {
      shardInitVecs[i] = new ShardInitVec(initSizePerShard);
    }
    using ShardWriter = ShardInitVec::Writer;
    std::vector<ShardWriter*> shardWriters(shardCount);
    for (uint32_t i = 0; i < shardCount; ++i) {
      shardWriters[i] =
          new ShardWriter(shardInitVecs[i]->begin(), shardInitVecs[i]->end());
    }
    EM::VirtualVector::VirtualWriter<ShuffleElement> shuffleWriter(
        shuffleVecSize, [&](uint64_t i, const ShuffleElement& elem) {
          uint32_t shardIdx = elem.shardIdx;
          shardWriters[shardIdx]->write(elem.kvPair);
        });
    // TODO: init each sub omap
  }

  template <class KeyIterator, class ShardIdxIterator>
  std::pair<std::vector<K>, std::vector<uint32_t>> extractKeyByShard(
      const KeyIterator keyInputBegin, const KeyIterator keyInputEnd,
      const ShardIdxIterator shardIdxBegin, uint32_t shardIdx,
      uint32_t shardSize) {
    uint64_t shardCount = shards.size();
    uint64_t batchSize = keyInputEnd - keyInputBegin;
    std::vector<K> keyOutput(keyInputBegin, keyInputEnd);
    std::vector<uint32_t> prefixSum(batchSize + 1);
    prefixSum[0] = 0;
    uint32_t accSum = 0;
    for (uint32_t i = 0; i < batchSize; ++i) {
      bool matchFlag = shardIdxBegin[i] == shardIdx;
      accSum += matchFlag;
      prefixSum[i + 1] = accSum;
    }
    EM::Algorithm::OrCompactSeparateMark(keyOutput.begin(), keyOutput.end(),
                                         prefixSum.begin());
    uint32_t numReal = prefixSum.back();
    shardSize =
        std::max(shardSize, numReal);  // in case num real > shard size, we have
                                       // to reveal the actual shard size
    if (shardSize < batchSize) {
      std::vector<K> keyOutputTmp(keyOutput.begin(),
                                  keyOutput.begin() + shardSize);
      return std::make_pair(keyOutputTmp, prefixSum);
    }
    return std::make_pair(keyOutput, prefixSum);
  }

  template <class KeyIterator>
  std::vector<uint32_t> getShardIndexVec(const KeyIterator keyBegin,
                                         const KeyIterator keyEnd,
                                         const std::vector<bool>& isDuplicate) {
    uint64_t shardCount = shards.size();
    uint64_t batchSize = keyEnd - keyBegin;
    std::vector<uint32_t> shardIndexVec(batchSize);
#pragma omp parallel for schedule(static)
    for (uint64_t i = 0; i < batchSize; ++i) {
      uint64_t keyHash =
          secure_hash_with_salt((uint8_t*)&keyBegin[i], sizeof(K), randSalt);
      shardIndexVec[i] = keyHash % shardCount;
      obliMove(isDuplicate[i], shardIndexVec[i], DUMMY<uint32_t>());
    }
    return shardIndexVec;
  }

  struct RecoveryInfo {
    std::vector<uint32_t> recoveryArray;
    std::vector<bool> isDuplicate;
    RecoveryInfo(uint32_t size) : recoveryArray(size), isDuplicate(size) {}
  };

  RecoveryInfo sortAndSuppressDuplicateKeysInFind(std::vector<K>& keyVec) {
    RecoveryInfo recoveryInfo(keyVec.size());
    for (uint32_t i = 0; i < keyVec.size(); ++i) {
      recoveryInfo.recoveryArray[i] = i;
    }
    EM::Algorithm::BitonicSortSepPayload(keyVec.begin(), keyVec.end(),
                                         recoveryInfo.recoveryArray.begin());
    recoveryInfo.isDuplicate[0] = false;
    for (uint32_t i = 1; i < keyVec.size(); ++i) {
      recoveryInfo.isDuplicate[i] = keyVec[i] == keyVec[i - 1];
    }
    return recoveryInfo;
  }

  RecoveryInfo sortAndSuppressDuplicateKeysInInsert(std::vector<K>& keyVec,
                                                    std::vector<V>& valueVec) {
    RecoveryInfo recoveryInfo(keyVec.size());
    for (uint32_t i = 0; i < keyVec.size(); ++i) {
      recoveryInfo.recoveryArray[i] = i;
    }
    auto customSwap = [&](bool cond, uint64_t idx0, uint64_t idx1) {
      obliSwap(cond, recoveryInfo.recoveryArray[idx0],
               recoveryInfo.recoveryArray[idx1]);
      obliSwap(cond, valueVec[idx0], valueVec[idx1]);
    };
    EM::Algorithm::BitonicSortCustomSwap(keyVec.begin(), keyVec.end(),
                                         customSwap);
    recoveryInfo.isDuplicate[0] = false;
    for (uint32_t i = 1; i < keyVec.size(); ++i) {
      recoveryInfo.isDuplicate[i] = keyVec[i] == keyVec[i - 1];
    }
    return recoveryInfo;
  }

  template <typename T, class OutIterator>
  void mergeShardOutput(const std::vector<std::vector<T>>& shardOutput,
                        const std::vector<uint32_t>& shardIndexVec,
                        OutIterator outBegin) {
    uint32_t batchSize = shardIndexVec.size();
#pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < batchSize; ++i) {
      for (uint32_t j = 0; j < shardOutput.size(); ++j) {
        bool matchFlag = shardIndexVec[i] == j;
        obliMove(matchFlag, *(outBegin + i), shardOutput[j][i]);
      }
    }
  }

  template <class KeyIterator, class ValueIterator>
  std::vector<uint8_t> findBatch(const KeyIterator keyBegin,
                                 const KeyIterator keyEnd,
                                 ValueIterator valueBegin) {
    uint64_t shardCount = shards.size();
    uint64_t batchSize = keyEnd - keyBegin;
    uint64_t shardSize = maxQueryPerShard(batchSize, shardCount);

    std::vector<K> keyVec(keyBegin, keyEnd);
    RecoveryInfo recoveryInfo = sortAndSuppressDuplicateKeysInFind(keyVec);
    std::vector<uint32_t> shardIndexVec = getShardIndexVec(
        keyVec.begin(), keyVec.end(), recoveryInfo.isDuplicate);
    std::vector<std::vector<V>> valueShard(shardCount,
                                           std::vector<V>(batchSize));
    std::vector<std::vector<uint8_t>> foundFlagsShard(
        shardCount, std::vector<uint8_t>(batchSize));

#pragma omp parallel for schedule(static)
    for (uint64_t i = 0; i < shardCount; ++i) {
      const auto& [keyShard, prefixSum] = extractKeyByShard(
          keyVec.begin(), keyVec.end(), shardIndexVec.begin(), i, shardSize);
      uint32_t numReal = prefixSum.back();
      //       omp_set_num_threads(2);
      // #pragma omp parallel for schedule(static)
      for (uint32_t j = 0; j < keyShard.size(); ++j) {
        foundFlagsShard[i][j] =
            shards[i].find(keyShard[j], valueShard[i][j], j >= numReal);
        // std::cout << "shard " << i << " find: " << keyShard[j] << " "
        //           << valueShard[j] << " real = " << (j < numReal)
        //           << " findRes = " << foundFlagsTmp.get()[j] <<
        // std::endl;
      }

      EM::Algorithm::OrDistributeSeparateMark(
          valueShard[i].begin(), valueShard[i].end(), prefixSum.begin());
      EM::Algorithm::OrDistributeSeparateMark(foundFlagsShard[i].begin(),
                                              foundFlagsShard[i].end(),
                                              prefixSum.begin());
    }
    mergeShardOutput(valueShard, shardIndexVec, valueBegin);
    std::vector<uint8_t> foundFlags(batchSize);
    mergeShardOutput(foundFlagsShard, shardIndexVec, foundFlags.begin());
    for (uint32_t i = 1; i < recoveryInfo.isDuplicate.size(); ++i) {
      obliMove(recoveryInfo.isDuplicate[i], *(valueBegin + i),
               *(valueBegin + i - 1));
    }
    for (uint32_t i = 1; i < recoveryInfo.isDuplicate.size(); ++i) {
      obliMove(recoveryInfo.isDuplicate[i], foundFlags[i], foundFlags[i - 1]);
    }
    auto customSwap = [&](bool cond, uint64_t idx0, uint64_t idx1) {
      obliSwap(cond, *(valueBegin + idx0), *(valueBegin + idx1));
      obliSwap(cond, foundFlags[idx0], foundFlags[idx1]);
    };
    EM::Algorithm::BitonicSortCustomSwap(recoveryInfo.recoveryArray.begin(),
                                         recoveryInfo.recoveryArray.end(),
                                         customSwap);
    return foundFlags;
  }

  template <class KeyIterator, class ValueIterator>
  std::vector<uint8_t> insertBatch(const KeyIterator keyBegin,
                                   const KeyIterator keyEnd,
                                   const ValueIterator valueBegin) {
    uint64_t shardCount = shards.size();
    uint64_t batchSize = keyEnd - keyBegin;
    uint64_t shardSize = maxQueryPerShard(batchSize, shardCount);
    std::vector<uint32_t> shardIndexVec =
        getShardIndexVec(keyBegin, keyEnd, std::vector<bool>(batchSize, false));
    std::vector<std::vector<uint8_t>> foundFlagsShard(
        shardCount, std::vector<uint8_t>(batchSize));
#pragma omp parallel for schedule(static)
    for (uint64_t i = 0; i < shardCount; ++i) {
      const auto& [keyShard, prefixSum] = extractKeyByShard(
          keyBegin, keyEnd, shardIndexVec.begin(), i, shardSize);
      std::vector<V> valueShard(valueBegin, valueBegin + batchSize);
      EM::Algorithm::OrCompactSeparateMark(valueShard.begin(), valueShard.end(),
                                           prefixSum.begin());
      uint32_t numReal = prefixSum.back();
      for (uint32_t j = 0; j < keyShard.size(); ++j) {
        foundFlagsShard[i][j] =
            shards[i].insert(keyShard[j], valueShard[j], j >= numReal);
        // std::cout << "shard " << i << " insert: " << keyShard[j] << " "
        //           << valueShard[j] << " real = " << (j < numReal)
        //           << " insertRes = " << foundFlagsTmp.get()[j] << std::endl;
      }
      EM::Algorithm::OrDistributeSeparateMark(foundFlagsShard[i].begin(),
                                              foundFlagsShard[i].end(),
                                              prefixSum.begin());
    }
    std::vector<uint8_t> foundFlags(batchSize);
    mergeShardOutput(foundFlagsShard, shardIndexVec, foundFlags.begin());
    return foundFlags;
  }
};
}  // namespace ODSL