#pragma once

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
    if (shardSize < batchSize) {
      std::vector<K> keyOutputTmp(keyOutput.begin(),
                                  keyOutput.begin() + shardSize);
      return std::make_pair(keyOutputTmp, prefixSum);
    }
    return std::make_pair(keyOutput, prefixSum);
  }

  template <class KeyIterator>
  std::vector<uint32_t> getShardIndexVec(const KeyIterator keyBegin,
                                         const KeyIterator keyEnd) {
    uint64_t shardCount = shards.size();
    uint64_t batchSize = keyEnd - keyBegin;
    std::vector<uint32_t> shardIndexVec(batchSize);
    for (uint64_t i = 0; i < batchSize; ++i) {
      uint64_t keyHash =
          secure_hash_with_salt((uint8_t*)&keyBegin[i], sizeof(K), randSalt);
      shardIndexVec[i] = keyHash % shardCount;
    }
    return shardIndexVec;
  }

  template <class KeyIterator, class ValueIterator>
  std::vector<bool> findBatch(const KeyIterator keyBegin,
                              const KeyIterator keyEnd,
                              ValueIterator valueBegin) {
    uint64_t shardCount = shards.size();
    uint64_t batchSize = keyEnd - keyBegin;
    uint64_t shardSize = maxQueryPerShard(batchSize, shardCount);
    std::vector<uint32_t> shardIndexVec = getShardIndexVec(keyBegin, keyEnd);
    std::vector<bool> foundFlags(batchSize);
    for (uint64_t i = 0; i < shardCount; ++i) {
      const auto& [keyShard, prefixSum] = extractKeyByShard(
          keyBegin, keyEnd, shardIndexVec.begin(), i, shardSize);
      std::vector<V> valueShard(shardSize);
      uint32_t numReal = prefixSum.back();
      std::unique_ptr<bool> foundFlagsTmp(new bool[batchSize]);
      for (uint32_t j = 0; j < shardSize; ++j) {
        foundFlagsTmp.get()[j] =
            shards[i].find(keyShard[j], valueShard[j], j >= numReal);
        // std::cout << "shard " << i << " find: " << keyShard[j] << " "
        //           << valueShard[j] << " real = " << (j < numReal)
        //           << " findRes = " << foundFlagsTmp.get()[j] << std::endl;
      }
      std::vector<V> valueTmp(batchSize);
      std::copy(valueShard.begin(), valueShard.begin() + batchSize,
                valueTmp.begin());
      EM::Algorithm::OrDistributeSeparateMark(valueTmp.begin(), valueTmp.end(),
                                              prefixSum.begin());
      EM::Algorithm::OrDistributeSeparateMark(foundFlagsTmp.get(),
                                              foundFlagsTmp.get() + batchSize,
                                              prefixSum.begin());
      for (uint32_t j = 0; j < batchSize; ++j) {
        bool isReal = prefixSum[j] != prefixSum[j + 1];
        obliMove(isReal, valueBegin[j], valueTmp[j]);
        bool prevFoundFlag = foundFlags[j];
        obliMove(isReal, prevFoundFlag, foundFlagsTmp.get()[j]);
        foundFlags[j] = prevFoundFlag;
      }
    }
    return foundFlags;
  }

  template <class KeyIterator, class ValueIterator>
  std::vector<bool> insertBatch(const KeyIterator keyBegin,
                                const KeyIterator keyEnd,
                                const ValueIterator valueBegin) {
    uint64_t shardCount = shards.size();
    uint64_t batchSize = keyEnd - keyBegin;
    uint64_t shardSize = maxQueryPerShard(batchSize, shardCount);
    std::vector<uint32_t> shardIndexVec = getShardIndexVec(keyBegin, keyEnd);
    std::vector<bool> foundFlags(batchSize);
    for (uint64_t i = 0; i < shardCount; ++i) {
      const auto& [keyShard, prefixSum] = extractKeyByShard(
          keyBegin, keyEnd, shardIndexVec.begin(), i, shardSize);
      std::vector<V> valueShard(valueBegin, valueBegin + batchSize);
      EM::Algorithm::OrCompactSeparateMark(valueShard.begin(), valueShard.end(),
                                           prefixSum.begin());
      valueShard.resize(shardSize);
      uint32_t numReal = prefixSum.back();
      std::unique_ptr<bool> foundFlagsTmp(new bool[batchSize]);
      for (uint32_t j = 0; j < shardSize; ++j) {
        foundFlagsTmp.get()[j] =
            shards[i].insert(keyShard[j], valueShard[j], j >= numReal);
        // std::cout << "shard " << i << " insert: " << keyShard[j] << " "
        //           << valueShard[j] << " real = " << (j < numReal)
        //           << " insertRes = " << foundFlagsTmp.get()[j] << std::endl;
      }
      EM::Algorithm::OrDistributeSeparateMark(foundFlagsTmp.get(),
                                              foundFlagsTmp.get() + batchSize,
                                              prefixSum.begin());
      for (uint32_t j = 0; j < batchSize; ++j) {
        bool isReal = prefixSum[j] != prefixSum[j + 1];
        bool prevFoundFlag = foundFlags[j];
        obliMove(isReal, prevFoundFlag, foundFlagsTmp.get()[j]);
        foundFlags[j] = prevFoundFlag;
      }
    }
    return foundFlags;
  }
};
}  // namespace ODSL