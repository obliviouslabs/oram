#pragma once

#include <omp.h>

#include "cuckoo.hpp"
#include "external_memory/algorithm/param_select.hpp"
#include "omap.hpp"
namespace ODSL {
template <typename K, typename V, typename PositionType = uint64_t>
struct ParOMap {
  // std::vector<OMap<K, V, 9, PositionType>> shards;
  using BaseMap = CuckooHashMap<K, V, true, PositionType, true>;
  std::vector<BaseMap> shards;
  uint64_t shardSize = 0;
  uint8_t randSalt[16];

  uint64_t shardHashRange;

  INLINE uint32_t getShardByHash(uint64_t hash) {
    return hash / shardHashRange;
  }

  static uint64_t maxQueryPerShard(uint64_t batchSize, uint64_t shardCount,
                                   double logFailProb = -40) {
    auto satisfy = [&](uint64_t n) {
      double logSf = EM::Algorithm::binomLogSf(n, batchSize, 1.0 / shardCount);
      return logSf < logFailProb;
    };
    return EM::Algorithm::lowerBound(divRoundUp(batchSize, shardCount),
                                     batchSize, satisfy);
  }

  static uint64_t numRealPerBucket(uint64_t bucketSize, uint64_t shardCount,
                                   double logFailProb = -60) {
    auto satisfy = [&](uint64_t numDummy) {
      double logSf = EM::Algorithm::binomLogSf(
          bucketSize, (bucketSize - numDummy) * shardCount, 1.0 / shardCount);
      return logSf < logFailProb;
    };
    return bucketSize - EM::Algorithm::lowerBound(1UL, bucketSize - 1, satisfy);
  }

  ParOMap() {}

  ParOMap(uint64_t mapSize, uint64_t shardCount) {
      SetSize(mapSize, shardCount);
  }

  void SetSize(uint64_t mapSize, uint64_t shardCount) {
    if (shardCount == 0) {
      throw std::runtime_error("shardCount should be positive");
    }
    if (shardCount > 128) {
      throw std::runtime_error("shardCount should be no more than 128");
    }
    shards.resize(shardCount);
    shardSize = maxQueryPerShard(mapSize, shardCount, -60);
    read_rand(randSalt, sizeof(randSalt));
    shardHashRange = (UINT64_MAX - 1) / shardCount + 1;
  }

  void Init(size_t cacheBytes = ((uint64_t)ENCLAVE_SIZE << 20) * 3UL / 4UL) {
#pragma omp parallel for
    for (auto& shard : shards) {
      shard.SetSize(shardSize, cacheBytes / shards.size());
      shard.Init();
    }
  }


//   template <class Reader>
//   void InitFromReaderInPlace(Reader& reader) {
//     uint64_t shardCount = shards.size();
//     std::vector<uint64_t> shardCounter(shardCount);
//     uint64_t initSize = reader.size();

//     struct Element {
//       uint32_t shardIdx;
//       std::pair<K, V> kvPair;
//       bool operator<(const Element& other) const {
//         return shardIdx < other.shardIdx;
//         ;
//       }
//     };

//     EM::VirtualVector::VirtualReader<Element> sortReader(
//         initSize, [&](uint64_t i) {
//           Element sortElement;

//           sortElement.kvPair = reader.read();
//           uint64_t hash = secure_hash_with_salt(
//               (uint8_t*)&sortElement.kvPair.first, sizeof(K), randSalt);
//           sortElement.shardIdx = getShardByHash(hash);
//           for (uint32_t j = 0; j < shardCount; ++j) {
//             bool matchFlag = sortElement.shardIdx == j;
//             shardCounter[j] += matchFlag;
//           }
//           return sortElement;
//         });

//     using SortVec = EM::NonCachedVector::Vector<Element>;
//     SortVec sortVec(initSize);
//     typename SortVec::Writer sortWriter(sortVec.begin(), sortVec.end());
//     EM::Algorithm::KWayButterflySort(sortReader, sortWriter,
//                                      ((uint64_t)ENCLAVE_SIZE << 20) / 5);
//     std::vector<PositionType> shardPrefixSum(shardCount + 1);
//     shardPrefixSum[0] = 0;
//     for (int i = 0; i < shardCount; ++i) {
//       shardPrefixSum[i + 1] = shardPrefixSum[i] + shardCounter[i];
//     }
// #pragma omp parallel for schedule(static)
//     for (uint32_t i = 0; i < shardCount; ++i) {
//       typename SortVec::Reader shardReader(
//           sortVec.begin() + shardPrefixSum[i],
//           sortVec.begin() + shardPrefixSum[i + 1]);
//       EM::VirtualVector::VirtualReader<std::pair<K, V>> shardKvReader(
//           shardCounter[i],
//           [&](uint64_t j) { return shardReader.read().kvPair; });
//       shards[i].InitFromReaderInPlace(shardKvReader);
//     }
//   }

  // TODO: support more flexible factorization
  std::vector<uint64_t> factorizeShardCount(uint64_t shardCount) {
    switch (shardCount) {
      case 2:
        return {2};
      case 4:
        return {4};
      case 8:
        return {8};
      case 16:
        return {4, 4};
      case 32:
        return {4, 8};
      case 64:
        return {8, 8};
      case 128:
        return {4, 4, 8};
      default:
        throw std::runtime_error(
            "shardCount should be a power of 2 between 2 and 128 for init with "
            "reader");
    }
  }

  struct KVPair {
    K key;
    V value;
  };

  template <class Reader>
  void InitFromReader(Reader& reader,
                      uint64_t cacheBytes = ((uint64_t)ENCLAVE_SIZE << 20) * 3 / 4) {
    uint64_t shardCount = shards.size();
    // std::vector<uint64_t> shardCounter(shardCount);
    uint64_t initSize = reader.size();
    using NonObliviousCuckooHashMap = CuckooHashMap<K, V, false, PositionType>;
    std::vector<NonObliviousCuckooHashMap> nonOMaps(shardCount, NonObliviousCuckooHashMap(shardSize, 0));
    using Element = EM::Algorithm::TaggedT<KVPair>;

    const size_t bktSize = 8192;
    size_t bktRealSize = numRealPerBucket(bktSize, shardCount, -60);
    uint64_t minBatchSize = bktSize * shardCount;
    uint64_t maxBatchSize = cacheBytes / (sizeof(Element) + 8) / 2;
    std::vector<uint64_t> factors = factorizeShardCount(shardCount);
    if (maxBatchSize < minBatchSize) {
      throw std::runtime_error("InitFromReader cache size too small");
    }

    uint64_t totalBucketNeeded = divRoundUp(initSize, bktRealSize);
    // make sure it's a multiple of shardCount
    totalBucketNeeded = divRoundUp(totalBucketNeeded, shardCount) * shardCount;
    bktRealSize = divRoundUp(
        initSize, totalBucketNeeded);  // split initial elements evenly
    std::cout << "bktRealSize = " << bktRealSize << std::endl;
    if (totalBucketNeeded * bktSize < maxBatchSize) {
      maxBatchSize = totalBucketNeeded * bktSize;  // don't waste space
    }
    uint64_t parBatchCount = maxBatchSize / minBatchSize;
    uint64_t bucketPerBatch = parBatchCount * shardCount;
    uint64_t batchSize = bucketPerBatch * bktSize;

    // for performing mergesplits
    Element* batch = new Element[batchSize];
    Element* tempElements = new Element[batchSize];
    uint8_t* tempMarks = new uint8_t[batchSize];
    int butterflyLevel = factors.size();

    while (!reader.eof()) {
      for (uint64_t bucketIdx = 0; bucketIdx < bucketPerBatch; ++bucketIdx) {
        uint64_t bucketOffset = bucketIdx * bktSize;
        for (uint64_t i = 0; i < bktSize; ++i) {
          Element& elem = batch[bucketOffset + i];
          if (i < bktRealSize && !reader.eof()) {
            std::pair<K, V> kvPair = reader.read();
            elem.v.key = kvPair.first;
            elem.v.value = kvPair.second;
            uint64_t hash = secure_hash_with_salt((uint8_t*)&elem.v.key,
                                                  sizeof(K), randSalt);
            elem.setTag(getShardByHash(hash));
          } else {
            elem.setDummy();
          }
        }
      }

      for (int level = 0, stride = parBatchCount; level < butterflyLevel;
           ++level) {
        uint64_t way = factors[level];
        std::cout << "level = " << level << " way = " << way << std::endl;
        uint64_t parCount = bucketPerBatch / way;
#pragma omp parallel for schedule(static)
        for (uint64_t parIdx = 0; parIdx < parCount; ++parIdx) {
          uint64_t groupIdx = parIdx / stride;
          uint64_t groupOffset = parIdx % stride;
          Element* KWayIts[8];
          for (uint64_t j = 0; j < way; ++j) {
            KWayIts[j] =
                batch + ((j + groupIdx * way) * stride + groupOffset) * bktSize;
          }
          size_t tempBucketsSize = way * bktSize * sizeof(Element);
          Element* localTempElements = tempElements + parIdx * way * bktSize;
          uint8_t* localTempMarks = tempMarks + parIdx * way * bktSize;
          MergeSplitKWay(KWayIts, way, bktSize, localTempElements,
                         localTempMarks);
        }
        stride *= way;
      }
#ifndef NDEBUG
      uint64_t realCount = 0;
      for (uint64_t i = 0; i < batchSize; ++i) {
        if (!batch[i].isDummy()) {
          uint64_t hash = secure_hash_with_salt((uint8_t*)&batch[i].v.key,
                                                sizeof(K), randSalt);
          uint64_t shardIdx = getShardByHash(hash);
          Assert(shardIdx == i / (bktSize * parBatchCount));
          ++realCount;
        }
      }
      Assert(realCount == initSize);
#endif

      uint64_t shardInitSize = parBatchCount * bktSize;
  #pragma omp parallel for schedule(static)
      for (uint32_t i = 0; i < shardCount; ++i) {
        for (uint64_t j = 0; j < shardInitSize; ++j) {
          const Element& elem = batch[i * shardInitSize + j];
          const KVPair& kvPair = elem.v;
          // insert dummies like normal elements
          nonOMaps[i].template insert<true>(kvPair.key, kvPair.value, elem.isDummy());
        }
      }
    }
    delete[] tempElements;
    delete[] tempMarks;
    delete[] batch;
    for (auto& shard : shards) {
      shard.SetSize(shardSize, cacheBytes / shardCount);
    }
    #pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < shardCount; ++i) {
      shards[i].InitFromNonOblivious(nonOMaps[i]);
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
    // #pragma omp for schedule(static)
    for (uint64_t i = 0; i < batchSize; ++i) {
      uint64_t keyHash =
          secure_hash_with_salt((uint8_t*)&keyBegin[i], sizeof(K), randSalt);
      shardIndexVec[i] = getShardByHash(keyHash);
      obliMove(isDuplicate[i], shardIndexVec[i], DUMMY<uint32_t>());
    }
    return shardIndexVec;
  }

  struct RecoveryInfo {
    std::vector<uint32_t> recoveryArray;
    std::vector<bool> isDuplicate;
    RecoveryInfo() {}
    RecoveryInfo(uint32_t size) : recoveryArray(size), isDuplicate(size) {}
  };

  RecoveryInfo sortAndSuppressDuplicateKeys(std::vector<K>& keyVec,
                                            int bitonicThread = 4) {
    RecoveryInfo recoveryInfo(keyVec.size());
    for (uint32_t i = 0; i < keyVec.size(); ++i) {
      recoveryInfo.recoveryArray[i] = i;
    }
    EM::Algorithm::ParBitonicSortSepPayload(keyVec.begin(), keyVec.end(),
                                            recoveryInfo.recoveryArray.begin(),
                                            bitonicThread);

    // }
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
    uint64_t shardCount = shards.size();
    // int threadNum = omp_get_num_threads();
    // int chunkSize = divRoundUp(batchSize, threadNum);
    // #pragma omp for schedule(static, chunkSize)
    for (uint32_t i = 0; i < batchSize; ++i) {
      for (uint32_t j = 0; j < shardOutput.size(); ++j) {
        bool matchFlag = shardIndexVec[i] == j;
        obliMove(matchFlag, *(outBegin + i), shardOutput[j][i]);
      }
    }
  }

  struct KeyInfo {
    K key;
    uint64_t hash;
    uint32_t shardIdx;
    bool operator<(const KeyInfo& other) const {
      return hash < other.hash | (hash == other.hash & key < other.key);
    }
#ifndef ENCLAVE_MODE
    friend std::ostream& operator<<(std::ostream& os, const KeyInfo& keyInfo) {
      os << "key = " << keyInfo.key << " hash = " << keyInfo.hash
         << " shardIdx = " << keyInfo.shardIdx;
      return os;
    }
#endif
  };

  struct KVInfo {
    K key;
    V value;
    uint64_t hash;
    uint32_t shardIdx;
    bool operator<(const KVInfo& other) const {
      return hash < other.hash | (hash == other.hash & key < other.key);
    }
#ifndef ENCLAVE_MODE
    friend std::ostream& operator<<(std::ostream& os, const KVInfo& keyInfo) {
      os << "key = " << keyInfo.key << " hash = " << keyInfo.hash
         << " shardIdx = " << keyInfo.shardIdx << " value = " << keyInfo.value;
      return os;
    }
#endif
  };

  template <class KeyIterator, class ValueIterator>
  std::vector<uint8_t> findParBatchDeferWriteBack(const KeyIterator keyBegin,
                                                  const KeyIterator keyEnd,
                                                  ValueIterator valueBegin) {
    uint64_t shardCount = shards.size();
    uint64_t batchSize = keyEnd - keyBegin;
    uint64_t shardSize = maxQueryPerShard(batchSize, shardCount);
    std::vector<KeyInfo> keyInfoVec(batchSize);
    for (uint32_t i = 0; i < batchSize; ++i) {
      keyInfoVec[i].key = *(keyBegin + i);
      keyInfoVec[i].hash = secure_hash_with_salt((uint8_t*)&keyInfoVec[i].key,
                                                 sizeof(K), randSalt);
      keyInfoVec[i].shardIdx = getShardByHash(keyInfoVec[i].hash);
    }
    std::vector<uint32_t> recoveryArr(batchSize);
    for (uint32_t i = 0; i < batchSize; ++i) {
      recoveryArr[i] = i;
    }
    EM::Algorithm::ParBitonicSortSepPayload(keyInfoVec.begin(), keyInfoVec.end(), recoveryArr.begin(), shards.size() * 2);
    std::vector<uint32_t> shardLoads(shardCount, 0);
    std::vector<uint32_t> prefixSumFirstCompaction(batchSize + 1);
    prefixSumFirstCompaction[0] = 0;
    prefixSumFirstCompaction[1] = 1;
    for (uint32_t j = 0; j < shardCount; ++j) {
      bool matchFlag = keyInfoVec[0].shardIdx == j;
      shardLoads[j] += matchFlag;
    }
    for (uint32_t i = 1; i < batchSize; ++i) {
      bool isDup = keyInfoVec[i].key == keyInfoVec[i - 1].key;
      prefixSumFirstCompaction[i + 1] = prefixSumFirstCompaction[i] + !isDup;
#pragma omp simd
      for (uint32_t j = 0; j < shardCount; ++j) {
        bool matchFlag = (keyInfoVec[i].shardIdx == j) & !isDup;
        shardLoads[j] += matchFlag;
      }
    }
    std::vector<K> keyVec(shardCount * shardSize);
    for (uint32_t i = 0; i < batchSize; ++i) {
      keyVec[i] = keyInfoVec[i].key;
    }
    EM::Algorithm::OrCompactSeparateMark(keyVec.begin(),
                                         keyVec.begin() + batchSize,
                                         prefixSumFirstCompaction.begin());

    std::vector<uint32_t> prefixSumSecondCompaction(shardCount * shardSize + 1);
    std::vector<uint32_t> shardLoadPrefixSum(shardCount);
    shardLoadPrefixSum[0] = 0;
    for (uint32_t i = 0; i < shardCount - 1; ++i) {
      shardLoadPrefixSum[i + 1] = shardLoadPrefixSum[i] + shardLoads[i];
    }
    prefixSumSecondCompaction[0] = 0;
    for (uint32_t i = 0; i < shardCount; ++i) {
      for (uint32_t j = 0; j < shardSize; ++j) {
        uint32_t rankInShard = j + 1;
        obliMove(rankInShard > shardLoads[i], rankInShard, shardLoads[i]);
        prefixSumSecondCompaction[i * shardSize + j + 1] =
            shardLoadPrefixSum[i] + rankInShard;
      }
    }

    EM::Algorithm::OrDistributeSeparateMark(keyVec.begin(), keyVec.end(),
                                            prefixSumSecondCompaction.begin());
    using ValResult = BaseMap::ValResult;
    std::vector<ValResult> resultVec(shardCount * shardSize);

#pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < shardCount; ++i) {
      shards[i].findBatchDeferWriteBack(keyVec.begin() + i * shardSize,
                                        keyVec.begin() + (i + 1) * shardSize,
                                        resultVec.begin() + i * shardSize);
    }

    EM::Algorithm::OrCompactSeparateMark(resultVec.begin(), resultVec.end(),
                                         prefixSumSecondCompaction.begin());

    EM::Algorithm::OrDistributeSeparateMark(resultVec.begin(),
                                            resultVec.begin() + batchSize,
                                            prefixSumFirstCompaction.begin());

    for (uint32_t i = 1; i < batchSize; ++i) {
      bool isDup = keyInfoVec[i].key == keyInfoVec[i - 1].key;
      obliMove(isDup, resultVec[i], resultVec[i - 1]);
    }

    EM::Algorithm::ParBitonicSortSepPayload(recoveryArr.begin(), recoveryArr.end(), resultVec.begin(), shards.size() * 2);
    std::vector<uint8_t> foundFlags(batchSize);
    for (uint32_t i = 0; i < batchSize; ++i) {
      foundFlags[i] = resultVec[i].found;
      *(valueBegin + i) = resultVec[i].value;
    }
    return foundFlags;
  }

  void ParWriteBack() {
    uint64_t shardCount = shards.size();
#pragma omp parallel for num_threads(shardCount * 2)
    for (uint64_t i = 0; i < shardCount * 2; ++i) {
      shards[i / 2].writeBackTable(i % 2);
    }
  }

  template <class KeyIterator, class ValueIterator>
  std::vector<uint8_t> findParBatch(const KeyIterator keyBegin,
                                    const KeyIterator keyEnd,
                                    ValueIterator valueBegin) {
    const auto& res = findParBatchDeferWriteBack(keyBegin, keyEnd, valueBegin);
    ParWriteBack();
    return res;
  }

  template <class KeyIterator>
  void requireDistinctKeys(const KeyIterator keyBegin,
                           const KeyIterator keyEnd) {
    std::vector<K> keyVec(keyBegin, keyEnd);
    EM::Algorithm::BitonicSort(keyVec.begin(), keyVec.end());
    bool distinctFlag = true;
    for (uint32_t i = 1; i < keyVec.size(); ++i) {
      distinctFlag &= keyVec[i] != keyVec[i - 1];
    }
    if (!distinctFlag) {
      throw std::runtime_error("keys should be distinct");
    }
  }

template <class KeyIterator, class ValueIterator>
  std::vector<uint8_t> insertParBatch(const KeyIterator keyBegin,
                                      const KeyIterator keyEnd,
                                      const ValueIterator valueBegin) {
  uint64_t shardCount = shards.size();
    uint64_t batchSize = keyEnd - keyBegin;
    uint64_t shardSize = maxQueryPerShard(batchSize, shardCount);
    std::vector<KVInfo> keyInfoVec(batchSize);
    for (uint32_t i = 0; i < batchSize; ++i) {
      keyInfoVec[i].key = *(keyBegin + i);
      keyInfoVec[i].value = *(valueBegin + i);
      keyInfoVec[i].hash = secure_hash_with_salt((uint8_t*)&keyInfoVec[i].key,
                                                 sizeof(K), randSalt);
      keyInfoVec[i].shardIdx = getShardByHash(keyInfoVec[i].hash);
    }
    std::vector<uint32_t> recoveryArr(batchSize);
    for (uint32_t i = 0; i < batchSize; ++i) {
      recoveryArr[i] = i;
    }
    EM::Algorithm::ParBitonicSortSepPayload(keyInfoVec.begin(), keyInfoVec.end(), recoveryArr.begin(), shards.size() * 2);
    std::vector<uint32_t> shardLoads(shardCount, 0);
    std::vector<uint32_t> prefixSumFirstCompaction(batchSize + 1);
    prefixSumFirstCompaction[0] = 0;
    prefixSumFirstCompaction[1] = 1;
    for (uint32_t j = 0; j < shardCount; ++j) {
      bool matchFlag = keyInfoVec[0].shardIdx == j;
      shardLoads[j] += matchFlag;
    }
    for (uint32_t i = 1; i < batchSize; ++i) {
      bool isDup = keyInfoVec[i].key == keyInfoVec[i - 1].key;
      prefixSumFirstCompaction[i + 1] = prefixSumFirstCompaction[i] + !isDup;
#pragma omp simd
      for (uint32_t j = 0; j < shardCount; ++j) {
        bool matchFlag = (keyInfoVec[i].shardIdx == j) & !isDup;
        shardLoads[j] += matchFlag;
      }
    }
    std::vector<KVPair> kvVec(shardCount * shardSize);
    for (uint32_t i = 0; i < batchSize; ++i) {
      kvVec[i].key = keyInfoVec[i].key;
      kvVec[i].value = keyInfoVec[i].value;
    }
    EM::Algorithm::OrCompactSeparateMark(kvVec.begin(),
                                         kvVec.begin() + batchSize,
                                         prefixSumFirstCompaction.begin());

    std::vector<uint32_t> prefixSumSecondCompaction(shardCount * shardSize + 1);
    std::vector<uint32_t> shardLoadPrefixSum(shardCount);
    shardLoadPrefixSum[0] = 0;
    for (uint32_t i = 0; i < shardCount - 1; ++i) {
      shardLoadPrefixSum[i + 1] = shardLoadPrefixSum[i] + shardLoads[i];
    }
    prefixSumSecondCompaction[0] = 0;
    for (uint32_t i = 0; i < shardCount; ++i) {
      for (uint32_t j = 0; j < shardSize; ++j) {
        uint32_t rankInShard = j + 1;
        obliMove(rankInShard > shardLoads[i], rankInShard, shardLoads[i]);
        prefixSumSecondCompaction[i * shardSize + j + 1] =
            shardLoadPrefixSum[i] + rankInShard;
      }
    }

    EM::Algorithm::OrDistributeSeparateMark(kvVec.begin(), kvVec.end(),
                                            prefixSumSecondCompaction.begin());
    // using ValResult = BaseMap::ValResult;
    std::vector<uint8_t> foundVec(shardCount * shardSize);
    #pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < shardCount; ++i) {
      uint32_t numReal = shardLoads[i];
      for (uint32_t j = 0; j < shardSize; ++j) {
          foundVec[i * shardSize + j] =
              shards[i].insertOblivious(kvVec[i * shardSize + j].key, kvVec[i * shardSize + j].value, j >= numReal);
        }
    }

    EM::Algorithm::OrCompactSeparateMark(foundVec.begin(), foundVec.end(),
                                         prefixSumSecondCompaction.begin());

    EM::Algorithm::OrDistributeSeparateMark(foundVec.begin(),
                                            foundVec.begin() + batchSize,
                                            prefixSumFirstCompaction.begin());

    for (uint32_t i = 1; i < batchSize; ++i) {
      bool isDup = keyInfoVec[i].key == keyInfoVec[i - 1].key;
      obliMove(isDup, foundVec[i], foundVec[i - 1]);
    }

    EM::Algorithm::ParBitonicSortSepPayload(recoveryArr.begin(), recoveryArr.end(), foundVec.begin(), shards.size() * 2);
    foundVec.resize(batchSize);
    return foundVec;

  }

  template <class KeyIterator>
  std::vector<uint8_t> eraseParBatch(const KeyIterator keyBegin,
                                     const KeyIterator keyEnd) {
    // requireDistinctKeys(keyBegin, keyEnd);
    uint64_t shardCount = shards.size();
    uint64_t batchSize = keyEnd - keyBegin;
    uint64_t shardSize = maxQueryPerShard(batchSize, shardCount);
    std::vector<KeyInfo> keyInfoVec(batchSize);
    for (uint32_t i = 0; i < batchSize; ++i) {
      keyInfoVec[i].key = *(keyBegin + i);
      keyInfoVec[i].hash = secure_hash_with_salt((uint8_t*)&keyInfoVec[i].key,
                                                 sizeof(K), randSalt);
      keyInfoVec[i].shardIdx = getShardByHash(keyInfoVec[i].hash);
    }
    std::vector<uint32_t> recoveryArr(batchSize);
    for (uint32_t i = 0; i < batchSize; ++i) {
      recoveryArr[i] = i;
    }
    EM::Algorithm::ParBitonicSortSepPayload(keyInfoVec.begin(), keyInfoVec.end(), recoveryArr.begin(), shards.size() * 2);
    std::vector<uint32_t> shardLoads(shardCount, 0);
    std::vector<uint32_t> prefixSumFirstCompaction(batchSize + 1);
    prefixSumFirstCompaction[0] = 0;
    prefixSumFirstCompaction[1] = 1;
    for (uint32_t j = 0; j < shardCount; ++j) {
      bool matchFlag = keyInfoVec[0].shardIdx == j;
      shardLoads[j] += matchFlag;
    }
    for (uint32_t i = 1; i < batchSize; ++i) {
      bool isDup = keyInfoVec[i].key == keyInfoVec[i - 1].key;
      prefixSumFirstCompaction[i + 1] = prefixSumFirstCompaction[i] + !isDup;
#pragma omp simd
      for (uint32_t j = 0; j < shardCount; ++j) {
        bool matchFlag = (keyInfoVec[i].shardIdx == j) & !isDup;
        shardLoads[j] += matchFlag;
      }
    }
    std::vector<K> keyVec(shardCount * shardSize);
    for (uint32_t i = 0; i < batchSize; ++i) {
      keyVec[i] = keyInfoVec[i].key;
    }
    EM::Algorithm::OrCompactSeparateMark(keyVec.begin(),
                                         keyVec.begin() + batchSize,
                                         prefixSumFirstCompaction.begin());

    std::vector<uint32_t> prefixSumSecondCompaction(shardCount * shardSize + 1);
    std::vector<uint32_t> shardLoadPrefixSum(shardCount);
    shardLoadPrefixSum[0] = 0;
    for (uint32_t i = 0; i < shardCount - 1; ++i) {
      shardLoadPrefixSum[i + 1] = shardLoadPrefixSum[i] + shardLoads[i];
    }
    prefixSumSecondCompaction[0] = 0;
    for (uint32_t i = 0; i < shardCount; ++i) {
      for (uint32_t j = 0; j < shardSize; ++j) {
        uint32_t rankInShard = j + 1;
        obliMove(rankInShard > shardLoads[i], rankInShard, shardLoads[i]);
        prefixSumSecondCompaction[i * shardSize + j + 1] =
            shardLoadPrefixSum[i] + rankInShard;
      }
    }

    EM::Algorithm::OrDistributeSeparateMark(keyVec.begin(), keyVec.end(),
                                            prefixSumSecondCompaction.begin());
    std::vector<uint8_t> foundVec(shardCount * shardSize);
    #pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < shardCount; ++i) {
      uint32_t numReal = shardLoads[i];
      for (uint32_t j = 0; j < shardSize; ++j) {
          foundVec[i * shardSize + j] =
              shards[i].eraseOblivious(keyVec[i * shardSize + j], j >= numReal);
        }
    }

    EM::Algorithm::OrCompactSeparateMark(foundVec.begin(), foundVec.end(),
                                         prefixSumSecondCompaction.begin());

    EM::Algorithm::OrDistributeSeparateMark(foundVec.begin(),
                                            foundVec.begin() + batchSize,
                                            prefixSumFirstCompaction.begin());

    for (uint32_t i = 1; i < batchSize; ++i) {
      bool isDup = keyInfoVec[i].key == keyInfoVec[i - 1].key;
      obliMove(isDup, foundVec[i], foundVec[i - 1]);
    }

    EM::Algorithm::ParBitonicSortSepPayload(recoveryArr.begin(), recoveryArr.end(), foundVec.begin(), shards.size() * 2);
    foundVec.resize(batchSize);
    return foundVec;
  }
};
}  // namespace ODSL