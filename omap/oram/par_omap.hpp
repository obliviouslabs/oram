#pragma once

#include <omp.h>

#include "external_memory/algorithm/merge_split.hpp"
#include "external_memory/algorithm/param_select.hpp"
#include "omap.hpp"

/// @brief A parallel oblivious map by sharding. Each shard is an oblivious map,
/// and we load balance a batch of queries to each shard obliviously.
/// Used the idea of https://eprint.iacr.org/2021/1280 but optimized the
/// algorithm.

namespace ODSL {
/**
 * @brief A parallel oblivious map by sharding. Each shard is an oblivious map,
 * and we load balance a batch of queries to each shard obliviously.
 *
 * @tparam K key type
 * @tparam V value type
 * @tparam PositionType The type of position used in the oblivious map, default
 * to uint64_t.
 */
template <typename K, typename V, typename PositionType = uint64_t>
struct ParOMap {
 private:
  using BaseMap = OHashMap<K, V, true, PositionType, true>;

  struct KVPair {
    K key;
    V value;
  };

  // data structure for sorting a batch of insertion requests
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

  // data structure for sorting a batch of queries / erasures
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

  // underlying shards
  std::vector<BaseMap> shards;
  // the size of each shard
  PositionType shardSize = 0;
  // random salt for hashing keys to shards
  uint8_t randSalt[16];

  // the size of hash range for each shard
  uint64_t shardHashRange;

  /**
   * @brief translate a hash to a shard index
   *
   * @param hash the hash value
   * @return uint32_t the shard index
   */
  INLINE uint32_t getShardByHash(uint64_t hash) {
    return (uint32_t)(hash / shardHashRange);
  }

  /**
   * @brief Calculate an upper bound of the number of queries each shard can
   * take, given an allowed failure probability.
   *
   * @param batchSize The number of queries in the batch
   * @param shardCount The number of shards
   * @param logFailProb Logarithm of the allowed failure probability
   * @return uint64_t Upper bound of the number of queries each shard can
   * take
   */
  static uint64_t maxQueryPerShard(uint64_t batchSize, uint64_t shardCount,
                                   double logFailProb = -40) {
    auto satisfy = [&](uint64_t n) {
      double logSf =
          EM::Algorithm::binomLogSf(n, batchSize, 1.0 / (double)shardCount);
      return logSf < logFailProb;
    };
    return EM::Algorithm::lowerBound(divRoundUp(batchSize, shardCount),
                                     batchSize, satisfy);
  }

  /**
   * @brief Calculate the maximum number of real elements each bkt can hold
   * initially when routing data through the butterfly network during the
   * initialization of the oblivious maps.
   *
   * @param bktSize The size of each bkt
   * @param shardCount The number of shards
   * @param logFailProb Logarithm of the allowed failure probability
   * @return uint64_t The maximum number of real elements each bkt may hold
   */
  static uint64_t numRealPerBkt(uint64_t bktSize, uint64_t shardCount,
                                double logFailProb = -60) {
    auto satisfy = [&](uint64_t numDummy) {
      double logSf = EM::Algorithm::binomLogSf(
          bktSize, (bktSize - numDummy) * shardCount, 1.0 / shardCount);
      return logSf < logFailProb;
    };
    return bktSize - EM::Algorithm::lowerBound(1UL, bktSize - 1, satisfy);
  }

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

 public:
  ParOMap() {}

  /**
   * @brief Construct a new parallel oblivious map object, does not allocate
   * resources. Must call Init or InitFromReader before using the map.
   *
   * @param mapSize The capacity of the map
   * @param shardCount The number of shards, should be a power of 2 to support
   * InitFromReader. Suggested to be the number of available threads / 2.
   */
  ParOMap(uint64_t mapSize, uint64_t shardCount) {
    SetSize(mapSize, shardCount);
  }

  /**
   * @brief Set the size of the map and the number of shards, does not allocate
   * resources.
   *
   * @param mapSize The capacity of the map
   * @param shardCount The number of shards, should be a power of 2 to support
   * InitFromReader. Suggested to be the number of available threads / 2.
   */
  void SetSize(uint64_t mapSize, uint64_t shardCount) {
    if (shardCount == 0) {
      throw std::runtime_error("shardCount should be positive");
    }
    if (shardCount > 128) {
      throw std::runtime_error("shardCount should be no more than 128");
    }
    shards.resize(shardCount);
    shardSize = (PositionType)maxQueryPerShard(mapSize, shardCount, -60);
    read_rand(randSalt, sizeof(randSalt));
    shardHashRange = (UINT64_MAX - 1) / shardCount + 1;
  }

  /**
   * @brief Initialize an empty map with a given cache size for all shards.
   *
   * @param cacheBytes The cache size in bytes (for all shards in total).
   */
  void Init(size_t cacheBytes = DEFAULT_HEAP_SIZE) {
#pragma omp parallel for
    for (auto& shard : shards) {
      shard.SetSize(shardSize, cacheBytes / shards.size());
      shard.Init();
    }
  }

  /**
   * @brief Initialize the map from a reader. The data is fetched in batches and
   * for each batch, we route the data in parallel through a multi-way butterfly
   * network to load balance them to each shard. Used techniques described in
   * https://eprint.iacr.org/2023/1258.
   *
   * @tparam Reader The type of the reader
   * @param reader The reader object
   * @param cacheBytes
   */
  template <class Reader>
  void InitFromReader(Reader& reader, uint64_t cacheBytes = DEFAULT_HEAP_SIZE) {
    uint64_t shardCount = shards.size();
    uint64_t initSize = reader.size();
    uint64_t maxInitSizePerShard = maxQueryPerShard(initSize, shardCount, -60);
    using NonObliviousOHashMap = OHashMap<K, V, false, PositionType>;
    std::vector<NonObliviousOHashMap> nonOMaps(shardCount);
    for (auto& nonOMap : nonOMaps) {
      nonOMap.SetSize(shardSize, 0);
    }

    using Element = EM::Algorithm::TaggedT<KVPair>;
    const size_t bktSize =
        std::min(8192UL, GetNextPowerOfTwo(maxInitSizePerShard));
    size_t bktRealSize = numRealPerBkt(bktSize, shardCount, -60);
    uint64_t minBatchSize = bktSize * shardCount;
    uint64_t maxBatchSize = cacheBytes / (sizeof(Element) + 8) / 2;
    std::vector<uint64_t> factors = factorizeShardCount(shardCount);
    if (maxBatchSize < minBatchSize) {
      throw std::runtime_error("InitFromReader cache size too small");
    }

    uint64_t totalBktNeeded = divRoundUp(initSize, bktRealSize);
    // make sure it's a multiple of shardCount
    uint64_t bktPerShard = divRoundUp(totalBktNeeded, shardCount);
    if (bktPerShard * bktSize > shardSize) {
      shardSize = bktPerShard *
                  bktSize;  // we need large oram to hold these many elements
    }
    totalBktNeeded = bktPerShard * shardCount;
    bktRealSize =
        divRoundUp(initSize, totalBktNeeded);  // split initial elements evenly
    if (totalBktNeeded * bktSize < maxBatchSize) {
      maxBatchSize = totalBktNeeded * bktSize;  // don't waste space
    }
    uint64_t parBatchCount = maxBatchSize / minBatchSize;
    uint64_t bktPerBatch = parBatchCount * shardCount;
    uint64_t batchSize = bktPerBatch * bktSize;
    uint64_t perBatchShardSize = parBatchCount * bktSize;

    // buffers for mergesplit
    Element* batch = new Element[batchSize];
    Element* tempElements = new Element[batchSize];
    uint8_t* tempMarks = new uint8_t[batchSize];
    int butterflyLevel = factors.size();

    while (!reader.eof()) {
      // set the first level of the butterfly network
      for (uint64_t bktIdx = 0; bktIdx < bktPerBatch; ++bktIdx) {
        uint64_t bktOffset = bktIdx * bktSize;
        for (uint64_t i = 0; i < bktSize; ++i) {
          Element& elem = batch[bktOffset + i];
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

      // route the data through the butterfly network in parallel
      for (int level = 0, stride = parBatchCount; level < butterflyLevel;
           ++level) {
        uint64_t way = factors[level];
        uint64_t parCount = bktPerBatch / way;
#pragma omp parallel for schedule(static)
        for (uint64_t parIdx = 0; parIdx < parCount; ++parIdx) {
          uint64_t groupIdx = parIdx / stride;
          uint64_t groupOffset = parIdx % stride;
          Element* KWayIts[8];
          for (uint64_t j = 0; j < way; ++j) {
            KWayIts[j] =
                batch + ((j + groupIdx * way) * stride + groupOffset) * bktSize;
          }
          size_t tempBktsSize = way * bktSize * sizeof(Element);
          Element* localTempElements = tempElements + parIdx * way * bktSize;
          uint8_t* localTempMarks = tempMarks + parIdx * way * bktSize;
          MergeSplitKWay(KWayIts, way, bktSize, localTempElements,
                         localTempMarks);
        }
        stride *= way;
      }

      // check if the data is routed correctly in debug mode
#ifndef NDEBUG
      uint64_t realCount = 0;
      for (uint64_t i = 0; i < batchSize; ++i) {
        if (!batch[i].IsDummy()) {
          uint64_t hash = secure_hash_with_salt((uint8_t*)&batch[i].v.key,
                                                sizeof(K), randSalt);
          uint64_t shardIdx = getShardByHash(hash);
          Assert(shardIdx == i / (bktSize * parBatchCount));
          ++realCount;
        }
      }
      Assert(realCount == initSize);
#endif

      // insert the data into the non-oblivious maps in parallel
#pragma omp parallel for schedule(static)
      for (uint32_t i = 0; i < shardCount; ++i) {
        for (uint64_t j = 0; j < perBatchShardSize; ++j) {
          const Element& elem = batch[i * perBatchShardSize + j];
          const KVPair& kvPair = elem.v;
          // insert dummies like normal elements
          nonOMaps[i].template Insert<true>(kvPair.key, kvPair.value,
                                            elem.IsDummy());
        }
      }
    }
    delete[] tempElements;
    delete[] tempMarks;
    delete[] batch;
    for (auto& shard : shards) {
      shard.SetSize(shardSize, cacheBytes / shardCount);
    }
    // initialize the oblivious maps from the non-oblivious maps in parallel
#pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < shardCount; ++i) {
      shards[i].InitFromNonOblivious(nonOMaps[i]);
    }
  }

  /**
   * @brief Find the values associated with a batch of keys. The keys may
   * contain duplicates and may be arranged in any order. The values are written
   * in the same order as the keys. Larger batch size will likely increase the
   * throughput. The writeback process of the internal orams are deferred until
   * ParWriteBack is called. The method may leak information with negligible
   * probability, which happens when the keys are very unbalanced. Since the
   * keys each shard receives are oblivious, an adversary cannot force such an
   * unbalanced distribution. This method is ~2x faster than the algorithm
   * described in the Snoopy's paper. Note: the method is not thread-safe, and
   * calling it in parallel can cause deadlock.
   *
   * @tparam KeyIterator The type of the key iterator
   * @tparam ValueIterator The type of the value iterator
   * @param keyBegin The beginning of the keys
   * @param keyEnd The end of the keys
   * @param valueBegin The beginning of the output values
   * @return std::vector<uint8_t> A vector of flags indicating whether each key
   * is found
   */
  template <class KeyIterator, class ValueIterator>
  std::vector<uint8_t> FindParBatchDeferWriteBack(const KeyIterator keyBegin,
                                                  const KeyIterator keyEnd,
                                                  ValueIterator valueBegin) {
    uint64_t shardCount = shards.size();
    uint64_t batchSize = keyEnd - keyBegin;
    // Calculate the max unique element per shard with high probability.
    uint64_t shardSize = maxQueryPerShard(batchSize, shardCount);
    std::vector<KeyInfo> keyInfoVec(batchSize);
    // hash the keys and determine the shard index for each key
    for (uint32_t i = 0; i < batchSize; ++i) {
      keyInfoVec[i].key = *(keyBegin + i);
      keyInfoVec[i].hash = secure_hash_with_salt((uint8_t*)&keyInfoVec[i].key,
                                                 sizeof(K), randSalt);
      keyInfoVec[i].shardIdx = getShardByHash(keyInfoVec[i].hash);
    }
    // Parallel bitonic sort the input batch first by their hash, and if
    // there's a collision, by the key (for obliviousness, always perform two
    // comparisons). Sort together an index array as payload for recovery.
    std::vector<uint32_t> recoveryArr(batchSize);
    for (uint32_t i = 0; i < batchSize; ++i) {
      recoveryArr[i] = i;
    }
    EM::Algorithm::ParBitonicSortSepPayload(
        keyInfoVec.begin(), keyInfoVec.end(), recoveryArr.begin(),
        (int)shards.size() * 2);

    // Obliviously count the number of unique keys for each shard (using
    // simd acceleration). Also generate a prefix sum of the number of unique
    // elements.
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
    // Compact the unique elements in the sorted array using the prefix
    // sum.
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
    // Distribute the compacted elements so that the elements of shard i starts
    // at offset shard_size * i.
    EM::Algorithm::OrDistributeSeparateMark(keyVec.begin(), keyVec.end(),
                                            prefixSumSecondCompaction.begin());
    using ValResult = BaseMap::ValResult;
    std::vector<ValResult> resultVec(shardCount * shardSize);

    // Parallel query each shard, and get the result values (along with flags of
    // existence).
#pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < shardCount; ++i) {
      shards[i].FindBatchDeferWriteBack(keyVec.begin() + i * shardSize,
                                        keyVec.begin() + (i + 1) * shardSize,
                                        resultVec.begin() + i * shardSize);
    }

    // Compact the result values in reverse order of the previous distribution.
    EM::Algorithm::OrCompactSeparateMark(resultVec.begin(), resultVec.end(),
                                         prefixSumSecondCompaction.begin());

    // Distribute the compacted result values in reverse order of the previous
    // compaction.
    EM::Algorithm::OrDistributeSeparateMark(resultVec.begin(),
                                            resultVec.begin() + batchSize,
                                            prefixSumFirstCompaction.begin());

    // Propagate the values of the duplicate elements.
    for (uint32_t i = 1; i < batchSize; ++i) {
      bool isDup = keyInfoVec[i].key == keyInfoVec[i - 1].key;
      obliMove(isDup, resultVec[i], resultVec[i - 1]);
    }

    // Parallel bitonic sort the values in reverse order the previous bitonic
    // sort using the recovery array
    EM::Algorithm::ParBitonicSortSepPayload(
        recoveryArr.begin(), recoveryArr.end(), resultVec.begin(),
        (int)shards.size() * 2);
    std::vector<uint8_t> foundFlags(batchSize);
    for (uint32_t i = 0; i < batchSize; ++i) {
      foundFlags[i] = resultVec[i].found;
      *(valueBegin + i) = resultVec[i].value;
    }
    return foundFlags;
  }

  /**
   * @brief Perform the deferred writeback of the internal orams in parallel,
   * and release the locks.
   *
   */
  void ParWriteBack() {
    uint64_t shardCount = shards.size();
#pragma omp parallel for num_threads(shardCount * 2)
    for (uint64_t i = 0; i < shardCount * 2; ++i) {
      shards[i / 2].WriteBackTable(i % 2);
    }
  }

  /**
   * @brief Combines FindParBatchDeferWriteBack and ParWriteBack.
   *
   * @tparam KeyIterator The type of the key iterator
   * @tparam ValueIterator The type of the value iterator
   * @param keyBegin The beginning of the keys
   * @param keyEnd The end of the keys
   * @param valueBegin The beginning of the output values
   * @return std::vector<uint8_t> A vector of flags indicating whether each key
   * is found
   */
  template <class KeyIterator, class ValueIterator>
  std::vector<uint8_t> FindParBatch(const KeyIterator keyBegin,
                                    const KeyIterator keyEnd,
                                    ValueIterator valueBegin) {
    const auto& res = FindParBatchDeferWriteBack(keyBegin, keyEnd, valueBegin);
    ParWriteBack();
    return res;
  }

  /**
   * @brief Insert a batch of key-value pairs into the map in parallel. For
   * duplicate keys, an arbitrary value will be inserted. The method may leak
   * information with negligible probability, similar to
   * FindParBatchDeferWriteBack. The method is not thread-safe.
   *
   * @tparam KeyIterator The type of the key iterator
   * @tparam ValueIterator The type of the value iterator
   * @param keyBegin The beginning of the keys
   * @param keyEnd The end of the keys
   * @param valueBegin The beginning of the input values
   * @return std::vector<uint8_t> A vector of flags indicating whether each key
   * already exists.
   */
  template <class KeyIterator, class ValueIterator>
  std::vector<uint8_t> InsertParBatch(const KeyIterator keyBegin,
                                      const KeyIterator keyEnd,
                                      const ValueIterator valueBegin) {
    // the overall procedure is similar to FindParBatchDeferWriteBack
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
    EM::Algorithm::ParBitonicSortSepPayload(
        keyInfoVec.begin(), keyInfoVec.end(), recoveryArr.begin(),
        shards.size() * 2);
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
    std::vector<uint8_t> foundVec(shardCount * shardSize);
#pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < shardCount; ++i) {
      uint32_t numReal = shardLoads[i];
      for (uint32_t j = 0; j < shardSize; ++j) {
        foundVec[i * shardSize + j] = shards[i].InsertOblivious(
            kvVec[i * shardSize + j].key, kvVec[i * shardSize + j].value,
            j >= numReal);
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

    EM::Algorithm::ParBitonicSortSepPayload(recoveryArr.begin(),
                                            recoveryArr.end(), foundVec.begin(),
                                            shards.size() * 2);
    foundVec.resize(batchSize);
    return foundVec;
  }

  /**
   * @brief Erase a batch of keys from the map in parallel. Keys may contain
   * duplicate and be arranged in any order. The method may leak information
   * with negligible probability, similar to FindParBatchDeferWriteBack. The
   * method is not thread-safe.
   *
   * @tparam KeyIterator The type of the key iterator
   * @param keyBegin The beginning of the keys to erase
   * @param keyEnd The end of the keys to erase
   * @return std::vector<uint8_t> A vector of flags indicating whether each key
   * existed in the map.
   */
  template <class KeyIterator>
  std::vector<uint8_t> EraseParBatch(const KeyIterator keyBegin,
                                     const KeyIterator keyEnd) {
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
    EM::Algorithm::ParBitonicSortSepPayload(
        keyInfoVec.begin(), keyInfoVec.end(), recoveryArr.begin(),
        shards.size() * 2);
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
            shards[i].EraseOblivious(keyVec[i * shardSize + j], j >= numReal);
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

    EM::Algorithm::ParBitonicSortSepPayload(recoveryArr.begin(),
                                            recoveryArr.end(), foundVec.begin(),
                                            shards.size() * 2);
    foundVec.resize(batchSize);
    return foundVec;
  }
};
}  // namespace ODSL