#pragma once

#include <omp.h>

#include "algorithm/merge_split.hpp"
#include "algorithm/param_select.hpp"
#include "omap_short_kv.hpp"

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
  using BaseMap = OHashMap<K, V, FULL_OBLIVIOUS, PositionType, true>;

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
      return (hash < other.hash) | ((hash == other.hash) & (key < other.key));
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
      return (hash < other.hash) | ((hash == other.hash) & (key < other.key));
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

  // the capacity of the map
  PositionType _size = 0;

  // whether the map is initialized
  bool inited = false;

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
          Algorithm::binomLogSf(n, batchSize, 1.0 / (double)shardCount);
      return logSf < logFailProb;
    };
    return Algorithm::lowerBound(divRoundUp(batchSize, shardCount), batchSize,
                                 satisfy);
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
      double logSf = Algorithm::binomLogSf(
          bktSize, (bktSize - numDummy) * shardCount, 1.0 / shardCount);
      return logSf < logFailProb;
    };
    return bktSize - Algorithm::lowerBound(1UL, bktSize - 1, satisfy);
  }

  /**
   * @brief Factorize the number of shards into a list of factors between 2
   * and 8. We want to have as few factors as possible and the factors should be
   * close to each other. The method will throw an error if the number of shards
   * doesn't meet the requirements.
   *
   * @param shardCount The number of shards
   * @return std::vector<uint64_t> The list of factors
   */
  static std::vector<uint64_t> factorizeShardCount(uint64_t shardCount) {
    if (shardCount == 1) {
      throw std::runtime_error(
          "shardCount should be at least 2 for init with "
          "reader");
    }
    if (shardCount > 64) {
      throw std::runtime_error(
          "shardCount should be no more than 64 for init with "
          "reader");
    }
    switch (shardCount) {
      case 2:
      case 3:
      case 4:
      case 5:
      case 6:
      case 7:
      case 8:
        return {shardCount};
      case 9:
        return {3, 3};
      case 10:
        return {2, 5};
      case 12:
        return {3, 4};
      case 14:
        return {2, 7};
      case 15:
        return {3, 5};
      case 16:
        return {4, 4};
      case 18:
        return {3, 6};
      case 20:
        return {4, 5};
      case 21:
        return {3, 7};
      case 24:
        return {3, 8};
      case 25:
        return {5, 5};
      case 27:
        return {3, 3, 3};
      case 28:
        return {4, 7};
      case 30:
        return {5, 6};
      case 32:
        return {4, 8};
      case 35:
        return {5, 7};
      case 36:
        return {6, 6};
      case 40:
        return {5, 8};
      case 42:
        return {6, 7};
      case 45:
        return {3, 3, 5};
      case 48:
        return {6, 8};
      case 49:
        return {7, 7};
      case 50:
        return {2, 5, 5};
      case 54:
        return {3, 3, 6};
      case 56:
        return {7, 8};
      case 60:
        return {3, 4, 5};
      case 63:
        return {3, 3, 7};
      case 64:
        return {8, 8};
      default:
        throw std::runtime_error(
            "shardCount should be a power of 2 between 2 and 64 for init with "
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
   * @param shardCount The number of shards, should be between 2 and 64 and be a
   * factor of numbers between 2 and 8 to support Initialization. Please use
   * GetSuitableShardCount method to get a suitable shard count.
   */
  ParOMap(PositionType mapSize, uint64_t shardCount) {
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
  void SetSize(PositionType mapSize, uint64_t shardCount) {
    if (shardCount == 0) {
      throw std::runtime_error("shardCount should be positive");
    }
    if (shardCount > 64) {
      throw std::runtime_error("shardCount should be no more than 64");
    }
    _size = mapSize;
    shards.resize(shardCount);
    shardSize = (PositionType)maxQueryPerShard(mapSize, shardCount, -60);
    read_rand(randSalt, sizeof(randSalt));
    shardHashRange = (UINT64_MAX - 1) / shardCount + 1;
  }

  /**
   * @brief Get a suitable shard count given the number of threads.
   *
   * @param numThreads The number of threads available for parallelization.
   * @param emptyInit Whether the omap will be initialized with empty data.
   * @return uint32_t The suitable shard count.
   */
  static uint32_t GetSuitableShardCount(uint32_t numThreads,
                                        bool emptyInit = false) {
    if (numThreads <= 4) {
      // too few threads, use 2 shards
      return 2;
    }
    uint32_t suitableShardCount = numThreads / 2;
    if (suitableShardCount > 64) {
      suitableShardCount = 64;
    }
    if (emptyInit) {
      // don't need to factorize
      return suitableShardCount;
    }
    for (; suitableShardCount >= 1; --suitableShardCount) {
      try {
        factorizeShardCount(suitableShardCount);
        // if no exception is thrown, then the shard count is valid
        return suitableShardCount;
      } catch (const std::runtime_error& e) {
        continue;
      }
    }
    return 2;
  }

  /**
   * @brief Return the capacity of the map.
   *
   * @return PositionType
   */
  PositionType size() const { return _size; }

  /**
   * @brief Initialize an empty map with a given cache size for all shards.
   *
   * @param cacheBytes The cache size in bytes (for all shards in total).
   */
  void Init(size_t cacheBytes = MAX_CACHE_SIZE) {
    if (_size == 0) {
      throw std::runtime_error("ParOMap size not set. Call SetSize first.");
    }
    if (inited) {
      throw std::runtime_error("ParOMap double initialization");
    }
    inited = true;
#pragma omp parallel for
    for (auto& shard : shards) {
      shard.SetSize(shardSize, cacheBytes / shards.size());
      shard.Init();
    }
  }

  /**
   * @brief An object that stores the curret state of initialization.
   * Faciliates initialization in a streaming fashion. We route the data in
   * parallel through a multi-way butterfly network to load balance them to each
   * shard. Used techniques described in https://eprint.iacr.org/2023/1258.
   * Might throw error with some negligible probability (when some bucket
   * overflows), in which case one may retry.
   *
   */
  struct InitContext {
    using Element = Algorithm::TaggedT<KVPair>;
    using NonObliviousOHashMap = OHashMap<K, V, NON_OBLIVIOUS, PositionType>;
    ParOMap& omap;                                // reference to the parent map
    std::vector<NonObliviousOHashMap*> nonOMaps;  // non-oblivious maps
    std::vector<uint64_t> factors;  // number of ways in each level of butterfly
    Element* batch;                 // buffer for the current batches
    Element* tempElements;          // buffer for mergesplit
    uint8_t* tempMarks;             // buffer for mergesplit
    uint64_t shardCount;            // number of shards
    uint64_t bktSize;               // size of each bucket
    uint64_t bktPerBatch;           // number of buckets in each batch
    uint64_t bktRealSize;           // number of real elements in each bucket
    uint64_t parBatchCount;         // number of batches in parallel
    uint64_t batchSize;          // size of each batch (in number of elements)
    uint64_t perBatchShardSize;  // size of each shard in each batch
    uint64_t cacheBytes;         // cache size available
    uint64_t currBktIdx;         // current bucket index in the current batch
    uint64_t currBktOffset;      // current offset in the current bucket
    uint64_t load;       // number of elements inserted in the current batch
    int butterflyLevel;  // number of levels in the butterfly network

    /**
     * @brief Route elements in the current batch through the butterfly network,
     * and insert into the non oblivious map.
     *
     */
    void route() {
      // route the data through the butterfly network in parallel
      bool initFailed = false;
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
          try {
            MergeSplitKWay(KWayIts, way, bktSize, localTempElements,
                           localTempMarks);
          } catch (const std::runtime_error& e) {
            initFailed = true;
          }
        }
        stride *= way;
      }
      if (initFailed) {
        throw std::runtime_error(
            "Initialization failed due to bucket overflow. Please retry.");
      }
      // insert the data into the non-oblivious maps in parallel
#pragma omp parallel for schedule(static)
      for (uint32_t i = 0; i < shardCount; ++i) {
        for (uint64_t j = 0; j < perBatchShardSize; ++j) {
          const Element& elem = batch[i * perBatchShardSize + j];
          const KVPair& kvPair = elem.v;
          // insert dummies like normal elements
          bool existFlag = nonOMaps[i]->template Insert<true>(
              kvPair.key, kvPair.value, elem.IsDummy());
          if (existFlag & !elem.IsDummy()) {
            initFailed = true;
            break;
          }
        }
      }
      if (initFailed) {
        throw std::runtime_error(
            "Initialization failed. Encountered duplicate keys.");
      }
    }

    /**
     * @brief Construct a new InitContext object
     *
     * @param omap The parent map
     * @param initSize The number of elements to be inserted
     * @param cacheBytes The cache size for all shards
     */
    InitContext(ParOMap& omap, PositionType initSize,
                uint64_t cacheBytes = MAX_CACHE_SIZE)
        : omap(omap),
          cacheBytes(cacheBytes),
          currBktIdx(0),
          currBktOffset(0),
          load(0) {
      if (omap.size() == 0) {
        throw std::runtime_error("ParOMap size not set. Call SetSize first.");
      }
      if (omap.inited) {
        throw std::runtime_error("ParOMap double initialization");
      }
      shardCount = omap.shards.size();
      uint64_t maxInitSizePerShard =
          maxQueryPerShard(initSize, shardCount, -60);
      using NonObliviousOHashMap = OHashMap<K, V, NON_OBLIVIOUS, PositionType>;
      nonOMaps.resize(shardCount);

      using Element = Algorithm::TaggedT<KVPair>;
      bktSize = std::min(8192UL, GetNextPowerOfTwo(maxInitSizePerShard));
      bktRealSize = numRealPerBkt(bktSize, shardCount, -60);
      uint64_t minBatchSize = bktSize * shardCount;
      uint64_t maxBatchSize = cacheBytes / (sizeof(Element) + 8) / 2;
      factors = omap.factorizeShardCount(shardCount);
      if (maxBatchSize < minBatchSize) {
        throw std::runtime_error("InitFromReader cache size too small");
      }

      uint64_t totalBktNeeded = divRoundUp(initSize, bktRealSize);
      // make sure it's a multiple of shardCount
      uint64_t bktPerShard = divRoundUp(totalBktNeeded, shardCount);
      if (bktPerShard * bktSize > omap.shardSize) {
        omap.shardSize =
            bktPerShard *
            bktSize;  // we need large oram to hold these many elements
      }
      for (auto& nonOMapPtr : nonOMaps) {
        nonOMapPtr = new NonObliviousOHashMap();
        nonOMapPtr->SetSize(omap.shardSize, 0);
      }
      totalBktNeeded = bktPerShard * shardCount;
      bktRealSize = divRoundUp(
          initSize, totalBktNeeded);  // split initial elements evenly
      if (totalBktNeeded * bktSize < maxBatchSize) {
        maxBatchSize = totalBktNeeded * bktSize;  // don't waste space
      }
      parBatchCount = maxBatchSize / minBatchSize;
      bktPerBatch = parBatchCount * shardCount;
      batchSize = bktPerBatch * bktSize;
      perBatchShardSize = parBatchCount * bktSize;

      // buffers for mergesplit
      batch = new Element[batchSize];
      tempElements = new Element[batchSize];
      tempMarks = new uint8_t[batchSize];
      butterflyLevel = factors.size();
    }

    InitContext(const InitContext&) = delete;
    InitContext(InitContext&&) = default;

    /**
     * @brief Insert a new key-value pair into the batch. The key must be
     * unique. The running time of this function can be unstable.
     * The method may throw error with some negligible probability (when some
     * bucket overflows), in which case one may retry the initialization.
     *
     * @param key
     * @param value
     */
    void Insert(const K& key, const V& value) {
      if (load >= omap.size()) {
        throw std::runtime_error("Too many elements inserted during init");
      }
      ++load;
      Element& elem = batch[currBktIdx * bktSize + currBktOffset];
      elem.v.key = key;
      elem.v.value = value;
      uint64_t hash = secure_hash_with_salt((uint8_t*)&elem.v.key, sizeof(K),
                                            omap.randSalt);
      elem.setTag(omap.getShardByHash(hash));
      if (++currBktOffset == bktRealSize) {
        uint64_t baseOffset = currBktIdx * bktSize;
        for (; currBktOffset < bktSize; ++currBktOffset) {
          batch[baseOffset + currBktOffset].setDummy();
        }
        currBktOffset = 0;
        if (++currBktIdx == bktPerBatch) {
          route();
          currBktIdx = 0;
        }
      }
    }

    void Insert(std::pair<K, V> kvPair) { Insert(kvPair.first, kvPair.second); }

    template <class Iterator>
    void InsertBatch(Iterator begin, Iterator end) {
      for (auto it = begin; it != end; ++it) {
        Insert(it->first, it->second);
      }
    }

    /**
     * @brief Finished inserting elements. Close the context and initialize
     * the map. The method may throw error with some negligible probability,
     * in which case one may retry the initialization.
     *
     */
    void Finalize() {
      // set the rest of buckets to dummy
      uint64_t batchOffset = currBktIdx * bktSize + currBktOffset;
      if (batchOffset != 0) {
        for (uint64_t i = batchOffset; i < batchSize; ++i) {
          batch[i].setDummy();
        }
        route();
      }
      delete[] tempElements;
      delete[] tempMarks;
      delete[] batch;
      tempElements = NULL;
      tempMarks = NULL;
      batch = NULL;
      for (auto& shard : omap.shards) {
        shard.SetSize(omap.shardSize, cacheBytes / shardCount);
      }
// initialize the oblivious maps from the non-oblivious maps in parallel
#pragma omp parallel for schedule(static)
      for (uint32_t i = 0; i < shardCount; ++i) {
        omap.shards[i].InitFromNonOblivious(*nonOMaps[i]);
        delete nonOMaps[i];
        nonOMaps[i] = NULL;
      }
    }

    ~InitContext() {
      if (tempElements) {
        delete[] tempElements;
      }
      if (tempMarks) {
        delete[] tempMarks;
      }
      if (batch) {
        delete[] batch;
      }
      for (auto& nonOMapPtr : nonOMaps) {
        if (nonOMapPtr) {
          delete nonOMapPtr;
        }
      }
    }
  };

  /**
   * @brief Obtain a new context to initialize this map. The initialization
   data
   * can be either private or public (the initialization and subsequent
   accesses
   * are oblivious). The data needs to be unique.
   * Example:
   *  auto* initContext = parOMap.NewInitContext(kvMap.size(), 1UL << 28);
      for (auto it = kvMap.begin(); it != kvMap.end(); ++it;) {
        initContext->Insert(it->first, it->second);
      }
      initContext->Finalize();
      delete initContext;
   *
   * @param initSize The number of elements to be inserted (could be an
   estimate, but should be no less than the actual number of elements)
   * @param cacheBytes The cache size for all shards.
   * @return InitContext The context object
   */
  InitContext* NewInitContext(PositionType initSize,
                              uint64_t cacheBytes = MAX_CACHE_SIZE) {
    return new InitContext(*this, initSize, cacheBytes);
  }

  /**
   * @brief Initialize the map from a reader. The data is fetched in batches
   * and for each batch, we route the data in parallel through a multi-way
   * butterfly network to load balance them to each shard. Used techniques
   * described in https://eprint.iacr.org/2023/1258. Might throw error with
   * some negligible probability (when some bucket overflows), in which case
   * one may retry.
   *
   * @tparam Reader The type of the reader
   * @param reader The reader object
   * @param cacheBytes The cache size for all shards
   */
  template <class Reader>
    requires Readable<Reader, std::pair<K, V>>
  void InitFromReader(Reader& reader, uint64_t cacheBytes = MAX_CACHE_SIZE) {
    InitContext context(*this, reader.size(), cacheBytes);
    while (!reader.eof()) {
      std::pair<K, V> kvPair = reader.read();
      context.Insert(kvPair.first, kvPair.second);
    }
    context.Finalize();
  }

  /**
   * @brief Find the values associated with a batch of keys. The keys may
   * contain duplicates and may be arranged in any order. The values are
   * written in the same order as the keys. Larger batch size will likely
   * increase the throughput. The writeback process of the internal orams are
   * deferred until WriteBack is called. The method may leak information with
   * negligible probability, which happens when the keys are very unbalanced.
   * Since the keys each shard receives are oblivious, an adversary cannot
   * force such an unbalanced distribution. This method is ~2x faster than the
   * algorithm described in the Snoopy's paper. Note: the method is not
   * thread-safe, and calling it in parallel can cause deadlock.
   *
   * @tparam KeyIterator The type of the key iterator
   * @tparam ValueIterator The type of the value iterator
   * @param keyBegin The beginning of the keys
   * @param keyEnd The end of the keys
   * @param valueBegin The beginning of the output values
   * @return std::vector<uint8_t> A vector of flags indicating whether each
   * key is found
   */
  template <class KeyIterator, class ValueIterator>
  std::vector<uint8_t> FindBatchDeferWriteBack(const KeyIterator keyBegin,
                                               const KeyIterator keyEnd,
                                               ValueIterator valueBegin) {
    uint64_t shardCount = shards.size();
    if (keyEnd - keyBegin > UINT32_MAX) {
      // it takes too long to load balance a super large batch, consider
      // breaking it into smaller ones
      throw std::runtime_error("FindBatchDeferWriteBack: batch size too large");
    }
    uint32_t batchSize = (uint32_t)(keyEnd - keyBegin);
    // Calculate the max unique element per shard with high probability.
    uint32_t batchSizePerShard =
        (uint32_t)maxQueryPerShard(batchSize, shardCount);
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
    Algorithm::ParBitonicSortSepPayload(keyInfoVec.begin(), keyInfoVec.end(),
                                        recoveryArr.begin(),
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
        bool matchFlag = (keyInfoVec[i].shardIdx == j) & (!isDup);
        shardLoads[j] += matchFlag;
      }
    }
    // With negligible probability, we need to enlarge the shard size for the
    // batch to avoid overflow
    for (uint32_t load : shardLoads) {
      obliMove(load > batchSizePerShard, batchSizePerShard, load);
    }

    std::vector<K> keyVec(shardCount * batchSizePerShard);
    for (uint32_t i = 0; i < batchSize; ++i) {
      keyVec[i] = keyInfoVec[i].key;
    }
    // Compact the unique elements in the sorted array using the prefix
    // sum.
    Algorithm::OrCompactSeparateMark(keyVec.begin(), keyVec.begin() + batchSize,
                                     prefixSumFirstCompaction.begin());
    std::vector<uint32_t> prefixSumSecondCompaction(
        shardCount * batchSizePerShard + 1);
    std::vector<uint32_t> shardLoadPrefixSum(shardCount);
    shardLoadPrefixSum[0] = 0;
    for (uint32_t i = 0; i < shardCount - 1; ++i) {
      shardLoadPrefixSum[i + 1] = shardLoadPrefixSum[i] + shardLoads[i];
    }
    prefixSumSecondCompaction[0] = 0;
    for (uint32_t i = 0; i < shardCount; ++i) {
      for (uint32_t j = 0; j < batchSizePerShard; ++j) {
        uint32_t rankInShard = j + 1;
        obliMove(rankInShard > shardLoads[i], rankInShard, shardLoads[i]);
        prefixSumSecondCompaction[i * batchSizePerShard + j + 1] =
            shardLoadPrefixSum[i] + rankInShard;
      }
    }
    // Distribute the compacted elements so that the elements of shard i
    // starts at offset shard_size * i.
    Algorithm::OrDistributeSeparateMark(keyVec.begin(), keyVec.end(),
                                        prefixSumSecondCompaction.begin());
    using ValResult = BaseMap::ValResult;
    std::vector<ValResult> resultVec(shardCount * batchSizePerShard);

    // Parallel query each shard, and get the result values (along with flags
    // of existence).
#pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < shardCount; ++i) {
      shards[i].FindBatchDeferWriteBack(
          keyVec.begin() + i * batchSizePerShard,
          keyVec.begin() + (i + 1) * batchSizePerShard,
          resultVec.begin() + i * batchSizePerShard);
    }

    // Compact the result values in reverse order of the previous
    // distribution.
    Algorithm::OrCompactSeparateMark(resultVec.begin(), resultVec.end(),
                                     prefixSumSecondCompaction.begin());

    // Distribute the compacted result values in reverse order of the previous
    // compaction.
    Algorithm::OrDistributeSeparateMark(resultVec.begin(),
                                        resultVec.begin() + batchSize,
                                        prefixSumFirstCompaction.begin());

    // Propagate the values of the duplicate elements.
    for (uint32_t i = 1; i < batchSize; ++i) {
      bool isDup = keyInfoVec[i].key == keyInfoVec[i - 1].key;
      obliMove(isDup, resultVec[i], resultVec[i - 1]);
    }

    // Parallel bitonic sort the values in reverse order the previous bitonic
    // sort using the recovery array
    Algorithm::ParBitonicSortSepPayload(recoveryArr.begin(), recoveryArr.end(),
                                        resultVec.begin(),
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
  void WriteBack() {
    uint64_t shardCount = shards.size();
#pragma omp parallel for num_threads(shardCount * 2)
    for (uint64_t i = 0; i < shardCount * 2; ++i) {
      shards[i / 2].WriteBackTable(i % 2);
    }
  }

  /**
   * @brief Combines FindBatchDeferWriteBack and WriteBack.
   *
   * @tparam KeyIterator The type of the key iterator
   * @tparam ValueIterator The type of the value iterator
   * @param keyBegin The beginning of the keys
   * @param keyEnd The end of the keys
   * @param valueBegin The beginning of the output values
   * @return std::vector<uint8_t> A vector of flags indicating whether each
   * key is found
   */
  template <class KeyIterator, class ValueIterator>
  std::vector<uint8_t> FindBatch(const KeyIterator keyBegin,
                                 const KeyIterator keyEnd,
                                 ValueIterator valueBegin) {
    const auto& res = FindBatchDeferWriteBack(keyBegin, keyEnd, valueBegin);
    WriteBack();
    return res;
  }

  /**
   * @brief Insert a batch of key-value pairs into the map in parallel. For
   * duplicate keys, an arbitrary value will be inserted. The method may leak
   * information with negligible probability, similar to
   * FindBatchDeferWriteBack. The method is not thread-safe.
   *
   * @tparam KeyIterator The type of the key iterator
   * @tparam ValueIterator The type of the value iterator
   * @param keyBegin The beginning of the keys
   * @param keyEnd The end of the keys
   * @param valueBegin The beginning of the input values
   * @return std::vector<uint8_t> A vector of flags indicating whether each
   * key already exists.
   */
  template <class KeyIterator, class ValueIterator>
  std::vector<uint8_t> InsertBatch(const KeyIterator keyBegin,
                                   const KeyIterator keyEnd,
                                   const ValueIterator valueBegin) {
    // the overall procedure is similar to FindBatchDeferWriteBack
    uint64_t shardCount = shards.size();
    if (keyEnd - keyBegin > UINT32_MAX) {
      // it takes too long to load balance a super large batch, consider
      // breaking it into smaller ones
      throw std::runtime_error("InsertBatch: batch size too large");
    }
    uint32_t batchSize = (uint32_t)(keyEnd - keyBegin);
    uint32_t batchSizePerShard =
        (uint32_t)maxQueryPerShard(batchSize, shardCount);
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
    Algorithm::ParBitonicSortSepPayload(keyInfoVec.begin(), keyInfoVec.end(),
                                        recoveryArr.begin(),
                                        (int)shards.size() * 2);
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
        bool matchFlag = (keyInfoVec[i].shardIdx == j) & (!isDup);
        shardLoads[j] += matchFlag;
      }
    }
    // With negligible probability, we need to enlarge the shard size for the
    // batch to avoid overflow
    for (uint32_t load : shardLoads) {
      obliMove(load > batchSizePerShard, batchSizePerShard, load);
    }

    std::vector<KVPair> kvVec(shardCount * batchSizePerShard);
    for (uint32_t i = 0; i < batchSize; ++i) {
      kvVec[i].key = keyInfoVec[i].key;
      kvVec[i].value = keyInfoVec[i].value;
    }
    Algorithm::OrCompactSeparateMark(kvVec.begin(), kvVec.begin() + batchSize,
                                     prefixSumFirstCompaction.begin());

    std::vector<uint32_t> prefixSumSecondCompaction(
        shardCount * batchSizePerShard + 1);
    std::vector<uint32_t> shardLoadPrefixSum(shardCount);
    shardLoadPrefixSum[0] = 0;
    for (uint32_t i = 0; i < shardCount - 1; ++i) {
      shardLoadPrefixSum[i + 1] = shardLoadPrefixSum[i] + shardLoads[i];
    }
    prefixSumSecondCompaction[0] = 0;
    for (uint32_t i = 0; i < shardCount; ++i) {
      for (uint32_t j = 0; j < batchSizePerShard; ++j) {
        uint32_t rankInShard = j + 1;
        obliMove(rankInShard > shardLoads[i], rankInShard, shardLoads[i]);
        prefixSumSecondCompaction[i * batchSizePerShard + j + 1] =
            shardLoadPrefixSum[i] + rankInShard;
      }
    }

    Algorithm::OrDistributeSeparateMark(kvVec.begin(), kvVec.end(),
                                        prefixSumSecondCompaction.begin());
    std::vector<uint8_t> foundVec(shardCount * batchSizePerShard);
#pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < shardCount; ++i) {
      uint32_t numReal = shardLoads[i];
      for (uint32_t j = 0; j < batchSizePerShard; ++j) {
        foundVec[i * batchSizePerShard + j] = shards[i].OInsert(
            kvVec[i * batchSizePerShard + j].key,
            kvVec[i * batchSizePerShard + j].value, j >= numReal);
      }
    }

    Algorithm::OrCompactSeparateMark(foundVec.begin(), foundVec.end(),
                                     prefixSumSecondCompaction.begin());

    Algorithm::OrDistributeSeparateMark(foundVec.begin(),
                                        foundVec.begin() + batchSize,
                                        prefixSumFirstCompaction.begin());

    for (uint32_t i = 1; i < batchSize; ++i) {
      bool isDup = keyInfoVec[i].key == keyInfoVec[i - 1].key;
      obliMove(isDup, foundVec[i], foundVec[i - 1]);
    }

    Algorithm::ParBitonicSortSepPayload(recoveryArr.begin(), recoveryArr.end(),
                                        foundVec.begin(),
                                        (int)shards.size() * 2);
    foundVec.resize(batchSize);
    return foundVec;
  }

  /**
   * @brief Erase a batch of keys from the map in parallel. Keys may contain
   * duplicate and be arranged in any order. The method may leak information
   * with negligible probability, similar to FindBatchDeferWriteBack. The
   * method is not thread-safe.
   *
   * @tparam KeyIterator The type of the key iterator
   * @param keyBegin The beginning of the keys to erase
   * @param keyEnd The end of the keys to erase
   * @return std::vector<uint8_t> A vector of flags indicating whether each
   * key existed in the map.
   */
  template <class KeyIterator>
  std::vector<uint8_t> EraseBatch(const KeyIterator keyBegin,
                                  const KeyIterator keyEnd) {
    if (keyEnd - keyBegin > UINT32_MAX) {
      // it takes too long to load balance a super large batch, consider
      // breaking it into smaller ones
      throw std::runtime_error("EraseBatch: batch size too large");
    }
    uint64_t shardCount = shards.size();
    uint32_t batchSize = (uint32_t)(keyEnd - keyBegin);
    uint32_t batchSizePerShard =
        (uint32_t)maxQueryPerShard(batchSize, shardCount);
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
    Algorithm::ParBitonicSortSepPayload(keyInfoVec.begin(), keyInfoVec.end(),
                                        recoveryArr.begin(), shards.size() * 2);
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
        bool matchFlag = (keyInfoVec[i].shardIdx == j) & (!isDup);
        shardLoads[j] += matchFlag;
      }
    }
    // With negligible probability, we need to enlarge the shard size for the
    // batch to avoid overflow
    for (uint32_t load : shardLoads) {
      obliMove(load > batchSizePerShard, batchSizePerShard, load);
    }

    std::vector<K> keyVec(shardCount * batchSizePerShard);
    for (uint32_t i = 0; i < batchSize; ++i) {
      keyVec[i] = keyInfoVec[i].key;
    }
    Algorithm::OrCompactSeparateMark(keyVec.begin(), keyVec.begin() + batchSize,
                                     prefixSumFirstCompaction.begin());

    std::vector<uint32_t> prefixSumSecondCompaction(
        shardCount * batchSizePerShard + 1);
    std::vector<uint32_t> shardLoadPrefixSum(shardCount);
    shardLoadPrefixSum[0] = 0;
    for (uint32_t i = 0; i < shardCount - 1; ++i) {
      shardLoadPrefixSum[i + 1] = shardLoadPrefixSum[i] + shardLoads[i];
    }
    prefixSumSecondCompaction[0] = 0;
    for (uint32_t i = 0; i < shardCount; ++i) {
      for (uint32_t j = 0; j < batchSizePerShard; ++j) {
        uint32_t rankInShard = j + 1;
        obliMove(rankInShard > shardLoads[i], rankInShard, shardLoads[i]);
        prefixSumSecondCompaction[i * batchSizePerShard + j + 1] =
            shardLoadPrefixSum[i] + rankInShard;
      }
    }

    Algorithm::OrDistributeSeparateMark(keyVec.begin(), keyVec.end(),
                                        prefixSumSecondCompaction.begin());
    std::vector<uint8_t> foundVec(shardCount * batchSizePerShard);
#pragma omp parallel for schedule(static)
    for (uint32_t i = 0; i < shardCount; ++i) {
      uint32_t numReal = shardLoads[i];
      for (uint32_t j = 0; j < batchSizePerShard; ++j) {
        foundVec[i * batchSizePerShard + j] =
            shards[i].OErase(keyVec[i * batchSizePerShard + j], j >= numReal);
      }
    }

    Algorithm::OrCompactSeparateMark(foundVec.begin(), foundVec.end(),
                                     prefixSumSecondCompaction.begin());

    Algorithm::OrDistributeSeparateMark(foundVec.begin(),
                                        foundVec.begin() + batchSize,
                                        prefixSumFirstCompaction.begin());

    for (uint32_t i = 1; i < batchSize; ++i) {
      bool isDup = keyInfoVec[i].key == keyInfoVec[i - 1].key;
      obliMove(isDup, foundVec[i], foundVec[i - 1]);
    }

    Algorithm::ParBitonicSortSepPayload(recoveryArr.begin(), recoveryArr.end(),
                                        foundVec.begin(), shards.size() * 2);
    foundVec.resize(batchSize);
    return foundVec;
  }
};
}  // namespace ODSL