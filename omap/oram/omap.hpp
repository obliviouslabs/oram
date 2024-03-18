#pragma once
#include "recursive_oram.hpp"

/// @brief This file implements a cuckoo hash map built on top of the
/// recursive ORAM. Specifically, it contains two hash tables implemented with
/// recursive oram. Each key is hashed to two positions, one for each table, and
/// at each position of a table is a bucket with two slots, the element can be
/// stored in any of these four slots. In addition, there is a stash for
/// elements that cannot be stored in the tables.

namespace ODSL {

template <typename K, typename PositionType = uint64_t>

/**
 * @brief Used to hash the key to two positions
 *
 */
struct OHashMapIndexer {
  static constexpr int saltLength = 16;  // 128 bits salt
  uint8_t salts[saltLength];             // salt for secure hash
  PositionType _size;                    // number of buckets in each hash table

  OHashMapIndexer() {}

  OHashMapIndexer(PositionType size) { SetSize(size); }

  /**
   * @brief Initialize the indexer
   *
   * @param size number of buckets in each hash table
   */
  void SetSize(PositionType size) {
    _size = size;
    read_rand(salts, saltLength);
  }

  struct HashIndices {
    uint64_t h0;
    uint64_t h1;
  };

  /**
   * @brief Get both hash indices
   *
   * @param key the key to hash
   * @param pos0 the first hash index
   * @param pos1 the second hash index
   */
  void getHashIndices(const K& key, PositionType& pos0,
                      PositionType& pos1) const {
    HashIndices hashIndices;
    secure_hash_with_salt(key, salts, (uint8_t*)&hashIndices,
                          sizeof(hashIndices));
    pos0 = hashIndices.h0 % _size;
    pos1 = hashIndices.h1 % _size;
  }

  /**
   * @brief Get hash index 0
   *
   * @param key the key to hash
   * @return PositionType the first hash index
   */
  PositionType getHashIdx0(const K& key) const {
    PositionType pos0, pos1;
    getHashIndices(key, pos0, pos1);
    return pos0;
  }

  /**
   * @brief Get hash index 1
   *
   * @param key the key to hash
   * @return PositionType the second hash index
   */
  PositionType getHashIdx1(const K& key) const {
    PositionType pos0, pos1;
    getHashIndices(key, pos0, pos1);
    return pos1;
  }
};

/**
 * @brief A single entry (slot) in the hash map. It contains a key and a value,
 * and two flags.
 *
 * @tparam K the key type
 * @tparam V the value type
 */
template <typename K, typename V>
struct OHashMapEntry {
  // whether the entry is valid. An invalid entry is an empty slot.
  bool valid = false;
  // whether the entry is a dummy. A dummy entry is used when the hash map is
  // not built on an ORAM, but we still want to hide whether an insert is real
  // or dummy. In this case, we set valid flag to true, but the dummy flag to
  // true.
  bool dummy = false;
  K key;
  V value;
#ifndef ENCLAVE_MODE
  // cout
  friend std::ostream& operator<<(std::ostream& os,
                                  const OHashMapEntry& entry) {
    os << "(" << entry.valid << ", " << entry.key << ", " << entry.value << ")";
    return os;
  }
#endif
};

/**
 * @brief A bucket can contain multiple entries that are fully associative. This
 * can significantly increase the load factor of the hash map.
 *
 * @tparam K the key type
 * @tparam V the value type
 * @tparam bucketSize the number of slots in the bucket
 */
template <typename K, typename V, const short bucketSize>
struct OHashMapBucket {
  OHashMapEntry<K, V> entries[bucketSize];
#ifndef ENCLAVE_MODE
  // cout
  friend std::ostream& operator<<(std::ostream& os,
                                  const OHashMapBucket& bucket) {
    os << "[";
    for (int i = 0; i < bucketSize; ++i) {
      os << bucket.entries[i];
      if (i != bucketSize - 1) {
        os << ", ";
      }
    }
    os << "]";
    return os;
  }
#endif
};

/**
 * @brief An cuckoo hash map built on top of either recursive ORAM or a standard
 * vector. It contains two hash tables, and a stash for elements that cannot be
 * stored in the tables.
 *
 * @tparam K the key type
 * @tparam V the value type
 * @tparam isOblivious whether the hash map is built on an ORAM
 * @tparam PositionType the type of the position
 * @tparam parallel_init whether to initialize the two hash tables in parallel
 */
template <typename K, typename V, const bool isOblivious = true,
          typename PositionType = uint64_t, const bool parallel_init = true>
struct OHashMap {
 private:
  // the ratio between the capacity and the total number of slots in
  // the hash tables
  static constexpr double loadFactor = 0.7;
  // number of slots in each bucket
  static constexpr short bucketSize = 2;
  // maximum number of elements in the stash
  static constexpr int stash_max_size = 10;
  // the capcity of the hash map
  PositionType _size;
  // the number of elements in the hash map
  PositionType load;
  // the size of each hash table
  PositionType tableSize;
  using BucketType = OHashMapBucket<K, V, bucketSize>;
  // for oblivious hash map, we use recursive ORAM
  using ObliviousTableType = RecursiveORAM<BucketType, PositionType>;
  // for non-oblivious hash map, we cache the front of the vector
  using NonObliviousTableType =
      EM::CacheFrontVector::Vector<BucketType, sizeof(BucketType), true, true,
                                   1024>;
  using TableType = std::conditional_t<isOblivious, ObliviousTableType,
                                       NonObliviousTableType>;
  TableType table0, table1;
  OHashMapIndexer<K, PositionType> indexer;
  using KVEntry = OHashMapEntry<K, V>;
  std::vector<KVEntry> stash;

  /**
   * @brief Helper function that accesses the hash table.
   *
   * @param addr the address to access
   * @param table the table to access
   * @param updateFunc the function to call on the accessed bucket
   */
  static void updateHelper(PositionType addr, TableType& table,
                           const auto& updateFunc) {
    if constexpr (isOblivious) {
      table.Access(addr, updateFunc);
    } else {
      updateFunc(table[addr]);
    }
  }

  /**
   * @brief Helper function that reads the hash table.
   *
   * @param addr the address to read
   * @param table the table to read
   * @param outputBucket the output bucket
   */
  static void readHelper(PositionType addr, TableType& table,
                         BucketType& outputBucket) {
    if constexpr (isOblivious) {
      table.Read(addr, outputBucket);
    } else {
      outputBucket = table[addr];
    }
  }

  /**
   * @brief A helper function that replaces an entry in a search range if it
   * exists. The function is not oblivious.
   *
   * @tparam hideDummy whether to hide the dummy flag of the inserted entry
   * and the existing entries in search range.
   * @param bucket the bucket to perform replace
   * @param entryToInsert the entry to insert
   * @return true if the entry is replaced, false otherwise
   *
   */
  template <const bool hideDummy = false, typename EntryIterator>
  static bool replaceIfExist(EntryIterator begin, EntryIterator end,
                             const KVEntry& entryToInsert) {
    for (auto it = begin; it != end; ++it) {
      auto& entry = *it;
      if constexpr (hideDummy) {
        if (entry.valid && ((entry.key == entryToInsert.key) &
                            (entryToInsert.dummy == entry.dummy))) {
          entry = entryToInsert;
          return true;
        }
      } else {
        if (entry.valid && entry.key == entryToInsert.key) {
          entry = entryToInsert;
          return true;
        }
      }
    }
    return false;
  }

  /**
   * @brief A helper function that replaces an entry in a bucket if it
   * exists. The function is not oblivious.
   *
   * @tparam hideDummy whether to hide the dummy flag of the inserted entry
   * and the existing entries in the bucket.
   * @param bucket the bucket to perform replace
   * @param entryToInsert the entry to insert
   * @return true if the entry is replaced, false otherwise
   *
   */
  template <const bool hideDummy = false>
  static bool replaceIfExist(BucketType& bucket, const KVEntry& entryToInsert) {
    return replaceIfExist<hideDummy>(
        bucket.entries, bucket.entries + bucketSize, entryToInsert);
  }

  /**
   * @brief A helper function that replaces an entry in the stash if it exists.
   * The function is not oblivious.
   *
   * @tparam hideDummy whether to hide the dummy flag of the inserted entry and
   * the existing entries in the stash.
   * @param bucket the bucket to perform replace
   * @param entryToInsert the entry to insert
   * @return true if the entry is replaced, false otherwise
   *
   */
  template <const bool hideDummy = false>
  static bool replaceIfExist(std::vector<KVEntry>& stash,
                             const KVEntry& entryToInsert) {
    return replaceIfExist<hideDummy>(stash.begin(), stash.end(), entryToInsert);
  }

  /**
   * @brief A helper function that inserts an entry to a bucket if a slot is
   * empty. The function is not oblivious.
   *
   * @param bucket the bucket to perform insert
   * @param entryToInsert the entry to insert
   * @return true if the entry is inserted, false otherwise
   */
  static bool insertIfEmpty(BucketType& bucket, const KVEntry& entryToInsert) {
    for (int i = 0; i < bucketSize; ++i) {
      auto& entry = bucket.entries[i];
      if (!entry.valid) {
        entry = entryToInsert;
        return true;
      }
    }
    return false;
  }

  /**
   * @brief A helper function that inserts an entry to the stash if a slot is
   * empty. The function is not oblivious.
   *
   * @param stash the stash to perform insert
   * @param entryToInsert the entry to insert
   * @return true if the entry is inserted, false otherwise
   */
  static bool insertToStash(const KVEntry& entryToInsert,
                            std::vector<KVEntry>& stash) {
    if (stash.size() < stash_max_size) {
      stash.push_back(entryToInsert);
      return true;
    }
    return false;
  }

  /**
   * @brief A helper function that replaces an entry in a search range if it
   * exists. The function is oblivious.
   *
   * @param bucket the bucket to perform replace
   * @param entryToInsert the entry to insert
   * @return true if the entry is replaced, false otherwise
   *
   */
  template <typename EntryIterator>
  static bool replaceIfExistOblivious(EntryIterator begin, EntryIterator end,
                                      const KVEntry& entryToInsert) {
    bool updated = !entryToInsert.valid;
    for (auto it = begin; it != end; ++it) {
      auto& entry = *it;
      bool matchFlag = entry.valid & entry.key == entryToInsert.key & !updated;
      obliMove(matchFlag, entry, entryToInsert);
      updated |= matchFlag;
    }
    return updated;
  }

  /**
   * @brief A helper function that replaces an entry in a bucket if it exists.
   * The function is oblivious.
   *
   * @param bucket the bucket to perform replace
   * @param entryToInsert the entry to insert
   * @return true if the entry is replaced, false otherwise
   */
  static bool replaceIfExistOblivious(BucketType& bucket,
                                      const KVEntry& entryToInsert) {
    return replaceIfExistOblivious(bucket.entries, bucket.entries + bucketSize,
                                   entryToInsert);
  }

  /**
   * @brief A helper function that replaces an entry in the stash if it exists.
   * The function is oblivious.
   *
   * @param bucket the bucket to perform replace
   * @param entryToInsert the entry to insert
   * @return true if the entry is replaced, false otherwise
   */
  static bool replaceIfExistOblivious(std::vector<KVEntry>& stash,
                                      const KVEntry& entryToInsert) {
    return replaceIfExistOblivious(stash.begin(), stash.end(), entryToInsert);
  }

  /**
   * @brief A helper function that inserts an entry to a bucket if a slot is
   * empty. The function is oblivious.
   *
   * @param bucket the bucket to perform insert
   * @param entryToInsert the entry to insert
   * @return true if the entry is inserted, false otherwise
   */
  static bool insertIfEmptyOblivious(BucketType& bucket,
                                     const KVEntry& entryToInsert) {
    bool updated = !entryToInsert.valid;
    for (int i = 0; i < bucketSize; ++i) {
      auto& entry = bucket.entries[i];
      bool isEmpty = !entry.valid;
      bool insertFlag = isEmpty & !updated;
      obliMove(insertFlag, entry, entryToInsert);
      updated |= insertFlag;
    }
    return updated;
  }

  /**
   * @brief A helper function that inserts an entry to the stash if a slot is
   * empty. The function is oblivious.
   *
   * @param entryToInsert the entry to insert
   * @param stash the stash to perform insert
   * @return true if the entry is inserted, false otherwise
   */
  static bool insertToStashOblivious(const KVEntry& entryToInsert,
                                     std::vector<KVEntry>& stash) {
    if (stash.size() < stash_max_size) {
      stash.resize(stash_max_size);
    }
    bool inserted = !entryToInsert.valid;
    for (KVEntry& entry : stash) {
      bool emptyFlag = !entry.valid;
      obliMove(emptyFlag & !inserted, entry, entryToInsert);
      inserted |= emptyFlag;
    }
    return inserted;
  }

  /**
   * @brief Helper function that searches a bucket for a key. Depending on
   * whether the hash map is oblivious, the function may reveal the number of
   * comparisons
   *
   * @param key the key to search
   * @param value the value to return
   * @param bucket the bucket to search
   * @return true if the key is found, false otherwise
   */
  static bool searchBucket(const K& key, V& value, BucketType& bucket) {
    if constexpr (isOblivious) {
      bool found = false;
      for (int i = 0; i < bucketSize; ++i) {
        const auto& entry = bucket.entries[i];
        bool match = entry.valid & entry.key == key;
        obliMove(match, value, entry.value);
        found |= match;
      }
      return found;
    } else {
      for (int i = 0; i < bucketSize; ++i) {
        const auto& entry = bucket.entries[i];
        if (entry.valid && entry.key == key) {
          value = entry.value;
          return true;
        }
      }
      return false;
    }
  }

 public:
  // result data structure for a query in a batch
  struct ValResult {
    V value;
    bool found;
  };

  OHashMap() {}

  /**
   * @brief Construct a new OHashMap object
   *
   * @param size capacity of the hash map
   */
  OHashMap(PositionType size) { SetSize(size); }

  /**
   * @brief Construct a new OHashMap object
   *
   * @param size capacity of the hash map
   * @param cacheBytes the size of available cache in bytes
   */
  OHashMap(PositionType size, uint64_t cacheBytes) {
    SetSize(size, cacheBytes);
  }

  /**
   * @brief Allocate resources for the hash map.
   *
   * @param size the capacity of the hash map
   */
  void SetSize(PositionType size) { SetSize(size, DEFAULT_HEAP_SIZE); }

  /**
   * @brief Allocate resources for the hash map.
   *
   * @param size
   * @param size the capacity of the hash map
   * @param cacheBytes the size of available cache in bytes
   */
  void SetSize(PositionType size, uint64_t cacheBytes) {
    _size = size;
    load = 0;
    tableSize = size / (2 * loadFactor * bucketSize) + 1;
    indexer.SetSize(tableSize);
    if constexpr (isOblivious) {
      // equally divide the cache between the two tables
      table0.SetSize(tableSize, cacheBytes / 2);
      table1.SetSize(tableSize, cacheBytes / 2);
    } else {
      // first give table0 all the cache, and give table1 the rest
      // because table0 is more likely to be accessed
      uint64_t minCacheBytes = TableType::GetMinMemoryUsage();
      if (cacheBytes < 2 * minCacheBytes) {
        cacheBytes = 2 * minCacheBytes;
      }
      table0.SetSizeByCacheBytes(tableSize, cacheBytes - minCacheBytes);
      if (cacheBytes - minCacheBytes < table0.GetMemoryUsage()) {
        throw std::runtime_error("Cache size too small");
      }
      uint64_t remainingCacheBytes = cacheBytes - table0.GetMemoryUsage();
      table1.SetSize(tableSize, remainingCacheBytes);
    }
  }

  /**
   * @brief Get the memory usage of this hash map
   *
   * @return uint64_t the heap memory usage in bytes
   */
  uint64_t GetMemoryUsage() const {
    return table0.GetMemoryUsage() + table1.GetMemoryUsage() +
           stash.size() * sizeof(KVEntry);
  }

  PositionType size() const { return _size; }

  PositionType GetLoad() const { return load; }

  PositionType GetTableSize() const { return tableSize; }

  const OHashMapIndexer<K, PositionType>& GetIndexer() const { return indexer; }

  const TableType& GetTable0() const { return table0; }

  const TableType& GetTable1() const { return table1; }

  const std::vector<KVEntry>& GetStash() const { return stash; }

  using NonObliviousOHashMap = OHashMap<K, V, false, PositionType>;

  /**
   * @brief Initialize the oblivious hash map from another non-oblivious hash
   * map
   *
   * @param other the non-oblivious hash map
   */
  void InitFromNonOblivious(NonObliviousOHashMap& other) {
    static_assert(isOblivious,
                  "Only oblivious hash map can call this function");
    if (_size != other.size()) {
      throw std::runtime_error("OHashMap InitFromNonOblivious failed");
    }
    load = other.GetLoad();
    indexer = other.GetIndexer();
    // change dummies to invalid
    EM::VirtualVector::VirtualReader<BucketType> reader0(
        other.GetTableSize(), [&](PositionType i) {
          BucketType bucket = other.GetTable0()[i];
          for (int j = 0; j < bucketSize; ++j) {
            bucket.entries[j].valid &= !bucket.entries[j].dummy;
            bucket.entries[j].dummy = false;
          }
          return bucket;
        });
    EM::VirtualVector::VirtualReader<BucketType> reader1(
        other.GetTableSize(), [&](PositionType i) {
          BucketType bucket = other.GetTable1()[i];
          for (int j = 0; j < bucketSize; ++j) {
            bucket.entries[j].valid &= !bucket.entries[j].dummy;
            bucket.entries[j].dummy = false;
          }
          return bucket;
        });
    if constexpr (parallel_init) {
#pragma omp task
      { table0.InitFromReader(reader0); }

      { table1.InitFromReader(reader1); }
#pragma omp taskwait
    } else {
      table0.InitFromReader(reader0);
      table1.InitFromReader(reader1);
    }
    for (const auto& entry : other.GetStash()) {
      if (entry.valid) {
        // valid is public but dummy is not
        // since we deleted all dummies, it's likely that we don't need to put
        // these entries in the stash
        insertOblivious(entry.key, entry.value, entry.dummy);
      }
    }
  }

  /**
   * @brief Initialize the hash map from a reader of key value pairs. If all the
   * keys are distinct, the operation is oblivious. Otherwise, it reveals
   * information about the data it reads, but does not affect the future
   * queries.
   *
   * @tparam Reader the type of the reader
   * @param reader A reader that reads std::pair<K, V>
   * @param additionalCacheBytes the size of additional available cache in
   * bytes
   */
  template <typename Reader>
  void InitFromReader(Reader& reader, uint64_t additionalCacheBytes = 0) {
    if constexpr (!isOblivious) {
      while (!reader.eof()) {
        std::pair<K, V> entry = reader.read();
        insert(entry.first, entry.second);
      }
    } else {
      additionalCacheBytes = std::max(
          additionalCacheBytes, NonObliviousTableType::GetMinMemoryUsage() * 2);
      NonObliviousOHashMap nonObliviousOHashMap(_size, additionalCacheBytes);
      nonObliviousOHashMap.InitFromReader(reader);
      InitFromNonOblivious(nonObliviousOHashMap);
    }
  }

  /**
   * @brief Initialize an empty hash map
   *
   */
  void InitDefault() {
    if constexpr (isOblivious) {
      if constexpr (parallel_init) {
#pragma omp task
        { table0.InitDefault(BucketType()); }

        { table1.InitDefault(BucketType()); }
#pragma omp taskwait
      } else {
        table0.InitDefault(BucketType());
        table1.InitDefault(BucketType());
      }
    }
  }

  /**
   * @brief Initialize an empty hash map
   *
   */
  void Init() { InitDefault(); }

  /**
   * @brief Non-oblivious insert, but may hide whether the insertion is dummy.
   *
   * @tparam hideDummy whether to hide the dummy flag of the inserted entry
   * @param key the key to insert
   * @param value the value to insert
   * @param isDummy whether the insertion is dummy
   * @return true if the key already exists, false otherwise
   */
  template <bool hideDummy = false>
  bool insert(const K& key, const V& value, bool isDummy = false) {
    static_assert(!(isOblivious && hideDummy),
                  "hideDummy is only useful for non oblivious");
    if constexpr (!hideDummy) {
      if (isDummy) {
        return false;
      }
    }

    KVEntry entryToInsert = {true, isDummy, key, value};
    if constexpr (hideDummy) {
      // generate a random key for dummy entry
      K randKey;
      read_rand((uint8_t*)&randKey, sizeof(K));
      obliMove(isDummy, entryToInsert.key, randKey);
    }
    PositionType idx0 = indexer.getHashIdx0(entryToInsert.key);

    bool inserted = false;
    bool exist = false;

    auto table0UpdateFunc = [&](BucketType& bucket0) {
      if (replaceIfExist<hideDummy>(bucket0, entryToInsert)) {
        inserted = true;
        exist = true;
        // found in table 0
        return;
      }
      PositionType idx1 = indexer.getHashIdx1(entryToInsert.key);
      auto table1UpdateFunc = [&](BucketType& bucket1) {
        if (replaceIfExist<hideDummy>(bucket1, entryToInsert)) {
          inserted = true;
          exist = true;
          // found in table 1
          return;
        }
        if (replaceIfExist<hideDummy>(stash, entryToInsert)) {
          exist = true;
          inserted = true;
          // found in stash
        }
        if (insertIfEmpty(bucket0, entryToInsert)) {
          inserted = true;
          // successfully inserted into table 0
          return;
        }
        if (insertIfEmpty(bucket1, entryToInsert)) {
          // successfully inserted into table 1
          inserted = true;
          return;
        }
        // replace existing entry
        std::swap(entryToInsert, bucket0.entries[0]);
      };
      updateHelper(idx1, table1, table1UpdateFunc);
    };
    updateHelper(idx0, table0, table0UpdateFunc);

    load += !exist;
    if (inserted) {
      return exist;
    }
    // the offset of the entry we intend to swap
    int offset;
    auto swapUpdateFunc = [&](BucketType& bucket) {
      if (insertIfEmpty(bucket, entryToInsert)) {
        inserted = true;
        return;
      }
      std::swap(entryToInsert, bucket.entries[offset]);
    };
    for (int r = 0; r < 15; ++r) {
      offset = (r + 1) % bucketSize;  // double check if it works
      PositionType idx1 = indexer.getHashIdx1(entryToInsert.key);
      updateHelper(idx1, table1, swapUpdateFunc);
      if (inserted) {
        return false;
      }
      PositionType idx0 = indexer.getHashIdx0(entryToInsert.key);
      updateHelper(idx0, table0, swapUpdateFunc);
      if (inserted) {
        return false;
      }
    }
    printf("Warning OHashMap uses stash\n");
    if (!insertToStash(entryToInsert, stash)) {
      throw std::runtime_error("OHashMap insert failed");
    }
    return false;
  }

  /**
   * @brief Insert obliviously. Hide the number of replacement and whether the
   * insertion is dummy
   *
   * @param key the key to insert
   * @param value the value to insert
   * @param isDummy whether the insertion is dummy
   * @return true if the key already exists, false otherwise
   */
  bool insertOblivious(const K& key, const V& value, bool isDummy = false) {
    PositionType idx0, idx1;
    indexer.getHashIndices(key, idx0, idx1);
    obliMove(isDummy, idx0, (PositionType)UniformRandom(tableSize - 1));
    obliMove(isDummy, idx1, (PositionType)UniformRandom(tableSize - 1));
    bool exist = false;
    KVEntry entryToInsert = {!isDummy, false, key, value};
    auto table0UpdateFunc = [&](BucketType& bucket0) {
      bool replaceSucceed = replaceIfExistOblivious(bucket0, entryToInsert);
      exist |= replaceSucceed;
      entryToInsert.valid &= !replaceSucceed;
      auto table1UpdateFunc = [&](BucketType& bucket1) {
        bool replaceSucceed = replaceIfExistOblivious(bucket1, entryToInsert);
        exist |= replaceSucceed;
        entryToInsert.valid &= !replaceSucceed;
        replaceSucceed = replaceIfExistOblivious(stash, entryToInsert);
        exist |= replaceSucceed;
        entryToInsert.valid &= !replaceSucceed;
        bool insertSucceed = insertIfEmptyOblivious(bucket0, entryToInsert);
        entryToInsert.valid &= !insertSucceed;
        insertSucceed = insertIfEmptyOblivious(bucket1, entryToInsert);
        entryToInsert.valid &= !insertSucceed;
        obliSwap(entryToInsert.valid, bucket0.entries[0], entryToInsert);
      };
      updateHelper(idx1, table1, table1UpdateFunc);
    };
    updateHelper(idx0, table0, table0UpdateFunc);

    load += !exist;
    int offset;
    auto swapUpdateFunc = [&](BucketType& bucket) {
      bool insertSucceed = insertIfEmptyOblivious(bucket, entryToInsert);
      entryToInsert.valid &= !insertSucceed;
      obliSwap(entryToInsert.valid, entryToInsert, bucket.entries[offset]);
    };
    for (int r = 0; r < 10; ++r) {
      offset = (r + 1) % bucketSize;  // double check if it works
      PositionType idx1 = indexer.getHashIdx1(entryToInsert.key);
      updateHelper(idx1, table1, swapUpdateFunc);  // modifies entryToInsert

      PositionType idx0 = indexer.getHashIdx0(entryToInsert.key);
      updateHelper(idx0, table0, swapUpdateFunc);
    }
    if (!insertToStashOblivious(entryToInsert, stash)) {
      throw std::runtime_error("OHashMap insert failed");
    }
    return exist;
  }

  /**
   * @brief Helper function that searches the stash for a key. Depending on
   * whether the hash map is oblivious, the function may reveal the number of
   * comparisons.
   *
   * @param key the key to search
   * @param value the value to return
   * @param stash the stash to search
   * @return true if the key is found, false otherwise
   */
  static bool searchStash(const K& key, V& value,
                          const std::vector<KVEntry>& stash) {
    if constexpr (isOblivious) {
      bool found = false;
      for (const auto& entry : stash) {
        bool match = entry.valid & entry.key == key;
        obliMove(match, value, entry.value);
        found |= match;
      }
      return found;
    } else {
      for (const auto& entry : stash) {
        if (entry.valid && entry.key == key) {
          value = entry.value;
          return true;
        }
      }
      return false;
    }
  }

  /**
   * @brief Find the value of a key. Depending on whether the hash map is
   * oblivious, the function may reveal the number of comparisons.
   *
   * @param key the key to search
   * @param value the value to return
   * @param isDummy whether the search is dummy
   * @return true if the key is found, false otherwise
   */
  bool find(const K& key, V& value, bool isDummy = false) {
    if constexpr (!isOblivious) {
      if (isDummy) {
        return false;
      }
      PositionType idx0 = indexer.getHashIdx0(key);
      BucketType bucket;
      readHelper(idx0, table0, bucket);
      bool found = searchBucket(key, value, bucket);
      if (found) {
        return true;
      }
      PositionType idx1 = indexer.getHashIdx1(key);
      readHelper(idx1, table1, bucket);
      found = searchBucket(key, value, bucket);
      if (found) {
        return true;
      }
      found = searchStash(key, value, stash);
      return found & !isDummy;
    } else {
      bool found = false;
      PositionType idx0, idx1;
      indexer.getHashIndices(key, idx0, idx1);

      obliMove(isDummy, idx0, (PositionType)UniformRandom(tableSize - 1));

      BucketType bucket;
      readHelper(idx0, table0, bucket);
      found = searchBucket(key, value, bucket);

      obliMove(isDummy, idx1, (PositionType)UniformRandom(tableSize - 1));
      readHelper(idx1, table1, bucket);
      found |= searchBucket(key, value, bucket);
      found |= searchStash(key, value, stash);
      return found & !isDummy;
    }
  }

  /**
   * @brief Find the value of a key in table 0 obliviously.
   *
   * @param key the key to search
   * @param value the value to return
   * @param isDummy whether the search is dummy
   * @return true if the key is found, false otherwise
   */
  bool findTable0(const K& key, V& value, bool isDummy = false) {
    static_assert(isOblivious);
    bool found = false;

    PositionType idx0 = indexer.getHashIdx0(key);

    obliMove(isDummy, idx0, (PositionType)UniformRandom(tableSize - 1));

    BucketType bucket;
    readHelper(idx0, table0, bucket);
    found = searchBucket(key, value, bucket);
    found |= searchStash(key, value, stash);

    return found & !isDummy;
  }

  /**
   * @brief Find the value of a key in table 1 obliviously.
   *
   * @param key the key to search
   * @param value the value to return
   * @param isDummy whether the search is dummy
   * @return true if the key is found, false otherwise
   */
  bool findTable1(const K& key, V& value, bool isDummy = false) {
    static_assert(isOblivious);
    bool found = false;

    PositionType idx1 = indexer.getHashIdx1(key);
    obliMove(isDummy, idx1, (PositionType)UniformRandom(tableSize - 1));
    BucketType bucket;
    readHelper(idx1, table1, bucket);
    found = searchBucket(key, value, bucket);

    return found & !isDummy;
  }

  /**
   * @brief Find the values of a batch of keys obliviously in one table, and
   * defer the writeback procedure of the recursive ORAM. Requires the hash map
   * to be oblivious. The input may contain duplicate keys and may be arranged
   * in any order.
   *
   * @tparam KeyIter the type of the key iterator
   * @tparam ValResIter the type of the value result iterator
   * @param tableNum the table to search, 0 or 1
   * @param keyBegin the begin iterator of the keys
   * @param keyEnd the end iterator of the keys
   * @param valResBegin the begin iterator of the value results
   * @param hashIndices the indices of the keys in the hash table
   */
  template <class KeyIter, class ValResIter>
  void findTableBatchDeferWriteBack(int tableNum, const KeyIter keyBegin,
                                    const KeyIter keyEnd,
                                    ValResIter valResBegin,
                                    std::vector<PositionType>& hashIndices) {
    static_assert(isOblivious);

    Assert(tableNum == 0 || tableNum == 1);

    size_t keySize = std::distance(keyBegin, keyEnd);
    Assert(hashIndices.size() == keySize);
    std::vector<PositionType> recoveryVec(keySize);
    for (PositionType i = 0; i < keySize; ++i) {
      recoveryVec[i] = i;
    }
    EM::Algorithm::BitonicSortSepPayload(hashIndices.begin(), hashIndices.end(),
                                         recoveryVec.begin());
    std::vector<BucketType> buckets(keySize);

    if (tableNum == 0) {
      table0.BatchReadDeferWriteBack(hashIndices, buckets);
    } else {
      table1.BatchReadDeferWriteBack(hashIndices, buckets);
    }

    EM::Algorithm::BitonicSortSepPayload(recoveryVec.begin(), recoveryVec.end(),
                                         buckets.begin());
    for (size_t i = 0; i < keySize; ++i) {
      (valResBegin + i)->found =
          searchBucket(*(keyBegin + i), (valResBegin + i)->value, buckets[i]);
    }
  }

  /**
   * @brief Find the values of a batch of keys obliviously in the two tables in
   * parallel. Defer the recursive ORAM writeback procedure. The input may have
   * duplicate keys and may be arranged in any order.
   *
   * @tparam KeyIter
   * @tparam ValResIter
   * @param keyBegin
   * @param keyEnd
   * @param valResBegin
   */
  template <class KeyIter, class ValResIter>
  void findBatchDeferWriteBack(const KeyIter keyBegin, const KeyIter keyEnd,
                               ValResIter valResBegin) {
    size_t keySize = std::distance(keyBegin, keyEnd);
    std::vector<std::vector<ValResult>> valResTables(
        2, std::vector<ValResult>(keySize));
    std::vector<std::vector<PositionType>> hashIndices(
        2, std::vector<PositionType>(keySize));

    for (size_t i = 0; i < keySize; ++i) {
      indexer.getHashIndices(*(keyBegin + i), hashIndices[0][i],
                             hashIndices[1][i]);
    }

#pragma omp parallel for num_threads(2) schedule(static, 1)
    for (int i = 0; i < 2; ++i) {
      findTableBatchDeferWriteBack(i, keyBegin, keyEnd, valResTables[i].begin(),
                                   hashIndices[i]);
    }

    for (size_t i = 0; i < keySize; ++i) {
      mergeValRes(valResTables[0][i], valResTables[1][i], *(valResBegin + i));
      (valResBegin + i)->found |=
          searchStash(*(keyBegin + i), (valResBegin + i)->value, stash);
    }
  }

  /**
   * @brief Perform writeback of the recursive ORAM.
   *
   * @param tableNum the table number to write back, 0 or 1
   */
  void writeBackTable(int tableNum) {
    static_assert(isOblivious);
    Assert(tableNum == 0 || tableNum == 1);
    if (tableNum == 0) {
      table0.WriteBack();
    } else {
      table1.WriteBack();
    }
  }

  /**
   * @brief Erase a key from the hash map. The function is not oblivious.
   *
   * @param key the key to erase
   * @param isDummy whether the erase is dummy operation
   * @return true if the key is found, false otherwise
   */
  bool erase(const K& key, bool isDummy = false) {
    if (isDummy) {
      return false;
    }
    bool erased = false;
    auto updateFunc = [&](BucketType& bucket) {
      for (int i = 0; i < bucketSize; ++i) {
        auto& entry = bucket.entries[i];
        if (entry.valid && entry.key == key) {
          entry.valid = false;
          erased = true;
          return;
        }
      }
    };
    PositionType idx0 = indexer.getHashIdx0(key);
    updateHelper(idx0, table0, updateFunc);
    if (erased) {
      --load;
      return true;
    }
    PositionType idx1 = indexer.getHashIdx1(key);
    updateHelper(idx1, table1, updateFunc);
    if (erased) {
      --load;
      return true;
    }
    for (auto it = stash.begin(); it != stash.end(); ++it) {
      if (it->key == key && it->valid) {
        --load;
        stash.erase(it);
        return true;
      }
    }
    return false;
  }

  /**
   * @brief Erase a key from the hash map. The function is oblivious.
   *
   * @param key the key to erase
   * @param isDummy whether the erase is dummy operation
   * @return true if the key is found, false otherwise
   */
  bool eraseOblivious(const K& key, bool isDummy = false) {
    if (isDummy) {
      return false;
    }
    bool erased = false;
    auto updateFunc = [&](BucketType& bucket) {
      for (int i = 0; i < bucketSize; ++i) {
        auto& entry = bucket.entries[i];
        bool matchFlag = entry.valid & (entry.key == key);
        entry.valid &= !matchFlag | erased;
        erased |= matchFlag;
      }
    };
    PositionType idx0, idx1;
    indexer.getHashIndices(key, idx0, idx1);
    updateHelper(idx0, table0, updateFunc);
    updateHelper(idx1, table1, updateFunc);
    for (auto& entry : stash) {
      bool matchFlag = entry.valid & (entry.key == key);
      entry.valid &= !matchFlag | erased;
      erased |= matchFlag;
    }
    load -= erased;
    return erased;
  }

 private:
  /**
   * @brief Helper function that merges the value results of two tables.
   *
   * @param valRes0 the value result of table 0
   * @param valRes1 the value result of table 1
   * @param valRes the merged value result
   */
  static void mergeValRes(const ValResult& valRes0, const ValResult& valRes1,
                          ValResult& valRes) {
    valRes.found = valRes0.found | valRes1.found;
    valRes.value = valRes0.value;
    obliMove(valRes1.found, valRes.value, valRes1.value);
  }

};  // namespace ODSL
}  // namespace ODSL