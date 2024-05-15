#pragma once
#include "page_oram.hpp"
#include "recursive_oram.hpp"

/// @brief This file implements a cuckoo hash map built on top of the
/// recursive ORAM. Specifically, it contains two hash tables implemented with
/// recursive oram. Each key is hashed to two positions, one for each table, and
/// at each position of a table is a bucket with two slots, the element can be
/// stored in any of these four slots. In addition, there is a stash for
/// elements that cannot be stored in the tables.

namespace ODSL {

/**
 * @brief Used to hash the key to two positions
 *
 */
template <typename K, typename PositionType = uint64_t>
struct OHashMapIndexer {
  static constexpr int saltLength = 16;  // 128 bits salt
  uint8_t salts[saltLength];             // salt for secure hash
  PositionType _size;                    // number of buckets in each hash table

  OHashMapIndexer() {}

  explicit OHashMapIndexer(PositionType size) { SetSize(size); }

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
    pos0 = (PositionType)(hashIndices.h0 % _size);
    pos1 = (PositionType)(hashIndices.h1 % _size);
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
 * @brief A stash storing elements that cannot be stored in the hash tables.
 * Helps deamortize the cost of oblivious insertion.
 *
 * @tparam K the key type
 * @tparam V the value type
 * @tparam stash_size the default stash size for oblivious insertion
 */
template <typename K, typename V, const uint64_t stash_size = 16>
struct LRUStash {
  using KVEntry = OHashMapEntry<K, V>;
  std::vector<KVEntry> stash;        // the stash data
  std::vector<uint64_t> timestamps;  // the timestamp each entry is inserted
  uint64_t currTime;                 // the current timestamp
  // as long as one oblivious insert occur, we cannot reveal the state of the
  // stash
  bool oInserted = false;  // whether an oblivious insert has occurred

  LRUStash() : currTime(0) {}

  explicit LRUStash(size_t capacity) : currTime(0) { SetSize(capacity); }

  /**
   * @brief Set the stash size for oblivious insertion
   *
   * @param capacity the stash size
   */
  void SetSize(size_t capacity) {
    stash.resize(capacity);
    timestamps.resize(capacity);
  }

  /**
   * @brief Obliviously insert an entry to the stash and record the timestamp.
   * If entry.valid is false, the insertion is dummy. If the stash overflows,
   * the method will enlarge the stash, which is not oblivious.
   *
   * @tparam highPriority whether the entry should be populated first
   * @param entry the entry to insert
   */
  template <const bool highPriority = false>
  void OInsert(const KVEntry& entry) {
    oInserted = true;
    if (stash.size() < stash_size) {
      stash.resize(stash_size);
      timestamps.resize(stash_size);
    }
    bool inserted = !entry.valid;
    uint64_t time = highPriority ? 0 : currTime;
    for (size_t i = 0; i < stash.size(); ++i) {
      bool isEmpty = !stash[i].valid;
      bool insertFlag = isEmpty & (!inserted);
      obliMove(insertFlag, stash[i], entry);
      obliMove(insertFlag, timestamps[i], time);
      inserted |= insertFlag;
    }
    bool overflowFlag = (!inserted) & entry.valid;
    if (overflowFlag) {
      PERFCTR_INCREMENT(OHMAP_DEAMORT_OVERFLOW);
      stash.push_back(entry);
      timestamps.push_back(time);
    }
    // reset the timestamps if the current time overflows
    if (currTime == UINT64_MAX) {
      std::fill(timestamps.begin(), timestamps.end(), 0);
    }
    ++currTime;
  }

  /**
   * @brief Insert an entry to the stash and record the timestamp.
   * If previously an oblivious insert has occurred, the method will insert
   * obliviously. Otherwise, it will directly push the entry to the stash.
   *
   * @param entry the entry to insert
   */
  void Insert(const KVEntry& entry) {
    stash.push_back(entry);
    timestamps.push_back(currTime);
  }

  /**
   * @brief Read the oldest entry in the stash and remove it. If the stash is
   * empty, entry.valid will be set to false.
   *
   * @param entry the oldest entry
   */
  void OPopOldest(KVEntry& entry) {
    uint64_t oldestTime = currTime;
    size_t oldestIdx = stash.size();
    for (size_t i = 0; i < stash.size(); ++i) {
      bool isOldest = (timestamps[i] <= oldestTime) & stash[i].valid;
      obliMove(isOldest, oldestTime, timestamps[i]);
      obliMove(isOldest, entry, stash[i]);
      obliMove(isOldest, oldestIdx, i);
    }
    for (size_t i = 0; i < stash.size(); ++i) {
      stash[i].valid &= (i != oldestIdx);
    }
    entry.valid &= (oldestIdx != stash.size());
  }

  using Iter = typename std::vector<KVEntry>::iterator;

  /**
   * @brief Non-obliviously erase an entry from the stash
   *
   * @param it the iterator to the entry to erase
   */
  void Erase(Iter it) {
    if (oInserted) {
      throw std::runtime_error("Cannot erase after oblivious insertion.");
    }
    uint64_t idx = std::distance(stash.begin(), it);
    stash.erase(it);
    timestamps.erase(timestamps.begin() + idx);
  }

  /**
   * @brief Returns the beginning iterator of the internal stash data
   *
   */
  Iter begin() { return stash.begin(); }

  /**
   * @brief Returns the end iterator of the internal stash data
   *
   */
  Iter end() { return stash.end(); }

  /**
   * @brief Returns the size of the stash
   *
   */
  size_t size() const { return stash.size(); }

  /**
   * @brief Bracket operator to access the stash
   *
   * @param idx the index to access
   * @return KVEntry& the entry at the index
   */
  KVEntry& operator[](size_t idx) { return stash[idx]; }

  /**
   * @brief Const bracket operator to access the stash
   *
   * @param idx the index to access
   * @return const KVEntry& the entry at the index
   */
  const KVEntry& operator[](size_t idx) const { return stash[idx]; }
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

enum ObliviousLevel { NON_OBLIVIOUS, PAGE_OBLIVIOUS, FULL_OBLIVIOUS };

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
template <typename K, typename V,
          const ObliviousLevel isOblivious = FULL_OBLIVIOUS,
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
  PositionType _size = 0;
  // the number of elements in the hash map
  PositionType load;
  // the size of each hash table
  PositionType tableSize;
  using BucketType = OHashMapBucket<K, V, bucketSize>;
  // for oblivious hash map, we use recursive ORAM
  using ObliviousTableType =
      std::conditional_t<isOblivious == FULL_OBLIVIOUS,
                         RecursiveORAM<BucketType, PositionType>,
                         PageORAM<BucketType, PositionType>>;
  // for non-oblivious hash map, we cache the front of the vector, for the
  // remaining data store it encrypted and authenticated in external memory,
  // check freshness when swapped in.
  using NonObliviousTableType = EM::CacheFrontVector::Vector<
      BucketType, sizeof(BucketType),
      EM::CacheFrontVector::EncryptType::ENCRYPT_AND_AUTH_FRESH, 1024>;
  using TableType = std::conditional_t<isOblivious, ObliviousTableType,
                                       NonObliviousTableType>;
  TableType table0, table1;
  OHashMapIndexer<K, PositionType> indexer;
  using KVEntry = OHashMapEntry<K, V>;
  using StashType = LRUStash<K, V>;
  StashType stash;
  bool inited = false;

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
  template <const bool hideDummy = false, class EntryIterator>
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
  static bool replaceIfExist(StashType& stash, const KVEntry& entryToInsert) {
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
   * @brief Helper function that searches the stash for a key. Depending on
   * whether the hash map is oblivious, the function may reveal the number of
   * comparisons.
   *
   * @param key the key to search
   * @param value the value to return
   * @param stash the stash to search
   * @return true if the key is found, false otherwise
   */
  static bool searchStash(const K& key, V& value, const StashType& stash) {
    if constexpr (isOblivious) {
      bool found = false;
      for (size_t i = 0; i < stash.size(); ++i) {
        const auto& entry = stash[i];
        bool match = entry.valid & (entry.key == key);
        obliMove(match, value, entry.value);
        found |= match;
      }
      return found;
    } else {
      for (size_t i = 0; i < stash.size(); ++i) {
        const auto& entry = stash[i];
        if (entry.valid && entry.key == key) {
          value = entry.value;
          return true;
        }
      }
      return false;
    }
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
  template <class EntryIterator>
  static bool replaceIfExistOblivious(EntryIterator begin, EntryIterator end,
                                      const KVEntry& entryToInsert) {
    bool updated = !entryToInsert.valid;
    for (auto it = begin; it != end; ++it) {
      auto& entry = *it;
      bool matchFlag =
          entry.valid & (entry.key == entryToInsert.key) & (!updated);
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
   * @param stash the stash to perform replace
   * @param entryToInsert the entry to insert
   * @return true if the entry is replaced, false otherwise
   */
  static bool replaceIfExistOblivious(StashType& stash,
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
      bool insertFlag = isEmpty & (!updated);
      obliMove(insertFlag, entry, entryToInsert);
      updated |= insertFlag;
    }
    return updated;
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
  static bool searchBucket(const K& key, V& value, const BucketType& bucket) {
    if constexpr (isOblivious) {
      bool found = false;
      for (int i = 0; i < bucketSize; ++i) {
        const auto& entry = bucket.entries[i];
        bool match = entry.valid & (entry.key == key);
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
    static_assert(isOblivious > NON_OBLIVIOUS);

    Assert(tableNum == 0 || tableNum == 1);

    size_t keySize = std::distance(keyBegin, keyEnd);
    Assert(hashIndices.size() == keySize);
    std::vector<PositionType> recoveryVec(keySize);
    for (PositionType i = 0; i < keySize; ++i) {
      recoveryVec[i] = i;
    }
    Algorithm::BitonicSortSepPayload(hashIndices.begin(), hashIndices.end(),
                                     recoveryVec.begin());
    std::vector<BucketType> buckets(keySize);

    if (tableNum == 0) {
      table0.BatchReadDeferWriteBack(hashIndices, buckets);
    } else {
      table1.BatchReadDeferWriteBack(hashIndices, buckets);
    }

    Algorithm::BitonicSortSepPayload(recoveryVec.begin(), recoveryVec.end(),
                                     buckets.begin());
    for (size_t i = 0; i < keySize; ++i) {
      (valResBegin + i)->found =
          searchBucket(*(keyBegin + i), (valResBegin + i)->value, buckets[i]);
    }
  }

  /**
   * @brief Try to insert entry into either table0 or table1 obliviously without
   * swapping existing elements. If the element already exists either in table
   * 0, table 1, or the stash, replace the existing element. If there's no
   * available slot, swap entryToInsert with a random element from the bucket in
   * table 0. If the insertion is successful, entryToInsert.valid will be set to
   * false, and dummy operations will be performed to ensure oblivousness.
   *
   * @param entryToInsert the entry to insert, and will be modified to the entry
   * swapped out if no slot is available.
   *
   * @return true if the key already exists, false otherwise
   */
  bool insertEntryOblivious(KVEntry& entryToInsert) {
    PositionType idx0, idx1;
    indexer.getHashIndices(entryToInsert.key, idx0, idx1);
    obliMove(!entryToInsert.valid, idx0, UniformRandom(tableSize - 1));
    obliMove(!entryToInsert.valid, idx1, UniformRandom(tableSize - 1));
    bool exist = false;
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
        int offset = (int)(UniformRandom32() % bucketSize);
        obliSwap(entryToInsert.valid, bucket0.entries[offset], entryToInsert);
      };
      updateHelper(idx1, table1, table1UpdateFunc);
    };
    updateHelper(idx0, table0, table0UpdateFunc);

    load += !exist;
    return exist;
  }

  /**
   * @brief Try to insert entry into table 1 obliviously, if table 1 is
   * occupied, swap a random element out of the bucket in table 1 and try to
   * insert this element into table 0. If table 0 is also occupied, swap a
   * random element out of the bucket in table 0 and save it in entryToInsert.
   * If either of the insertion succeeds, entryToInsert.valid will become false,
   * and dummy operations will be performed to ensure oblivoiusness. The method
   * may retry multiple times.
   *
   * @param entryToInsert the entry to insert, and will be modified to the entry
   * swapped out if no slot is available.
   * @param maxRetry the maximum number of retries, default to 1
   */
  void insertEntryObliviousRetry(KVEntry& entryToInsert, int maxRetry = 1) {
    auto swapUpdateFunc = [&](BucketType& bucket) {
      int offset = (int)(UniformRandom32() % bucketSize);
      bool insertSucceed = insertIfEmptyOblivious(bucket, entryToInsert);
      entryToInsert.valid &= !insertSucceed;
      obliSwap(entryToInsert.valid, entryToInsert, bucket.entries[offset]);
    };
    for (int r = 0; r < maxRetry; ++r) {
      PositionType idx1 = indexer.getHashIdx1(entryToInsert.key);
      updateHelper(idx1, table1, swapUpdateFunc);  // modifies entryToInsert

      PositionType idx0 = indexer.getHashIdx0(entryToInsert.key);
      updateHelper(idx0, table0, swapUpdateFunc);
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
  explicit OHashMap(PositionType size) { SetSize(size); }

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
  void SetSize(PositionType size) { SetSize(size, MAX_CACHE_SIZE); }

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
    tableSize =
        (PositionType)((double)size / (2 * loadFactor * bucketSize) + 1);
    indexer.SetSize(tableSize);
    if constexpr (isOblivious) {
      // equally divide the cache between the two tables
      table0.SetSize(tableSize, cacheBytes / 2);
      table1.SetSize(tableSize, cacheBytes / 2);
      stash.SetSize(16);
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

  const StashType& GetStash() const { return stash; }

  using NonObliviousHashMap = OHashMap<K, V, NON_OBLIVIOUS, PositionType>;

  /**
   * @brief Initialize the oblivious hash map from another non-oblivious hash
   * map
   *
   * @param other the non-oblivious hash map
   */
  void InitFromNonOblivious(NonObliviousHashMap& other) {
    static_assert(isOblivious > NON_OBLIVIOUS,
                  "Only oblivious hash map can call this function");
    if (_size != other.size()) {
      throw std::runtime_error("OHashMap InitFromNonOblivious failed");
    }
    if (_size == 0) {
      throw std::runtime_error(
          "OHashMap size not set. Call SetSize before initialization.");
    }
    PositionType load0 = 0;
    PositionType load1 = 0;
    indexer = other.GetIndexer();
    // change dummies to invalid
    EM::VirtualVector::VirtualReader<BucketType> reader0(
        other.GetTableSize(), [&](PositionType i) {
          BucketType bucket = other.GetTable0()[i];
          for (int j = 0; j < bucketSize; ++j) {
            bucket.entries[j].valid &= !bucket.entries[j].dummy;
            load0 += bucket.entries[j].valid;
            bucket.entries[j].dummy = false;
          }
          return bucket;
        });
    EM::VirtualVector::VirtualReader<BucketType> reader1(
        other.GetTableSize(), [&](PositionType i) {
          BucketType bucket = other.GetTable1()[i];
          for (int j = 0; j < bucketSize; ++j) {
            bucket.entries[j].valid &= !bucket.entries[j].dummy;
            load1 += bucket.entries[j].valid;
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
    load = load0 + load1;
    for (const auto& entry : other.GetStash().stash) {
      if (entry.valid) {
        // valid is public but dummy is not
        // since we deleted all dummies, it's likely that we don't need to put
        // these entries in the stash
        OInsert(entry.key, entry.value, entry.dummy);
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
    requires Readable<Reader, std::pair<K, V>>
  void InitFromReader(Reader& reader, uint64_t additionalCacheBytes = 0) {
    if (inited) {
      throw std::runtime_error("OHashMap initialized twice.");
    }
    inited = true;
    if (_size == 0) {
      throw std::runtime_error(
          "OHashMap size not set. Call SetSize before initialization.");
    }
    if constexpr (!isOblivious) {
      while (!reader.eof()) {
        const std::pair<K, V>& entry = reader.read();
        Insert(entry.first, entry.second);
      }
    } else {
      additionalCacheBytes = std::max(
          additionalCacheBytes, NonObliviousTableType::GetMinMemoryUsage() * 2);
      NonObliviousHashMap nonObliviousHashMap(_size, additionalCacheBytes);
      nonObliviousHashMap.InitFromReader(reader);
      InitFromNonOblivious(nonObliviousHashMap);
    }
  }

  /**
   * @brief An object that stores the curret state of initialization. Faciliates
   * initialization in a streaming fashion.
   *
   */
  struct InitContext {
   private:
    NonObliviousHashMap* nonObliviousHashMap;
    OHashMap& oHashMap;  // the parent map

   public:
    /**
     * @brief Construct a new InitContext object
     *
     * @param omap The parent map
     * @param additionalCacheBytes the size of additional available cache for
     * initialization
     */
    explicit InitContext(OHashMap& map, uint64_t additionalCacheBytes = 0)
        : oHashMap(map),
          nonObliviousHashMap(
              new NonObliviousHashMap(map.size(), additionalCacheBytes)) {
      if (map.inited) {
        throw std::runtime_error("OHashMap initialized twice.");
      }
      map.inited = true;
      if (map.size() == 0) {
        throw std::runtime_error(
            "OHashMap size not set. Call SetSize before initialization.");
      }
    }

    /**
     * @brief Move constructor
     * @param other the other InitContext
     *
     */
    InitContext(const InitContext&& other)
        : oHashMap(other.oHashMap),
          nonObliviousHashMap(other.nonObliviousHashMap) {}

    InitContext(const InitContext& other) = delete;

    /**
     * @brief Insert a new key value pair for initialization. The method will
     * reveal whether the key has been inserted before. But if all the keys are
     * distinct, the operation is oblivious. The method will throw an exception
     * if too many keys are inserted.
     *
     * @param key
     * @param value
     */
    void Insert(const K& key, const V& value) {
      nonObliviousHashMap->Insert(key, value);
    }

    void Insert(const std::pair<K, V>& entry) {
      nonObliviousHashMap->Insert(entry.first, entry.second);
    }

    template <class Iterator>
    void InsertBatch(Iterator begin, Iterator end) {
      for (auto it = begin; it != end; ++it) {
        nonObliviousHashMap->Insert(it->first, it->second);
      }
    }

    /**
     * @brief Finalize the initialization. The method will copy the data from
     * the non-oblivious hash map to the oblivious hash map.
     *
     */
    void Finalize() {
      oHashMap.InitFromNonOblivious(*nonObliviousHashMap);
      delete nonObliviousHashMap;
    }
  };

  /**
   * @brief Obtain a new context to initialize this map. The initialization data
   * can be either private or public (the initialization and subsequent accesses
   * are oblivious).
   * Example:
   *  auto* initContext = oMap.NewInitContext(1UL << 28);
      for (auto it = kvMap.begin(); it != kvMap.end(); ++it;) {
        initContext->Insert(it->first, it->second);
      }
      initContext->Finalize();
      delete initContext;
   *
   * @param additionalCacheBytes
   * @return InitContext
   */
  InitContext* NewInitContext(uint64_t additionalCacheBytes = 0) {
    static_assert(isOblivious > NON_OBLIVIOUS,
                  "Only oblivious hash map can call this function");
    return new InitContext(*this, additionalCacheBytes);
  }

  /**
   * @brief Initialize an empty hash map
   *
   */
  void Init() {
    if (inited) {
      throw std::runtime_error("OHashMap initialized twice.");
    }
    inited = true;
    if (_size == 0) {
      throw std::runtime_error(
          "OHashMap size not set. Call SetSize before initialization.");
    }
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
   * @brief Non-oblivious insert, but may hide whether the insertion is dummy.
   *
   * @tparam hideDummy whether to hide the dummy flag of the inserted entry
   * @param key the key to insert
   * @param value the value to insert
   * @param isDummy whether the insertion is dummy
   * @return true if the key already exists, false otherwise
   */
  template <bool hideDummy = false>
  bool Insert(const K& key, const V& value, bool isDummy = false) {
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
    PositionType idx1;
    bool inserted = false;
    bool exist = false;

    auto table0UpdateFunc = [&](BucketType& bucket0) {
      if (replaceIfExist<hideDummy>(bucket0, entryToInsert)) {
        inserted = true;
        exist = true;
        // found in table 0
        return;
      }
      idx1 = indexer.getHashIdx1(entryToInsert.key);
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
          return;
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
      offset = UniformRandom32() % bucketSize;
      idx1 = indexer.getHashIdx1(entryToInsert.key);
      updateHelper(idx1, table1, swapUpdateFunc);
      if (inserted) {
        return false;
      }
      idx0 = indexer.getHashIdx0(entryToInsert.key);
      offset = UniformRandom32() % bucketSize;
      updateHelper(idx0, table0, swapUpdateFunc);
      if (inserted) {
        return false;
      }
    }
    stash.Insert(entryToInsert);
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
  bool OInsert(const K& key, const V& value, bool isDummy = false) {
    KVEntry entryToInsert = {!isDummy, false, key, value};
    bool exist = insertEntryOblivious(entryToInsert);
    // the element just swapped out is more likely to get inserted to somewhere
    // else
    stash.template OInsert<true>(entryToInsert);

    for (int i = 0; i < 2; ++i) {
      // use FIFO order so that we won't get stuck by loops in the random graph
      // of cuckoo hashing
      stash.OPopOldest(entryToInsert);
      insertEntryObliviousRetry(entryToInsert, 1);
      stash.OInsert(entryToInsert);
    }

    return exist;
  }

  /**
   * @brief Update the value of a key. No operation is performed if the key does
   * not exist. Depending on whether the hash map is oblivious, the function may
   * reveal the number of comparisons.
   *
   * @param key the key to update
   * @param value the value to update
   * @param isDummy whether the search is dummy
   * @return true if the key is found, false otherwise
   */
  bool Update(const K& key, const V& value, bool isDummy = false) {
    bool found = false;
    if constexpr (!isOblivious) {
      if (isDummy) {
        return false;
      }
      PositionType idx0 = indexer.getHashIdx0(key);
      auto updateFunc = [&](BucketType& bucket) {
        for (int i = 0; i < bucketSize; ++i) {
          auto& entry = bucket.entries[i];
          if (entry.valid && entry.key == key) {
            entry.value = value;
            found = true;
            return;
          }
        }
      };
      updateHelper(idx0, table0, updateFunc);
      if (found) {
        return true;
      }
      PositionType idx1 = indexer.getHashIdx1(key);
      updateHelper(idx1, table1, updateFunc);
      if (found) {
        return true;
      }
      // found = searchStash(key, value, stash);
      for (auto it = stash.begin(); it != stash.end(); ++it) {
        if (it->key == key && it->valid) {
          it->value = value;
          return true;
        }
      }
      return false;
    } else {
      PositionType idx0, idx1;
      indexer.getHashIndices(key, idx0, idx1);

      obliMove(isDummy, idx0, UniformRandom(tableSize - 1));

      obliMove(isDummy, idx1, UniformRandom(tableSize - 1));

      auto updateFunc = [&](BucketType& bucket) {
        for (int i = 0; i < bucketSize; ++i) {
          auto& entry = bucket.entries[i];
          bool matchFlag = entry.valid & (entry.key == key);
          obliMove(matchFlag, entry.value, value);
          found |= matchFlag;
        }
      };
      updateHelper(idx0, table0, updateFunc);
      updateHelper(idx1, table1, updateFunc);
      for (auto it = stash.begin(); it != stash.end(); ++it) {
        bool matchFlag = (it->key == key) & it->valid;
        obliMove(matchFlag, it->value, value);
        found |= matchFlag;
      }
      return found & (!isDummy);
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
  bool Find(const K& key, V& value, bool isDummy = false) {
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
      return found & (!isDummy);
    } else {
      bool found = false;
      PositionType idx0, idx1;
      indexer.getHashIndices(key, idx0, idx1);

      obliMove(isDummy, idx0, UniformRandom(tableSize - 1));

      BucketType bucket;
      readHelper(idx0, table0, bucket);
      found = searchBucket(key, value, bucket);

      obliMove(isDummy, idx1, UniformRandom(tableSize - 1));
      readHelper(idx1, table1, bucket);
      found |= searchBucket(key, value, bucket);
      found |= searchStash(key, value, stash);
      return found & (!isDummy);
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
  void FindBatchDeferWriteBack(const KeyIter keyBegin, const KeyIter keyEnd,
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
  void WriteBackTable(int tableNum) {
    static_assert(isOblivious > NON_OBLIVIOUS,
                  "Only oblivious hash map can call this function");
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
  bool Erase(const K& key, bool isDummy = false) {
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
        stash.Erase(it);
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
  bool OErase(const K& key, bool isDummy = false) {
    if (isDummy) {
      return false;
    }
    bool erased = false;
    auto updateFunc = [&](BucketType& bucket) {
      for (int i = 0; i < bucketSize; ++i) {
        auto& entry = bucket.entries[i];
        bool matchFlag = entry.valid & (entry.key == key);
        entry.valid &= (!matchFlag) | erased;
        erased |= matchFlag;
      }
    };
    PositionType idx0, idx1;
    indexer.getHashIndices(key, idx0, idx1);
    updateHelper(idx0, table0, updateFunc);
    updateHelper(idx1, table1, updateFunc);
    for (auto& entry : stash) {
      bool matchFlag = entry.valid & (entry.key == key);
      entry.valid &= (!matchFlag) | erased;
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