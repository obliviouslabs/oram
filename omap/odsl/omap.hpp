#pragma once
#include "recursive_oram.hpp"

/// @brief This file implements an oblivious unordered map. The position map is
/// implemented with a cuckoo hash table on top of a recursive ORAM. The values
/// are stored separatedly in a circuit oram. Separating the values to a
/// separate oram reduces storage overhead and speeds up the query, insertion,
/// and erasure especially when keys and values are moderately large (> 32
/// bytes).
namespace ODSL {
/**
 * @brief Hash the key to two positions and use remaining hash bits to
 * distinguish entries.
 * @tparam K the key type
 * @tparam H the type of the hash value stored in the position map
 * @tparam HExtra the type of the extra hash value stored in the main oram
 * @tparam PositionType the type of the position,
 *
 */
template <typename K, typename H, typename HExtra, typename PositionType>
struct OPosMapIndexer {
  static constexpr int saltLength = 16;  // 128 bits salt
  uint8_t salts[saltLength];             // salt for secure hash
  PositionType _size;                    // number of buckets in each hash table
  static constexpr uint64_t uid_length =
      std::min(sizeof(PositionType) * 2 + sizeof(H), 16UL);
  using UidType = Bytes<uid_length>;  // type of the unique identifier
  OPosMapIndexer() {}

  explicit OPosMapIndexer(PositionType size) { SetSize(size); }

  /**
   * @brief Initialize the indexer
   *
   * @param size number of buckets in each hash table
   */
  void SetSize(PositionType size) {
    _size = size;
    read_rand(salts, saltLength);
  }

  /**
   * @brief Data structure that receives the output of SHA256
   *
   */
  struct HashIndices {
    uint64_t h0;
    uint64_t h1;
    H keyHash;
    HExtra extraHash;
  };

  /**
   * @brief Get both hash indices
   *
   * @param key the key to hash
   * @param pos0 the first hash index
   * @param pos1 the second hash index
   */
  void getHashIndices(const K& key, PositionType& pos0, PositionType& pos1,
                      H& keyHash, UidType& uid, HExtra& extraHash) const {
    HashIndices hashIndices;
    secure_hash_with_salt(key, salts, (uint8_t*)&hashIndices,
                          sizeof(hashIndices));
    // position in table 0
    pos0 = (PositionType)(hashIndices.h0 % _size);
    // position in table 1
    pos1 = (PositionType)(hashIndices.h1 % _size);
    // keyHash is used to distinguish entries that have the same positions in
    // both tables. Set the least significant bit to 1 to indicate a valid
    // entry. In the position map, uid = (pos0 || pos1 || keyHash) must be
    // unique. If a newly inserted key-value pair has the same uid as an
    // existing entry, we will insert it to the stash.
    keyHash = hashIndices.keyHash | (H)1;
    // extraHash is used to further distinguish entries that shares the same
    // uid. It is stored in the main oram, so that we can check the
    // existence of an element with very high confidence.
    extraHash = hashIndices.extraHash;
    // generate a unique identifier that includes the key hash and the
    // positions, and with additional bits to distinguish elements
    std::memset(uid.data, 0, sizeof(UidType));
    int sizeMaxLen = GetLogBaseTwo(_size - 1) / 8 + 1;
    Assert(2 * sizeMaxLen + sizeof(H) <= uid_length);
    for (int i = 0; i < sizeMaxLen; ++i) {
      uid.data[i] = (uint8_t)(pos0 >> (i * 8));
      uid.data[i + sizeMaxLen] = (uint8_t)(pos1 >> (i * 8));
    }
    std::memcpy(uid.data + 2 * sizeMaxLen, &keyHash, sizeof(H));
    // std::cout << "pos0: " << pos0 << ", pos1: " << pos1
    //           << ", keyHash: " << keyHash << ", extraHash: " << extraHash
    //           << ", uid: " << uid << std::endl;
  }
};

/**
 * @brief A single entry (slot) in the position map. To compress the size of the
 * position map, we do not store the key, but store 1) the index of the
 * key in the other table, which provides some entropy to distinguish entries in
 * the same slot, and moreover, allows us to swap the entry with the other table
 * when the entry is evicted 2) an extra hash value "keyHash" to further
 * distinguish the entries that maps to same same slots in both tables. The
 * least significant bit of keyHash is used to indicate whether the entry is
 * valid. In case both the indices and the keyHash collides, we store the entry
 * in a stash.
 *
 * @tparam K the key type
 * @tparam V the value type, which is the position of the value in the oram
 */
template <typename K, typename V, typename H, typename PositionType>
struct GenericOPosMapEntry {
  H keyHash = 0;
  PositionType otherIdx;  // the index in the other table
  V value;  // the value of the entry, i.e., the position of the entry in the
            // oram that stores the value

  /**
   * @brief Returns whether the entry is valid
   */
  bool valid() const { return keyHash & 0x1; }
  /**
   * @brief Set the entry to valid if real is true
   */
  void setValid(bool real) { keyHash |= real; }
  /**
   * @brief Set the entry to invalid if real is true
   */
  void setInvalid(bool real) { keyHash &= ~((H)real); }
#ifndef ENCLAVE_MODE
  // cout
  friend std::ostream& operator<<(std::ostream& os,
                                  const GenericOPosMapEntry& entry) {
    os << "(" << entry.valid() << ", " << entry.keyHash << ", "
       << entry.otherIdx << ", " << entry.value << ")";
    return os;
  }
#endif
};

/**
 * @brief A bucket can contain multiple entries that are fully associative.
 * This can significantly increase the load factor of the cuckoo hash table.
 *
 * @tparam K the key type
 * @tparam V the value type
 * @tparam bucketSize the number of slots in the bucket
 */
template <const short bucketSize, typename K, typename V, typename H,
          typename PositionType>
struct OPosMapBucket {
  using OPosMapEntry = GenericOPosMapEntry<K, V, H, PositionType>;
  OPosMapEntry entries[bucketSize];
#ifndef ENCLAVE_MODE
  // cout
  friend std::ostream& operator<<(std::ostream& os,
                                  const OPosMapBucket& bucket) {
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
template <typename K, typename PositionType = uint64_t,
          const bool isOblivious = true, const bool parallel_init = true>
struct OPosMap {
  using V = PositionType;   // value type is the position
  using H = uint32_t;       // additional hash bits to distinguish elements
  using HExtra = uint64_t;  // extra hash bits stored in the main oram to check
                            // for the existence of an element
  using Indexer = OPosMapIndexer<K, H, HExtra, PositionType>;

  // assume that the oram size < 2^48
  using UidType = typename Indexer::UidType;  // type of the unique identifier
  using OPosMapEntry = GenericOPosMapEntry<K, V, H, PositionType>;
  static constexpr bool isObliviousPosMap = isOblivious;
  /**
   * @brief A stash storing position map entries that cannot be stored in the
   * hash tables because both buckets are full. The stash works as a queue and
   * helps deamortize the cost of oblivious insertion. Notice that this stash is
   * different from the stash in the main oram. The latter is used to store
   * indistinguishable keys.
   *
   * @tparam K the key type
   * @tparam V the value type
   * @tparam stash_size the default stash size for oblivious insertion
   */
  template <const uint64_t stash_size = 16>
  struct LRUStash {
    struct StashEntry {
      // the entry for table 0
      OPosMapEntry entry;
      // the index of the entry in table 0
      PositionType idx0;
      // the timestamp when the entry is inserted
      uint64_t timestamp;
    };
    using StashVec = std::vector<StashEntry>;
    StashVec stash;     // the stash data
    uint64_t currTime;  // the current timestamp
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
    void SetSize(size_t capacity) { stash.resize(capacity); }

    /**
     * @brief Obliviously insert an entry to the stash and record the timestamp.
     * If entry.valid is false, the insertion is dummy. If the stash overflows,
     * the method will enlarge the stash, which is not oblivious.
     *
     * @tparam highPriority whether the entry should be populated first
     * @param entry the entry to insert
     */
    template <const bool highPriority = false>
    void OInsert(const OPosMapEntry& entry, PositionType idx0) {
      oInserted = true;
      if (stash.size() < stash_size) {
        stash.resize(stash_size);
      }
      bool inserted = !entry.valid();
      uint64_t time = highPriority ? 0 : currTime;
      for (size_t i = 0; i < stash.size(); ++i) {
        bool isEmpty = !stash[i].entry.valid();
        bool insertFlag = isEmpty & (!inserted);
        obliMove(insertFlag, stash[i], {entry, idx0, time});

        inserted |= insertFlag;
      }
      bool overflowFlag = (!inserted) & entry.valid();
      if (overflowFlag) {
        PERFCTR_INCREMENT(OHMAP_DEAMORT_OVERFLOW);
        stash.push_back({entry, idx0, time});
      }
      // reset the timestamps if the current time overflows
      if (currTime == UINT64_MAX) {
        for (size_t i = 0; i < stash.size(); ++i) {
          stash[i].timestamp = 0;
        }
      }
      ++currTime;
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
    void Insert(const OPosMapEntry& entry, PositionType idx0) {
      if (oInserted) {
        throw std::runtime_error(
            "Cannot insert non-obliviously to the stash after an oblivious "
            "insert");
      }
      if (!entry.valid()) {
        return;
      }
      if (stash.size() < stash_size) {
        stash.resize(stash_size);
      }
      bool inserted = false;
      uint64_t time = highPriority ? 0 : currTime;
      for (size_t i = 0; i < stash.size(); ++i) {
        bool isEmpty = !stash[i].entry.valid();
        if (isEmpty) {
          stash[i] = {entry, idx0, time};
          inserted = true;
          break;
        }
      }
      bool overflowFlag = !inserted;
      if (overflowFlag) {
        PERFCTR_INCREMENT(OHMAP_DEAMORT_OVERFLOW);
        stash.push_back({entry, idx0, time});
      }
      // reset the timestamps if the current time overflows
      if (currTime == UINT64_MAX) {
        for (size_t i = 0; i < stash.size(); ++i) {
          stash[i].timestamp = 0;
        }
      }
      ++currTime;
    }

    /**
     * @brief Read the oldest entry in the stash and remove it. If the stash is
     * empty, entry.valid will be set to false.
     *
     * @param entry the oldest entry
     */
    void OPopOldest(OPosMapEntry& entry, PositionType& idx0) {
      StashEntry oldestEntry;
      idx0 = 0;
      oldestEntry.idx0 = 0;
      uint64_t oldestTime = currTime;
      size_t oldestIdx = stash.size();
      for (size_t i = 0; i < stash.size(); ++i) {
        bool isOldest =
            (stash[i].timestamp <= oldestTime) & stash[i].entry.valid();
        obliMove(isOldest, oldestEntry, stash[i]);
        obliMove(isOldest, oldestIdx, i);
      }
      for (size_t i = 0; i < stash.size(); ++i) {
        stash[i].entry.setInvalid(i == oldestIdx);
      }
      entry = oldestEntry.entry;
      idx0 = oldestEntry.idx0;
      entry.setInvalid(oldestIdx == stash.size());
    }

    /**
     * @brief Read the oldest entry in the stash and remove it. If the stash is
     * empty, entry.valid will be set to false. This function is not oblivious.
     *
     * @param entry the oldest entry
     */
    void PopOldest(OPosMapEntry& entry, PositionType& idx0) {
      StashEntry oldestEntry;
      idx0 = 0;
      oldestEntry.idx0 = 0;
      uint64_t oldestTime = currTime;
      size_t oldestIdx = stash.size();
      for (size_t i = 0; i < stash.size(); ++i) {
        bool isOldest =
            (stash[i].timestamp <= oldestTime) && stash[i].entry.valid();
        if (isOldest) {
          oldestIdx = i;
        }
      }
      if (oldestIdx != stash.size()) {
        entry = stash[oldestIdx].entry;
        idx0 = stash[oldestIdx].idx0;
        stash[oldestIdx].entry.setInvalid(true);
      } else {
        entry.setInvalid(true);
      }
    }

    using Iter = StashVec::iterator;

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
     * @return OPosMapEntry& the entry at the index
     */
    StashEntry& operator[](size_t idx) { return stash[idx]; }

    /**
     * @brief Const bracket operator to access the stash
     *
     * @param idx the index to access
     * @return const OPosMapEntry& the entry at the index
     */
    const StashEntry& operator[](size_t idx) const { return stash[idx]; }
  };

 private:
  // the ratio between the capacity and the total number of slots in
  // the hash tables
  static constexpr double loadFactor = 0.8;
  // number of slots in each bucket
  static constexpr short bucketSize = 4;
  // maximum number of elements in the stash
  static constexpr int stash_max_size = 10;
  // the capcity of the hash map
  PositionType _size = 0;
  // the number of elements in the hash map
  PositionType load;
  // the size of each hash table
  PositionType tableSize;
  using BucketType = OPosMapBucket<bucketSize, K, V, H, PositionType>;
  // for oblivious hash map, we use recursive ORAM
  using ObliviousTableType = RecursiveORAM<BucketType, PositionType>;
  // for non-oblivious hash map, we cache the front of the vector, for the
  // remaining data store it encrypted and authenticated in external memory,
  // check freshness when swapped in.
  using NonObliviousTableType = EM::CacheFrontVector::Vector<
      BucketType, sizeof(BucketType),
      EM::CacheFrontVector::EncryptType::ENCRYPT_AND_AUTH_FRESH, 1024>;
  using TableType = std::conditional_t<isOblivious, ObliviousTableType,
                                       NonObliviousTableType>;
  // the two hash tables
  TableType table0, table1;
  // the indexer to hash the key
  Indexer indexer;

  using StashType = LRUStash<16>;
  // the stash for elements that cannot be stored in the hash tables due to
  // bucket size limitation
  StashType stash;

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
   * @brief Return a random bucket offset to replace
   *
   * @return int
   */
  static int getRandBucketOffset() {
    return (int)(UniformRandom32() % bucketSize);
  }

  /**
   * @brief A helper function that replaces an entry in a bucket if it exists,
   * and change the input entry to the old entry. The function is not oblivious.
   *
   * @param bucket the bucket to perform replace
   * @param entryToInsert the entry to insert
   * @return true if the entry is replaced, false otherwise
   */
  static bool replaceIfExist(BucketType& bucket, OPosMapEntry& entryToInsert) {
    bool updated = !entryToInsert.valid();
    for (short i = 0; i < bucketSize; ++i) {
      OPosMapEntry& entry = bucket.entries[i];
      bool matchFlag = (entry.otherIdx == entryToInsert.otherIdx) &&
                       (entry.keyHash == entryToInsert.keyHash) && (!updated);
      if (matchFlag) {
        std::swap(entry, entryToInsert);
        return true;
      }
    }
    return updated;
  }

  /**
   * @brief A helper function that replaces an entry in a bucket if it exists,
   * and change the input entry to the old entry. The function is oblivious.
   *
   * @param bucket the bucket to perform replace
   * @param entryToInsert the entry to insert
   * @return true if the entry is replaced, false otherwise
   */
  static bool replaceIfExistOblivious(BucketType& bucket,
                                      OPosMapEntry& entryToInsert) {
    bool updated = !entryToInsert.valid();
    for (short i = 0; i < bucketSize; ++i) {
      OPosMapEntry& entry = bucket.entries[i];
      bool matchFlag = (entry.otherIdx == entryToInsert.otherIdx) &
                       (entry.keyHash == entryToInsert.keyHash) & (!updated);
      obliSwap(matchFlag, entry, entryToInsert);
      updated |= matchFlag;
    }
    return updated;
  }

  /**
   * @brief A helper function that replaces an entry in the stash if it exists,
   * and change the input entry to the old entry. The function is not oblivious.
   *
   * @param stash the stash to perform replace
   * @param entryToInsert the entry to insert
   * @return true if the entry is replaced, false otherwise
   */
  static bool replaceIfExist(StashType& stash, OPosMapEntry& entryToInsert,
                             const PositionType& idx0) {
    bool updated = !entryToInsert.valid();
    for (size_t i = 0; i < stash.size(); ++i) {
      auto& entry = stash[i].entry;
      bool matchFlag = (entry.otherIdx == entryToInsert.otherIdx) &&
                       (stash[i].idx0 == idx0) &&
                       (entry.keyHash == entryToInsert.keyHash) && (!updated);
      if (matchFlag) {
        std::swap(entry, entryToInsert);
        return true;
      }
    }
    return updated;
  }

  /**
   * @brief A helper function that replaces an entry in the stash if it exists,
   * and change the input entry to the old entry. The function is oblivious.
   *
   * @param stash the stash to perform replace
   * @param entryToInsert the entry to insert
   * @return true if the entry is replaced, false otherwise
   */
  static bool replaceIfExistOblivious(StashType& stash,
                                      OPosMapEntry& entryToInsert,
                                      const PositionType& idx0) {
    bool updated = !entryToInsert.valid();
    for (size_t i = 0; i < stash.size(); ++i) {
      auto& entry = stash[i].entry;
      bool matchFlag = (entry.otherIdx == entryToInsert.otherIdx) &
                       (stash[i].idx0 == idx0) &
                       (entry.keyHash == entryToInsert.keyHash) & (!updated);
      obliSwap(matchFlag, entry, entryToInsert);
      updated |= matchFlag;
    }
    return updated;
  }

  /**
   * @brief A helper function that inserts an entry to a bucket if a slot is
   * empty. The function is not oblivious.
   *
   * @param bucket the bucket to perform insert
   * @param entryToInsert the entry to insert
   * @return true if the entry is inserted, false otherwise
   */
  static bool insertIfEmpty(BucketType& bucket,
                            const OPosMapEntry& entryToInsert) {
    bool updated = !entryToInsert.valid();
    for (int i = 0; i < bucketSize; ++i) {
      auto& entry = bucket.entries[i];
      bool isEmpty = !entry.valid();
      bool insertFlag = isEmpty && (!updated);

      if (insertFlag) {
        entry = entryToInsert;
        return true;
      }
    }
    return updated;
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
                                     const OPosMapEntry& entryToInsert) {
    bool updated = !entryToInsert.valid();
    for (int i = 0; i < bucketSize; ++i) {
      auto& entry = bucket.entries[i];
      bool isEmpty = !entry.valid();
      bool insertFlag = isEmpty & (!updated);
      obliMove(insertFlag, entry, entryToInsert);
      updated |= insertFlag;
    }
    return updated;
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
  bool insertEntry(OPosMapEntry& entryToInsert, PositionType& idx0) {
    if (!entryToInsert.valid()) {
      return true;
    }
    bool exist = false;
    auto table0UpdateFunc = [&](BucketType& bucket0) {
      bool replaceSucceed = replaceIfExist(bucket0, entryToInsert);
      if (replaceSucceed) {
        entryToInsert.setInvalid(true);
        exist = true;
        return;
      }
      std::swap(entryToInsert.otherIdx, idx0);
      auto table1UpdateFunc = [&](BucketType& bucket1) {
        bool replaceSucceed = replaceIfExist(bucket1, entryToInsert);
        if (replaceSucceed) {
          entryToInsert.setInvalid(true);
          exist = true;
          return;
        }
        std::swap(entryToInsert.otherIdx, idx0);
        replaceSucceed = replaceIfExist(stash, entryToInsert, idx0);
        if (replaceSucceed) {
          entryToInsert.setInvalid(true);
          exist = true;
          return;
        }
        ++load;

        bool insertSucceed = insertIfEmptyOblivious(bucket0, entryToInsert);
        if (insertSucceed) {
          entryToInsert.setInvalid(true);
          return;
        }
        std::swap(entryToInsert.otherIdx, idx0);
        insertSucceed = insertIfEmptyOblivious(bucket1, entryToInsert);
        if (insertSucceed) {
          entryToInsert.setInvalid(true);
          return;
        }
        int offset = getRandBucketOffset();
        std::swap(bucket1.entries[offset], entryToInsert);
        std::swap(entryToInsert.otherIdx, idx0);
      };
      updateHelper(idx0, table1, table1UpdateFunc);
    };
    updateHelper(idx0, table0, table0UpdateFunc);

    return exist;
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
  bool insertEntryOblivious(OPosMapEntry& entryToInsert, PositionType& idx0) {
    bool exist = false;
    auto table0UpdateFunc = [&](BucketType& bucket0) {
      bool replaceSucceed = replaceIfExistOblivious(bucket0, entryToInsert);
      exist |= replaceSucceed;
      entryToInsert.setInvalid(replaceSucceed);
      obliSwap(entryToInsert.valid(), entryToInsert.otherIdx, idx0);
      auto table1UpdateFunc = [&](BucketType& bucket1) {
        bool replaceSucceed = replaceIfExistOblivious(bucket1, entryToInsert);
        exist |= replaceSucceed;
        entryToInsert.setInvalid(replaceSucceed);
        obliSwap(entryToInsert.valid(), entryToInsert.otherIdx, idx0);
        replaceSucceed = replaceIfExistOblivious(stash, entryToInsert, idx0);
        exist |= replaceSucceed;
        entryToInsert.setInvalid(replaceSucceed);

        bool insertSucceed = insertIfEmptyOblivious(bucket0, entryToInsert);
        entryToInsert.setInvalid(insertSucceed);
        obliSwap(entryToInsert.valid(), entryToInsert.otherIdx, idx0);
        insertSucceed = insertIfEmptyOblivious(bucket1, entryToInsert);
        entryToInsert.setInvalid(insertSucceed);
        int offset = getRandBucketOffset();
        obliSwap(entryToInsert.valid(), bucket1.entries[offset], entryToInsert);
        obliSwap(entryToInsert.valid(), entryToInsert.otherIdx, idx0);
      };
      updateHelper(idx0, table1, table1UpdateFunc);
    };
    updateHelper(idx0, table0, table0UpdateFunc);

    load += !exist;
    return exist;
  }

  /**
   * @brief Try to insert entry into table 1, if table 1 is
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
  void insertEntryRetry(OPosMapEntry& entryToInsert, PositionType& idx0,
                        int maxRetry = 1) {
    auto swapUpdateFunc = [&](BucketType& bucket) {
      int offset = getRandBucketOffset();
      bool insertSucceed = insertIfEmptyOblivious(bucket, entryToInsert);
      entryToInsert.setInvalid(insertSucceed);
      if (entryToInsert.valid()) {
        std::swap(entryToInsert, bucket.entries[offset]);
        std::swap(entryToInsert.otherIdx, idx0);
      }
    };
    for (int r = 0; r < maxRetry; ++r) {
      updateHelper(idx0, table0, swapUpdateFunc);

      if (!entryToInsert.valid()) {
        break;
      }

      updateHelper(idx0, table1, swapUpdateFunc);  // modifies entryToInsert

      if (!entryToInsert.valid()) {
        break;
      }
    }
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
  void insertEntryObliviousRetry(OPosMapEntry& entryToInsert,
                                 PositionType& idx0, int maxRetry = 1) {
    auto swapUpdateFunc = [&](BucketType& bucket) {
      int offset = getRandBucketOffset();
      bool insertSucceed = insertIfEmptyOblivious(bucket, entryToInsert);
      entryToInsert.setInvalid(insertSucceed);
      obliSwap(entryToInsert.valid(), entryToInsert, bucket.entries[offset]);
      obliSwap(entryToInsert.valid(), entryToInsert.otherIdx, idx0);
    };
    for (int r = 0; r < maxRetry; ++r) {
      updateHelper(idx0, table0, swapUpdateFunc);
      updateHelper(idx0, table1, swapUpdateFunc);  // modifies entryToInsert
    }
  }

 public:
  // result data structure for a query in a batch
  struct ValResult {
    V value;
    bool found;
  };

  OPosMap() {}

  /**
   * @brief Construct a new OPosMap object
   *
   * @param size capacity of the hash map
   */
  explicit OPosMap(PositionType size) { SetSize(size); }

  /**
   * @brief Construct a new OPosMap object
   *
   * @param size capacity of the hash map
   * @param cacheBytes the size of available cache in bytes
   */
  OPosMap(PositionType size, uint64_t cacheBytes) { SetSize(size, cacheBytes); }

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
    stash.SetSize(16);
  }

  /**
   * @brief Get the memory usage of this hash map
   *
   * @return uint64_t the heap memory usage in bytes
   */
  uint64_t GetMemoryUsage() const {
    return table0.GetMemoryUsage() + table1.GetMemoryUsage() +
           stash.size() * sizeof(OPosMapEntry);
  }

  PositionType size() const { return _size; }

  PositionType GetLoad() const { return load; }

  PositionType GetTableSize() const { return tableSize; }

  const Indexer& GetIndexer() const { return indexer; }

  void SetIndexer(const Indexer& _indexer) { indexer = _indexer; }

  TableType& GetTable0() { return table0; }

  TableType& GetTable1() { return table1; }

  const StashType& GetStash() const { return stash; }

  using NonObliviousPosMap = OPosMap<K, V, false, true>;

  /**
   * @brief Initialize the oblivious position map from another non-oblivious
   * position map. This function is more efficient than inserting elements one
   * by one.
   *
   * @param other the non-oblivious position map
   */
  void InitFromNonOblivious(NonObliviousPosMap& other) {
    static_assert(isOblivious,
                  "Only oblivious position map can call this function");
    if (_size != other.size()) {
      throw std::runtime_error(
          "OPosMap InitFromNonOblivious failed because the size of the two "
          "maps are different");
    }
    if (_size == 0) {
      throw std::runtime_error(
          "OPosMap size not set. Call SetSize before initialization.");
    }
    indexer = other.GetIndexer();

    typename NonObliviousTableType::Reader reader0(other.GetTable0().begin(),
                                                   other.GetTable0().end());
    typename NonObliviousTableType::Reader reader1(other.GetTable1().begin(),
                                                   other.GetTable1().end());
    if constexpr (parallel_init) {
#pragma omp task
      { table0.InitFromReader(reader0); }

      { table1.InitFromReader(reader1); }
#pragma omp taskwait
    } else {
      table0.InitFromReader(reader0);
      table1.InitFromReader(reader1);
    }
    load = other.GetLoad();
    for (const auto& entry : other.GetStash().stash) {
      const OPosMapEntry& posmapEntry = entry.entry;
      if (posmapEntry.valid()) {
        stash.OInsert(posmapEntry, entry.idx0);
      }
    }
  }

  /**
   * @brief Initialize an empty hash map
   *
   */
  void Init() {
    if (_size == 0) {
      throw std::runtime_error(
          "OPosMap size not set. Call SetSize before initialization.");
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
   * @brief Insert obliviously. Hide the number of replacement and whether the
   * insertion is dummy
   *
   * @param key the key to insert
   * @param value the value to insert
   * @param[out] uid outputs the unique identifier of the key
   * @param[out] extraHash outputs the extra hash value
   * @return true if the key already exists, false otherwise
   */
  bool Insert(const K& key, V& value, UidType& uid, HExtra& extraHash) {
    PositionType idx0, idx1;
    H keyHash;
    indexer.getHashIndices(key, idx0, idx1, keyHash, uid, extraHash);
    OPosMapEntry entryToInsert = {keyHash, idx1, value};
    bool exist = insertEntry(entryToInsert, idx0);
    obliMove(exist, value, entryToInsert.value);
    // the element just swapped out is more likely to get inserted to somewhere
    // else
    stash.template Insert<true>(entryToInsert, idx0);
    for (int i = 0; i < 10; ++i) {
      // use FIFO order so that we won't get stuck by loops in the random graph
      // of cuckoo hashing
      stash.PopOldest(entryToInsert, idx0);
      if (!entryToInsert.valid()) {
        break;
      }
      insertEntryRetry(entryToInsert, idx0, 1);
      stash.Insert(entryToInsert, idx0);
    }
    return exist;
  }

  /**
   * @brief Insert obliviously. Hide the number of replacement and whether the
   * insertion is dummy
   *
   * @param key the key to insert
   * @param value the value to insert
   * @param[out] uid outputs the unique identifier of the key
   * @param[out] extraHash outputs the extra hash value
   * @param isDummy whether the insertion is dummy
   * @return true if the key already exists, false otherwise
   */
  bool OInsert(const K& key, V& value, UidType& uid, HExtra& extraHash,
               bool isDummy = false) {
    PositionType idx0, idx1;
    H keyHash;
    K keyToHash = key;
    // if the underlying ram is not oblivious, generate a random key for dummy
    // entry
    if constexpr (!isOblivious) {
      K randKey;
      read_rand((uint8_t*)&randKey, sizeof(K));
      obliMove(isDummy, keyToHash, randKey);
    }
    indexer.getHashIndices(keyToHash, idx0, idx1, keyHash, uid, extraHash);
    OPosMapEntry entryToInsert = {keyHash, idx1, value};
    entryToInsert.setInvalid(isDummy);
    bool exist = insertEntryOblivious(entryToInsert, idx0);
    obliMove(exist, value, entryToInsert.value);
    // the element just swapped out is more likely to get inserted to somewhere
    // else
    stash.template OInsert<true>(entryToInsert, idx0);
    for (int i = 0; i < 2; ++i) {
      // use FIFO order so that we won't get stuck by loops in the random graph
      // of cuckoo hashing
      stash.OPopOldest(entryToInsert, idx0);
      const bool isDummy = !entryToInsert.valid();
      if constexpr (!isOblivious) {
        // the underlying ram is not oblivious, generate a random access
        // position
        obliMove(isDummy, idx0, (PositionType)UniformRandom(tableSize - 1));
        obliMove(isDummy, entryToInsert.otherIdx,
                 (PositionType)UniformRandom(tableSize - 1));
      }
      insertEntryObliviousRetry(entryToInsert, idx0, 1);
      stash.OInsert(entryToInsert, idx0);
    }
    return exist;
  }

  /**
   * @brief Update the value of a key if it exists and return the old value.
   *
   * @param key the key to search
   * @param[in, out] value the input new value and value to return
   * @param[out] uid outputs the unique identifier of the key
   * @param[out] extraHash outputs the extra hash value
   * @param isDummy whether the search is dummy
   * @return true if there is a matching entry in the position map, false
   * otherwise
   */
  bool Update(const K& key, V& value, UidType& uid, HExtra& extraHash,
              bool isDummy = false) {
    PositionType idx0, idx1;
    H keyHash;
    K keyToHash = key;
    bool found = false;
    // generate a random key for dummy entry if the underlying ram is not
    // oblivious
    if constexpr (!isOblivious) {
      K randKey;
      read_rand((uint8_t*)&randKey, sizeof(K));
      obliMove(isDummy, keyToHash, randKey);
    }
    indexer.getHashIndices(keyToHash, idx0, idx1, keyHash, uid, extraHash);
    OPosMapEntry entryToFind = {keyHash, idx1, value};
    entryToFind.setInvalid(isDummy);
    auto bucketAccessor = [&](BucketType& bucket) {
      for (int i = 0; i < bucketSize; ++i) {
        auto& entry = bucket.entries[i];
        bool matchFlag =
            (entry.keyHash == entryToFind.keyHash) & (entry.otherIdx == idx1);
        found |= matchFlag;
        obliSwap(matchFlag, entry.value, value);
      }
    };
    updateHelper(idx0, table0, bucketAccessor);
    std::swap(idx0, idx1);
    updateHelper(idx0, table1, bucketAccessor);
    std::swap(idx0, idx1);
    for (size_t i = 0; i < stash.size(); ++i) {
      auto& stashEntry = stash[i];
      auto& entry = stashEntry.entry;
      bool match = (entry.keyHash == entryToFind.keyHash) &
                   (entry.otherIdx == idx1) & stashEntry.idx0 == idx0;
      obliSwap(match, value, entry.value);
      found |= match;
    }

    return found & (!isDummy);
  }

  /**
   * @brief Erase a key from the position map. The function is not oblivious.
   * Note that there's chance that the key is not actually present even though
   * the indices and key hash match. We need to check the extra hash in the main
   * oram before erasing the entry from the position map.
   *
   * @param key the key to erase
   * @param value the new value for the entry in case the key is mismatched in
   * the position map
   * @param mainMapErase the function to erase the key in the main oram, should
   * return true if the key is found
   * @param isDummy whether the erase is dummy operation
   * @return true if there is an entry in the position map that matches the key,
   * false otherwise
   */
  template <class MainMapEraseFunc>
  bool Erase(const K& key, const V& value,
             const MainMapEraseFunc& mainMapErase) {
    PositionType idx0, idx1;
    UidType uid;
    HExtra extraHash;
    H keyHash;
    indexer.getHashIndices(key, idx0, idx1, keyHash, uid, extraHash);
    bool erased = false;
    OPosMapEntry entryToErase = {keyHash, idx1, DUMMY<V>()};
    auto eraseTable0Func = [&](BucketType& bucket0) {
      PositionType entryToErasePos = value;
      for (int i = 0; i < bucketSize; ++i) {
        OPosMapEntry& entry = bucket0.entries[i];
        bool matchFlag =
            (entry.keyHash == entryToErase.keyHash) && (entry.otherIdx == idx1);
        if (matchFlag) {
          std::swap(entry.value, entryToErasePos);
          erased = mainMapErase(entryToErasePos, uid, extraHash);
          entry.setInvalid(erased);
          return;
        }
      }
      auto eraseTable1Func = [&](BucketType& bucket1) {
        for (int i = 0; i < bucketSize; ++i) {
          OPosMapEntry& entry = bucket1.entries[i];
          bool matchFlag = (entry.keyHash == entryToErase.keyHash) &&
                           (entry.otherIdx == idx0);
          if (matchFlag) {
            std::swap(entry.value, entryToErasePos);
            erased = mainMapErase(entryToErasePos, uid, extraHash);
            entry.setInvalid(erased);
            return;
          }
        }

        for (size_t i = 0; i < stash.size(); ++i) {
          auto& stashEntry = stash[i];
          OPosMapEntry& entry = stashEntry.entry;
          bool match = (entry.keyHash == entryToErase.keyHash) &&
                       (entry.otherIdx == idx1) && (stashEntry.idx0 == idx0);
          if (match) {
            std::swap(entry.value, entryToErasePos);
            erased = mainMapErase(entryToErasePos, uid, extraHash);
            entry.setInvalid(erased);
            return;
          }
        }
      };
      updateHelper(idx1, table1, eraseTable1Func);
    };

    updateHelper(idx0, table0, eraseTable0Func);
    load -= erased;
    return erased;
  }

  /**
   * @brief Erase a key from the position map. The function is oblivious. Note
   * that there's chance that the key is not actually present even though the
   * indices and key hash match. We need to check the extra hash in the main
   * oram before erasing the entry from the position map.
   *
   * @param key the key to erase
   * @param value the new value for the entry in case the key is mismatched in
   * the position map
   * @param mainMapErase the function to erase the key in the main oram, should
   * return true if the key is found
   * @param isDummy whether the erase is dummy operation
   * @return true if there is an entry in the position map that matches the key,
   * false otherwise
   */
  template <class MainMapEraseFunc>
  bool OErase(const K& key, const V& value,
              const MainMapEraseFunc& mainMapErase, bool isDummy = false) {
    PositionType idx0, idx1;
    UidType uid;
    HExtra extraHash;
    H keyHash;
    K keyToHash = key;
    if constexpr (!isOblivious) {
      K randKey;
      read_rand((uint8_t*)&randKey, sizeof(K));
      obliMove(isDummy, keyToHash, randKey);
    }
    indexer.getHashIndices(keyToHash, idx0, idx1, keyHash, uid, extraHash);
    bool erased = false;
    OPosMapEntry entryToErase = {keyHash, idx1, DUMMY<V>()};
    entryToErase.setInvalid(isDummy);
    auto eraseTable0Func = [&](BucketType& bucket0) {
      PositionType entryToErasePos = value;
      int eraseTable0Idx = -1;
      for (int i = 0; i < bucketSize; ++i) {
        OPosMapEntry& entry = bucket0.entries[i];
        bool matchFlag =
            (entry.keyHash == entryToErase.keyHash) & (entry.otherIdx == idx1);
        obliSwap(matchFlag, entryToErasePos, entry.value);
        obliMove(matchFlag, eraseTable0Idx, i);
      }
      auto eraseTable1Func = [&](BucketType& bucket1) {
        int eraseTable1Idx = -1;
        for (int i = 0; i < bucketSize; ++i) {
          OPosMapEntry& entry = bucket1.entries[i];
          bool matchFlag = (entry.keyHash == entryToErase.keyHash) &
                           (entry.otherIdx == idx0);
          obliSwap(matchFlag, entryToErasePos, entry.value);
          obliMove(matchFlag, eraseTable1Idx, i);
        }

        size_t eraseStashIdx = -1;
        for (size_t i = 0; i < stash.size(); ++i) {
          auto& stashEntry = stash[i];
          OPosMapEntry& entry = stashEntry.entry;
          bool match = (entry.keyHash == entryToErase.keyHash) &
                       (entry.otherIdx == idx1) & (stashEntry.idx0 == idx0);
          obliSwap(match, entryToErasePos, entry.value);
          obliMove(match, eraseStashIdx, i);
        }
        bool notFoundFlag = (eraseTable0Idx == -1) & (eraseTable1Idx == -1) &
                            (eraseStashIdx == -1);
        obliMove(notFoundFlag | isDummy, uid, DUMMY<UidType>());
        // we check the main oram to ensure that the element is actually present
        // (rather than a false positive) before erasing it
        erased = mainMapErase(entryToErasePos, uid, extraHash);
        for (int i = 0; i < bucketSize; ++i) {
          auto& entry = bucket0.entries[i];
          entry.setInvalid(erased & (i == eraseTable0Idx));
        }
        for (int i = 0; i < bucketSize; ++i) {
          auto& entry = bucket1.entries[i];
          entry.setInvalid(erased & (i == eraseTable1Idx));
        }
        for (size_t i = 0; i < stash.size(); ++i) {
          auto& stashEntry = stash[i];
          stashEntry.entry.setInvalid(erased & (i == eraseStashIdx));
        }
      };
      updateHelper(idx1, table1, eraseTable1Func);
    };

    updateHelper(idx0, table0, eraseTable0Func);
    load -= erased;
    return erased;
  }

};  // class OPosMap

template <typename K, typename V, typename PositionType = uint64_t>
struct OMap {
  using PosMapType = OPosMap<K, PositionType, true>;
  PosMapType keyPosMap;
  bool inited = false;
  using UidType = typename PosMapType::UidType;
  using HExtra = typename PosMapType::HExtra;
  struct ORAMEntry {
    HExtra hash;
    V value;
#ifndef ENCLAVE_MODE
    friend std::ostream& operator<<(std::ostream& os, const ORAMEntry& pair) {
      os << "(key hash: " << pair.hash << " value: " << pair.value << ")";
      return os;
    }
#endif
  };

  struct StashEntry {
    bool valid = false;
    K key;
    V value;
  };

  std::vector<StashEntry> stash;
  // additional stash to put elements with duplicate
  // info in position map but has different keys

  // TODO: freshness check if swap is needed
  CircuitORAM::ORAM<ORAMEntry, 2, 20, PositionType, UidType, 4096, false> oram;

  OMap() {}
  OMap(PositionType size) { SetSize(size); }
  OMap(PositionType size, uint64_t cacheBytes) { SetSize(size, cacheBytes); }

  void SetSize(PositionType size) {
    keyPosMap.SetSize(size);
    oram.SetSize(size);
    stash.resize(4);
  }

  void SetSize(PositionType size, uint64_t cacheBytes) {
    double posMapCacheRatio = 0.7;
    keyPosMap.SetSize(size, (uint64_t)(cacheBytes * posMapCacheRatio));
    cacheBytes -= keyPosMap.GetMemoryUsage();
    oram.SetSize(size, cacheBytes);
    stash.resize(4);
  }

  void Init() {
    if (inited) {
      throw std::runtime_error("OPosMap initialized twice.");
    }
    inited = true;
    keyPosMap.Init();
  }

  PositionType size() const { return keyPosMap.size(); }

  /**
   * @brief An object that stores the curret state of initialization.
   Faciliates
   * initialization in a streaming fashion.
   *
   */
  struct InitContext {
    using NonObliviousPosMap = OPosMap<K, PositionType, false, true>;

   private:
    NonObliviousPosMap* nonObliviousPosMap;
    OMap& omap;  // the parent map

   public:
    /**
     * @brief Construct a new InitContext object
     *
     * @param omap The parent map
     * @param additionalCacheBytes the size of additional available cache
     for
     * initialization
     */
    explicit InitContext(OMap& map, uint64_t additionalCacheBytes = 0)
        : omap(map) {
      if (map.inited) {
        throw std::runtime_error("OPosMap initialized twice.");
      }
      map.inited = true;
      if (map.size() == 0) {
        throw std::runtime_error(
            "OPosMap size not set. Call SetSize before initialization.");
      }
      if constexpr (omap.keyPosMap.isObliviousPosMap) {
        nonObliviousPosMap =
            new NonObliviousPosMap(map.size(), additionalCacheBytes);
      }
      map.keyPosMap.SetIndexer(nonObliviousPosMap->GetIndexer());
    }

    /**
     * @brief Move constructor
     * @param other the other InitContext
     *
     */
    InitContext(const InitContext&& other)
        : omap(other.omap), nonObliviousPosMap(other.nonObliviousPosMap) {}

    InitContext(const InitContext& other) = delete;

    /**
     * @brief Insert a new key value pair for initialization. The method
     will
     * reveal whether the key has been inserted before. But if all the keys
     are
     * distinct, the operation is oblivious. The method will throw an
     exception
     * if too many keys are inserted.
     *
     * @param key
     * @param value
     */
    void Insert(const K& key, const V& value) {
      if constexpr (!omap.keyPosMap.isObliviousPosMap) {
        omap.Insert(key, value);
      } else {
        omap.InsertWithCustomPosMap(key, value, *nonObliviousPosMap);
      }
    }

    void Insert(const std::pair<K, V>& entry) {
      Insert(entry.first, entry.second);
    }

    template <class Iterator>
    void InsertBatch(Iterator begin, Iterator end) {
      for (auto it = begin; it != end; ++it) {
        Insert(it->first, it->second);
      }
    }

    /**
     * @brief Finalize the initialization. The method will copy the data
     from
     * the non-oblivious hash map to the oblivious hash map.
     *
     */
    void Finalize() {
      if constexpr (omap.keyPosMap.isObliviousPosMap) {
        omap.keyPosMap.InitFromNonOblivious(*nonObliviousPosMap);
        delete nonObliviousPosMap;
      }
    }
  };

  /**
 * @brief Obtain a new context to initialize this map. The initialization
 data
 * can be either private or public (the initialization and subsequent
 accesses
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
    return new InitContext(*this, additionalCacheBytes);
  }

  template <typename Reader>
    requires Readable<Reader, std::pair<K, V>>
  void InitFromReader(Reader& reader) {
    InitContext initContext(*this);
    while (!reader.eof()) {
      initContext.Insert(reader.read());
    }
    initContext.Finalize();
  }

  template <typename Reader>
    requires Readable<Reader, std::pair<K, V>>
  void InitFromReader(Reader& reader, uint64_t additionalCacheBytes) {
    InitContext initContext(*this, additionalCacheBytes);
    while (!reader.eof()) {
      initContext.Insert(reader.read());
    }
    initContext.Finalize();
  }

  template <class PosMap>
  bool InsertWithCustomPosMap(const K& key, const V& value, PosMap& posMap) {
    // first pass, check if key exists in stash, if so, replace the value and we
    // are done
    for (StashEntry& pair : stash) {
      bool matchFlag = pair.key == key && pair.valid;
      if (matchFlag) {
        pair.value = value;
        return true;
      }
    }
    PositionType newPos = oram.GetRandPos();
    PositionType pos = newPos;
    UidType uid;
    HExtra extraHash;
    bool existInCuckoo = posMap.Insert(key, pos, uid, extraHash);
    if (!existInCuckoo) {
      pos = oram.GetRandPos();
    }

    bool insertToStashFlag = false;
    oram.Update(pos, uid, newPos, [&](ORAMEntry& prevPair) {
      bool sameKey = prevPair.hash == extraHash;
      insertToStashFlag = existInCuckoo && !sameKey;
      bool modifyPairFlag = !insertToStashFlag;
      existInCuckoo = existInCuckoo && sameKey;
      if (modifyPairFlag) {
        prevPair = {extraHash, value};
      }
      return true;
    });
    // case 1: exist and same key -> don't change stash
    // case 2: exist and not same key -> put to the stash (if the key exists in
    // stash, replace, otherwise, insert)
    // case 3: not exist -> if the key exists in stash, remove it and set exist
    // flag to true, otherwise don't change the stash

    // second pass, insert to stash if needed
    if (insertToStashFlag) {
      for (StashEntry& pair : stash) {
        bool insertFlag = !pair.valid;
        if (insertFlag) {
          pair = {true, key, value};
          insertToStashFlag = false;
          break;
        }
      }
      if (insertToStashFlag) {
        PERFCTR_INCREMENT(OHMAP_DEAMORT_OVERFLOW);
        stash.push_back({true, key, value});
      }
    }

    return existInCuckoo;
  }

  /**
   * @brief Insert a key value pair into the map, return true if the key is
   * already in the map, in which case the value is updated. The function is
   * not oblivious and reveals whether the key exists. It is faster than
   * OInsert, and should be used when the database is public.
   *
   */
  bool Insert(const K& key, const V& value) {
    return InsertWithCustomPosMap(key, value, keyPosMap);
  }

  bool OInsert(const K& key, const V& value) {
    PositionType newPos = oram.GetRandPos();
    PositionType pos = newPos;
    UidType uid;
    HExtra extraHash;
    bool existInCuckoo = keyPosMap.OInsert(key, pos, uid, extraHash, false);

    obliMove(!existInCuckoo, pos, oram.GetRandPos());
    bool existInStash = false;

    // first pass, check if key exists in stash, if so, replace the value and we
    // are done
    for (StashEntry& pair : stash) {
      bool matchFlag = pair.valid & (pair.key == key);
      obliMove(matchFlag, pair.value, value);
      existInStash |= matchFlag;
    }
    bool insertToStashFlag = false;
    bool updateKeepFlag = existInCuckoo | !existInStash;
    oram.Update(pos, uid, newPos, [&](ORAMEntry& prevPair) {
      bool sameKey = prevPair.hash == extraHash;
      insertToStashFlag = existInCuckoo & !sameKey & !existInStash;
      bool modifyPairFlag = !existInStash & !insertToStashFlag;
      existInCuckoo &= sameKey;
      obliMove(modifyPairFlag, prevPair, {extraHash, value});

      return updateKeepFlag;
    });
    // case 1: exist and same key -> don't change stash
    // case 2: exist and not same key -> put to the stash (if the key exists in
    // stash, replace, otherwise, insert)
    // case 3: not exist -> if the key exists in stash, remove it and set exist
    // flag to true, otherwise don't change the stash

    // second pass, insert to stash if needed
    for (StashEntry& pair : stash) {
      bool insertFlag = !pair.valid & insertToStashFlag;
      obliMove(insertFlag, pair, {true, key, value});
      insertToStashFlag &= !insertFlag;
    }
    if (insertToStashFlag) {
      PERFCTR_INCREMENT(OHMAP_DEAMORT_OVERFLOW);
      stash.push_back({true, key, value});
    }
    return existInCuckoo | existInStash;
  }

  /**
   * @brief Find the value of a key in the map. The function is oblivious.
   *
   * @param key the key to search
   * @param[out] value the value to return
   * @return true if the key is found, false otherwise
   */
  bool Find(const K& key, V& value) {
    PositionType newPos = oram.GetRandPos();
    PositionType pos = newPos;
    UidType uid;
    HExtra extraHash;

    bool exist = keyPosMap.Update(key, pos, uid, extraHash, false);
    // if (exist) {
    //   std::cout << "Found position " << pos << " for key " << key << "\n";
    // } else {
    //   std::cout << "Doesn't exist in cuckoo\n";
    // }
    // std::cout << "Write to new position " << newPos << "\n";
    ORAMEntry pair;

    obliMove(!exist, uid, DUMMY<UidType>());
    oram.Read(pos, uid, pair, newPos);
    exist &= pair.hash == extraHash;
    value = pair.value;

    for (StashEntry& pair : stash) {
      bool matchFlag = pair.valid & (pair.key == key);
      // if (matchFlag) {
      //   std::cout << "find key " << key << " in stash\n";
      // }
      obliMove(matchFlag, value, pair.value);
      exist |= matchFlag;
    }
    // oram.printState();
    return exist;
  }

  /**
   * @brief Erase a key from the map. The function is not oblivious and reveals
   * whether the key exists in the map. The function is faster than OErase and
   * should be used when the database is public.
   *
   * @param key the key to erase
   * @return true if there is an entry in the position map that matches the key,
   * false otherwise
   */
  bool Erase(const K& key) {
    PositionType newPos = oram.GetRandPos();  // possibly the erasure is dummy
    auto mainMapErase = [&](PositionType pos, UidType uid,
                            const HExtra& extraHash) {
      bool exist = false;
      oram.Update(pos, uid, newPos, [&](ORAMEntry& prevPair) {
        bool sameKey = prevPair.hash == extraHash;
        exist = sameKey;
        // return false to erase the entry
        return !sameKey;
      });

      return exist;
    };

    bool erased = keyPosMap.Erase(key, newPos, mainMapErase);
    if (erased) {
      return true;
    }
    for (StashEntry& pair : stash) {
      bool matchFlag = pair.valid && (pair.key == key);
      if (matchFlag) {
        pair.valid = false;
        return true;
      }
    }
    return false;
  }

  /**
   * @brief Erase a key from the map. The function is oblivious.
   *
   * @param key the key to erase
   * @return true if there is an entry in the position map that matches the key,
   * false otherwise
   */
  bool OErase(const K& key) {
    PositionType newPos = oram.GetRandPos();  // possibly the erasure is dummy
    auto mainMapErase = [&](PositionType pos, UidType uid,
                            const HExtra& extraHash) {
      // if the position map doesn't find the key, perform a dummy erase
      bool dummyFlag = uid == DUMMY<UidType>();
      obliMove(dummyFlag, pos, oram.GetRandPos());
      bool exist = false;
      oram.Update(pos, uid, newPos, [&](ORAMEntry& prevPair) {
        bool sameKey = prevPair.hash == extraHash;
        exist = sameKey;
        // return false to erase the entry
        return !sameKey;
      });

      exist &= !dummyFlag;
      return exist;
    };

    bool erased = keyPosMap.OErase(key, newPos, mainMapErase, false);

    for (StashEntry& pair : stash) {
      bool matchFlag = pair.valid & (pair.key == key);
      pair.valid &= !matchFlag;
      erased |= matchFlag;
    }
    return erased;
  }

  /**
   * @brief Get the number of elements in the map
   */
  PositionType GetLoad() const { return keyPosMap.GetLoad(); }
};
}  // namespace ODSL