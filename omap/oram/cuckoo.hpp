#pragma once
#include "recursive_oram.hpp"
namespace ODSL {

template <typename K, typename PositionType = uint64_t>
struct CuckooHashMapIndexer {
  static constexpr int saltLength = 16;
  uint8_t salts[saltLength];
  PositionType _size;

  CuckooHashMapIndexer() {}

  CuckooHashMapIndexer(PositionType size) { SetSize(size); }

  void SetSize(PositionType size) {
    _size = size;
    read_rand(salts, saltLength);
  }

  struct HashIndices {
    uint64_t h0;
    uint64_t h1;
  };

  void getHashIndices(const K& key, PositionType& pos0,
                      PositionType& pos1) const {
    HashIndices hashIndices;
    secure_hash_with_salt(key, salts, (uint8_t*)&hashIndices,
                          sizeof(hashIndices));
    pos0 = hashIndices.h0 % _size;
    pos1 = hashIndices.h1 % _size;
  }

  PositionType getHashIdx0(const K& key) const {
    PositionType pos0, pos1;
    getHashIndices(key, pos0, pos1);
    return pos0;
  }

  PositionType getHashIdx1(const K& key) const {
    PositionType pos0, pos1;
    getHashIndices(key, pos0, pos1);
    return pos1;
  }
};

template <typename K, typename V>
struct CuckooHashMapEntry {
  bool valid = false;
  bool dummy = false;
  K key;
  V value;
#ifndef ENCLAVE_MODE
  // cout
  friend std::ostream& operator<<(std::ostream& os,
                                  const CuckooHashMapEntry& entry) {
    os << "(" << entry.valid << ", " << entry.key << ", " << entry.value << ")";
    return os;
  }
#endif
};

template <typename K, typename V, const short bucketSize>
struct CuckooHashMapBucket {
  CuckooHashMapEntry<K, V> entries[bucketSize];
#ifndef ENCLAVE_MODE
  // cout
  friend std::ostream& operator<<(std::ostream& os,
                                  const CuckooHashMapBucket& bucket) {
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

template <typename K, typename V, const bool isOblivious,
          typename PositionType = uint64_t, const bool parallel_init = true>
struct CuckooHashMap {
  static constexpr int saltLength = 16;
  static constexpr double loadFactor = 0.7;
  static constexpr short bucketSize = 2;
  static constexpr int stash_max_size = 10;
  PositionType _size;
  PositionType load;
  PositionType tableSize;
  using BucketType = CuckooHashMapBucket<K, V, bucketSize>;
  using TableType =
      std::conditional_t<isOblivious, RecursiveORAM<BucketType, PositionType>,
                         StdVector<BucketType>>;
  // EM::CacheFrontVector::Vector<BucketType, sizeof(BucketType), true, true,
  //  1024>>;
  TableType table0, table1;
  CuckooHashMapIndexer<K, PositionType> indexer;
  using KVEntry = CuckooHashMapEntry<K, V>;
  std::vector<KVEntry> stash;

  struct ValResult {
    V value;
    bool found;
  };

  CuckooHashMap() {}

  CuckooHashMap(PositionType size) { SetSize(size); }

  CuckooHashMap(PositionType size, uint64_t cacheBytes) {
    SetSize(size, cacheBytes);
  }

  CuckooHashMap(PositionType size, uint64_t cacheBytes, int maxThreads) {
    SetSize(size, cacheBytes, maxThreads);
  }

  void SetSize(PositionType size) { SetSize(size, DEFAULT_HEAP_SIZE); }

  void SetSize(PositionType size, uint64_t cacheBytes) {
    _size = size;
    load = 0;
    tableSize = size / (2 * loadFactor * bucketSize) + 1;
    indexer.SetSize(tableSize);
    table0.SetSize(tableSize, cacheBytes / 2);
    table1.SetSize(tableSize, cacheBytes / 2);
  }

  void SetIndexer(const CuckooHashMapIndexer<K, PositionType>& indexer) {
    this->indexer = indexer;
  }

  using NonObliviousCuckooHashMap = CuckooHashMap<K, V, false, PositionType>;

  void InitFromNonOblivious(NonObliviousCuckooHashMap& other) {
    static_assert(isOblivious);
    if (tableSize != other.tableSize || _size != other._size) {
      throw std::runtime_error("CuckooHashMap InitFromNonOblivious failed");
    }
    load = other.load;
    // stash = other.stash;
    indexer = other.indexer;
    // change dummies to invalid
    EM::VirtualVector::VirtualReader<BucketType> reader0(
        other.tableSize, [&](PositionType i) {
          BucketType bucket = other.table0[i];
          for (int j = 0; j < bucketSize; ++j) {
            bucket.entries[j].valid &= !bucket.entries[j].dummy;
            bucket.entries[j].dummy = false;
          }
          return bucket;
        });
    EM::VirtualVector::VirtualReader<BucketType> reader1(
        other.tableSize, [&](PositionType i) {
          BucketType bucket = other.table1[i];
          for (int j = 0; j < bucketSize; ++j) {
            bucket.entries[j].valid &= !bucket.entries[j].dummy;
            bucket.entries[j].dummy = false;
          }
          return bucket;
        });
    if constexpr (parallel_init) {
      // seems to have concurrency bug
#pragma omp task
      { table0.InitFromReaderInPlace(reader0); }

      { table1.InitFromReaderInPlace(reader1); }
#pragma omp taskwait
    } else {
      table0.InitFromReaderInPlace(reader0);
      table1.InitFromReaderInPlace(reader1);
    }
    for (const auto& entry : other.stash) {
      if (entry.valid) {
        // valid is public but dummy is not
        // since we deleted all dummies, it's likely that we don't need to put
        // these entries in the stash
        insertOblivious(entry.key, entry.value, entry.dummy);
      }
    }
  }

  template <typename Reader>
  void InitFromReaderInPlace(Reader& reader) {
    if constexpr (!isOblivious) {
      while (!reader.eof()) {
        std::pair<K, V> entry = reader.read();
        insert(entry.first, entry.second);
      }
    } else {
      NonObliviousCuckooHashMap nonObliviousCuckooHashMap(_size, 0);
      nonObliviousCuckooHashMap.InitFromReaderInPlace(reader);
      InitFromNonOblivious(nonObliviousCuckooHashMap);
    }
  }

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

  void Init() { InitDefault(); }

  static void updateHelper(PositionType addr, TableType& table,
                           const auto& updateFunc) {
    if constexpr (isOblivious) {
      table.Access(addr, updateFunc);
    } else {
      updateFunc(table[addr]);
    }
  }

  static void readHelper(PositionType addr, TableType& table,
                         BucketType& outputBucket) {
    if constexpr (isOblivious) {
      table.Read(addr, outputBucket);
    } else {
      outputBucket = table[addr];
    }
  }

  template <const bool hideDummy = false>
  static bool replaceIfExist(BucketType& bucket, const KVEntry& entryToInsert) {
    for (int i = 0; i < bucketSize; ++i) {
      auto& entry = bucket.entries[i];
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

  template <const bool hideDummy = false>
  static bool replaceIfExist(std::vector<KVEntry>& stash,
                             const KVEntry& entryToInsert) {
    for (KVEntry& entry : stash) {
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

  static bool insertToStash(const KVEntry& entryToInsert,
                            std::vector<KVEntry>& stash) {
    if (stash.size() < stash_max_size) {
      stash.push_back(entryToInsert);
      return true;
    }
    return false;
  }

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

  static bool replaceIfExistOblivious(BucketType& bucket,
                                      const KVEntry& entryToInsert) {
    bool updated = !entryToInsert.valid;
    for (int i = 0; i < bucketSize; ++i) {
      auto& entry = bucket.entries[i];
      bool matchFlag = entry.valid & entry.key == entryToInsert.key & !updated;
      obliMove(matchFlag, entry, entryToInsert);
      updated |= matchFlag;
    }
    return updated;
  }

  static bool replaceIfExistOblivious(std::vector<KVEntry>& stash,
                                      const KVEntry& entryToInsert) {
    bool updated = !entryToInsert.valid;
    for (KVEntry& entry : stash) {
      bool matchFlag = entry.valid & entry.key == entryToInsert.key & !updated;
      obliMove(matchFlag, entry, entryToInsert);
      updated |= matchFlag;
    }
    return updated;
  }

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
        return;
      }
      PositionType idx1 = indexer.getHashIdx1(entryToInsert.key);
      auto table1UpdateFunc = [&](BucketType& bucket1) {
        if (replaceIfExist<hideDummy>(bucket1, entryToInsert)) {
          inserted = true;
          exist = true;
          return;
        }
        if (replaceIfExist<hideDummy>(stash, entryToInsert)) {
          exist = true;
          inserted = true;
        }
        if (insertIfEmpty(bucket0, entryToInsert)) {
          inserted = true;
          return;
        }
        if (insertIfEmpty(bucket1, entryToInsert)) {
          inserted = true;
          return;
        }
        std::swap(entryToInsert, bucket0.entries[0]);
      };
      updateHelper(idx1, table1, table1UpdateFunc);
    };
    updateHelper(idx0, table0, table0UpdateFunc);

    load += !exist;
    if (inserted) {
      return exist;
    }
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
    printf("Warning CuckooHashMap uses stash\n");
    if (!insertToStash(entryToInsert, stash)) {
      throw std::runtime_error("CuckooHashMap insert failed");
    }
    return false;
  }

  // hide the number of swaps in cuckoo hash table and whether the insertion is
  // dummy
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
        // obliMove(!entryToInsert.valid, entryToInsert.key,
        //          bucket0.entries[0].key);
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
      throw std::runtime_error("CuckooHashMap insert failed");
    }
    return exist;
  }

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

  static void mergeValRes(const ValResult& valRes0, const ValResult& valRes1,
                          ValResult& valRes) {
    valRes.found = valRes0.found | valRes1.found;
    valRes.value = valRes0.value;
    obliMove(valRes1.found, valRes.value, valRes1.value);
  }

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

  void writeBackTable(int tableNum) {
    static_assert(isOblivious);
    Assert(tableNum == 0 || tableNum == 1);
    if (tableNum == 0) {
      table0.WriteBack();
    } else {
      table1.WriteBack();
    }
  }

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

};  // namespace ODSL
}  // namespace ODSL