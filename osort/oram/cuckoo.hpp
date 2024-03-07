#pragma once
#include "recursive_oram.hpp"
namespace ODSL {

template <typename K, typename PositionType = uint64_t>
struct CuckooHashMapIndexer {
  static constexpr int saltLength = 16;
  uint8_t salts[2][saltLength];
  PositionType _size;

  CuckooHashMapIndexer() {}

  CuckooHashMapIndexer(PositionType size) { SetSize(size); }

  void SetSize(PositionType size) {
    _size = size;
    for (int i = 0; i < 2; ++i) {
      read_rand(salts[i], saltLength);
    }
  }

  PositionType getHashIdx0(const K& key) const {
    uint64_t h0 = secure_hash_with_salt((uint8_t*)&key, sizeof(K), salts[0]);
    return h0 % _size;
  }

  PositionType getHashIdx1(const K& key) const {
    uint64_t h1 = secure_hash_with_salt((uint8_t*)&key, sizeof(K), salts[1]);
    return h1 % _size;
  }
};

template <typename K, typename V>
struct CuckooHashMapEntry {
  bool valid = false;
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
          typename PositionType = uint64_t>
struct CuckooHashMap {
  static constexpr int saltLength = 16;
  static constexpr double loadFactor = 0.7;
  static constexpr short bucketSize = 2;
  PositionType _size;
  PositionType load;
  PositionType tableSize;
  using BucketType = CuckooHashMapBucket<K, V, bucketSize>;
  using TableType =
      std::conditional_t<isOblivious, RecursiveORAM<BucketType, PositionType>,
                         EM::CacheFrontVector::Vector<
                             BucketType, sizeof(BucketType), true, true, 1024>>;
  TableType table0, table1;
  CuckooHashMapIndexer<K, PositionType> indexer;
  std::vector<CuckooHashMapEntry<K, V>> stash;

  CuckooHashMap() {}

  CuckooHashMap(PositionType size) { SetSize(size); }

  CuckooHashMap(PositionType size, uint64_t cacheBytes) {
    SetSize(size, cacheBytes);
  }

  void SetSize(PositionType size,
               uint64_t cacheBytes = ((uint64_t)ENCLAVE_SIZE << 20) * 3UL /
                                     4UL) {
    _size = size;
    load = 0;
    tableSize = size / (2 * loadFactor * bucketSize);
    indexer.SetSize(tableSize);
    table0.SetSize(tableSize, cacheBytes / 2);
    table1.SetSize(tableSize, cacheBytes / 2);
  }
  using NonObliviousCuckooHashMap = CuckooHashMap<K, V, false, PositionType>;

  void InitFromNonOblivious(NonObliviousCuckooHashMap& other) {
    static_assert(isOblivious);
    if (tableSize != other.tableSize || _size != other._size) {
      throw std::runtime_error("CuckooHashMap InitFromNonOblivious failed");
    }

    load = other.load;
    stash = other.stash;
    indexer = other.indexer;
    table0.InitFromVector(other.table0);
    table1.InitFromVector(other.table1);
  }

  template <typename Reader>
  void InitFromReaderInPlace(Reader& reader) {
    if constexpr (!isOblivious) {
      while (!reader.eof()) {
        std::pair<K, V> entry = reader.read();
        insert(entry.first, entry.second);
      }
      /**
      std::cout << "CuckooHashMap load: " << load << std::endl;
      std::cout << "CuckooHashMap table0: " << std::endl;
      for (int i = 0; i < tableSize; ++i) {
        std::cout << table0[i] << std::endl;
      }
      std::cout << "CuckooHashMap table1: " << std::endl;
      for (int i = 0; i < tableSize; ++i) {
        std::cout << table1[i] << std::endl;
      }
      std::cout << "CuckooHashMap stash: " << std::endl;
      for (const auto& entry : stash) {
        std::cout << entry << " ";
      }
      std::cout << std::endl;
      */
    } else {
      NonObliviousCuckooHashMap nonObliviousCuckooHashMap(_size, 0);
      nonObliviousCuckooHashMap.InitFromReaderInPlace(reader);
      InitFromNonOblivious(nonObliviousCuckooHashMap);
    }
  }

  void InitDefault() {
    if constexpr (isOblivious) {
      table0.InitDefault(BucketType());
      table1.InitDefault(BucketType());
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

  static bool replaceIfExist(BucketType& bucket, const K& key, const V& value) {
    for (int i = 0; i < bucketSize; ++i) {
      auto& entry = bucket.entries[i];
      if (entry.valid && entry.key == key) {
        entry = {true, key, value};
        return true;
      }
    }
    return false;
  }

  static bool insertIfEmpty(BucketType& bucket, const K& key, const V& value) {
    for (int i = 0; i < bucketSize; ++i) {
      auto& entry = bucket.entries[i];
      if (!entry.valid) {
        entry = {true, key, value};
        return true;
      }
    }
    return false;
  }

  bool insert(const K& key, const V& value, bool isDummy = false) {
    if (isDummy) {
      return false;
    }
    PositionType idx0 = indexer.getHashIdx0(key);

    bool inserted = false;
    bool exist = false;
    K k = key;
    V v = value;
    auto table0UpdateFunc = [&](BucketType& bucket0) {
      if (replaceIfExist(bucket0, k, v)) {
        inserted = true;
        exist = true;
        return;
      }
      PositionType idx1 = indexer.getHashIdx1(key);
      auto table1UpdateFunc = [&](BucketType& bucket1) {
        if (replaceIfExist(bucket1, k, v)) {
          inserted = true;
          exist = true;
          return;
        }
        if (insertIfEmpty(bucket0, k, v)) {
          inserted = true;
          return;
        }
        if (insertIfEmpty(bucket1, k, v)) {
          inserted = true;
          return;
        }
        std::swap(k, bucket0.entries[0].key);
        std::swap(v, bucket0.entries[0].value);
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
      if (insertIfEmpty(bucket, k, v)) {
        inserted = true;
        return;
      }
      std::swap(k, bucket.entries[offset].key);
      std::swap(v, bucket.entries[offset].value);
    };
    for (int r = 0; r < 15; ++r) {
      offset = (r + 1) % bucketSize;  // double check if it works
      PositionType idx1 = indexer.getHashIdx1(k);
      updateHelper(idx1, table1, swapUpdateFunc);
      if (inserted) {
        return false;
      }
      PositionType idx0 = indexer.getHashIdx0(k);
      updateHelper(idx0, table0, swapUpdateFunc);
      if (inserted) {
        return false;
      }
    }
    // printf("Warning CuckooHashMap uses stash\n");
    stash.push_back({true, k, v});
    if (stash.size() > 10) {
      throw std::runtime_error("CuckooHashMap insert failed");
    }
    return false;
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
                          const std::vector<CuckooHashMapEntry<K, V>>& stash) {
    if constexpr (isOblivious) {
      bool found = false;
      for (const auto& entry : stash) {
        bool match = entry.key == key;
        obliMove(match, value, entry.value);
        found |= match;
      }
      return found;
    } else {
      for (const auto& entry : stash) {
        if (entry.key == key) {
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
    }
    PositionType idx0 = indexer.getHashIdx0(key);
    if constexpr (isOblivious) {
      obliMove(isDummy, idx0, (PositionType)UniformRandom(tableSize - 1));
    }
    BucketType bucket;
    readHelper(idx0, table0, bucket);
    bool found = searchBucket(key, value, bucket);
    if constexpr (!isOblivious) {
      if (found) {
        return true;
      }
    }
    PositionType idx1 = indexer.getHashIdx1(key);
    if constexpr (isOblivious) {
      obliMove(isDummy, idx1, (PositionType)UniformRandom(tableSize - 1));
    }
    readHelper(idx1, table1, bucket);
    found |= searchBucket(key, value, bucket);
    if constexpr (!isOblivious) {
      if (found) {
        return true;
      }
    }
    found |= searchStash(key, value, stash);
    return found & !isDummy;
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
      if (it->key == key) {
        --load;
        stash.erase(it);
        return true;
      }
    }
    return false;
  }
};
}  // namespace ODSL