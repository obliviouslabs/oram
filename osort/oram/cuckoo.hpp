#pragma once
#include "recursive_oram.hpp"
namespace ODSL {

template <typename K, typename PositionType = uint64_t>
struct CuckooHashMapIndexer {
  static constexpr int saltLength = 16;
  uint8_t salts[2][saltLength];
  PositionType _size;
  CuckooHashMapIndexer(PositionType size) : _size(size) {
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
};

template <typename K, typename V, const short bucketSize>
struct CuckooHashMapBucket {
  CuckooHashMapEntry<K, V> entries[bucketSize];
};

template <typename K, typename V, typename PositionType = uint64_t>
struct NonObliviousCuckooHashMap {
  static constexpr int saltLength = 16;
  static constexpr double loadFactor = 0.7;
  static constexpr short bucketSize = 2;
  PositionType _size;
  PositionType tableSize;
  StdVector<CuckooHashMapBucket<K, V, bucketSize>> table0, table1;
  CuckooHashMapIndexer<K, PositionType> indexer;
  NonObliviousCuckooHashMap(PositionType size)
      : _size(size),
        indexer(size / (2 * loadFactor * bucketSize)),
        tableSize(size / (2 * loadFactor * bucketSize)),
        table0(size / (2 * loadFactor * bucketSize)),
        table1(size / (2 * loadFactor * bucketSize)) {}

  void insert(const K& key, const V& value) {
    PositionType idx0 = indexer.getHashIdx0(key);
    for (int i = 0; i < bucketSize; ++i) {
      if (!(table0[idx0].entries[i].valid &&
            table0[idx0].entries[i].key != key)) {
        table0[idx0].entries[i] = {true, key, value};
        return;
      }
    }
    PositionType idx1 = indexer.getHashIdx1(key);
    for (int i = 0; i < bucketSize; ++i) {
      if (!(table1[idx1].entries[i].valid &&
            table1[idx1].entries[i].key != key)) {
        table1[idx1].entries[i] = {true, key, value};
        return;
      }
    }
    K k = key;
    V v = value;
    for (int r = 0; r < 100; ++r) {
      int offset = r % bucketSize;
      std::swap(k, table0[idx0].entries[offset].key);
      std::swap(v, table0[idx0].entries[offset].value);
      PositionType idx1 = indexer.getHashIdx1(k);
      for (int i = 0; i < bucketSize; ++i) {
        if (!table1[idx1].entries[i].valid) {
          table1[idx1].entries[i] = {true, k, v};
          return;
        }
      }
      std::swap(k, table1[idx1].entries[offset].key);
      std::swap(v, table1[idx1].entries[offset].value);
      idx0 = indexer.getHashIdx0(k);
      for (int i = 0; i < bucketSize; ++i) {
        if (!table0[idx0].entries[i].valid) {
          table0[idx0].entries[i] = {true, k, v};
          return;
        }
      }
    }

    throw std::runtime_error("CuckooHashMap insert failed");
  }

  bool find(const K& key, V& value) const {
    PositionType idx0 = indexer.getHashIdx0(key);
    for (int i = 0; i < bucketSize; ++i) {
      if (table0[idx0].entries[i].valid && table0[idx0].entries[i].key == key) {
        value = table0[idx0].entries[i].value;
        return true;
      }
    }
    PositionType idx1 = indexer.getHashIdx1(key);
    for (int i = 0; i < bucketSize; ++i) {
      if (table1[idx1].entries[i].valid && table1[idx1].entries[i].key == key) {
        value = table1[idx1].entries[i].value;
        return true;
      }
    }
    return false;
  }
};
}  // namespace ODSL