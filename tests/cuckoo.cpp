#include "oram/cuckoo.hpp"

#include <gtest/gtest.h>

#include "unordered_map"

using namespace ODSL;

template <bool isOblivious>
void testCuckooHashMap() {
  for (int r = 0; r < 1e2; ++r) {
    int mapSize = 5000;
    int keySpace = 10000;
    CuckooHashMap<int, int, isOblivious> map(mapSize);
    map.InitDefault();
    std::unordered_map<int, int> std_map;
    for (int r = 0; r < 2 * keySpace; ++r) {
      if (std_map.size() < mapSize) {
        int key = rand() % keySpace;
        int value = rand();
        map.insert(key, value);
        std_map[key] = value;
      }

      int key = rand() % keySpace;
      int value;
      bool foundFlag = map.find(key, value);
      auto it = std_map.find(key);
      if (it != std_map.end()) {
        ASSERT_TRUE(foundFlag);
        ASSERT_EQ(value, it->second);
      } else {
        ASSERT_FALSE(foundFlag);
      }
    }
  }
}

template <bool isOblivious>
void testCuckooHashMapInitFromReader() {
  for (int r = 0; r < 1e2; ++r) {
    int mapSize = 5000;
    int keySpace = 10000;
    CuckooHashMap<int, int, isOblivious> map(mapSize);
    std::unordered_map<int, int> std_map;
    for (int r = 0; r < mapSize; ++r) {
      int key = rand() % keySpace;
      int value = rand();
      std_map[key] = value;
    }
    auto it = std_map.begin();
    EM::VirtualVector::VirtualReader<std::pair<int, int>> reader(
        std_map.size(), [&it](uint64_t) { return *it++; });
    map.InitFromReaderInPlace(reader);
    for (int r = 0; r < 2 * keySpace; ++r) {
      if (std_map.size() < mapSize) {
        int key = rand() % keySpace;
        int value = rand();
        map.insert(key, value);
        std_map[key] = value;
      }

      int key = rand() % keySpace;
      int value;
      bool foundFlag = map.find(key, value);
      auto it = std_map.find(key);
      if (it != std_map.end()) {
        ASSERT_TRUE(foundFlag);
        ASSERT_EQ(value, it->second);
      } else {
        ASSERT_FALSE(foundFlag);
      }
    }
  }
}

template <bool isOblivious>
void testCuckooHashMapFindParBatch() {
  for (int r = 0; r < 10; ++r) {
    int mapSize = 23456;
    int keySpace = mapSize * 5;
    int batchSize = 1000;
    int numBatch = 10;
    int numThreads = 4;
    CuckooHashMap<int, int, isOblivious, uint64_t, true, true> map(
        mapSize, 1UL << 62, numThreads);
    std::unordered_map<int, int> std_map;
    for (int r = 0; r < mapSize; ++r) {
      int key = rand() % keySpace;
      int value = rand();
      std_map[key] = value;
    }
    auto it = std_map.begin();
    EM::VirtualVector::VirtualReader<std::pair<int, int>> reader(
        std_map.size(), [&it](uint64_t) { return *it++; });
    map.InitFromReaderInPlace(reader);
    for (int r = 0; r < numBatch; ++r) {
      std::vector<int> keys(batchSize);
      std::vector<int> vals(batchSize);
      for (int i = 0; i < batchSize; ++i) {
        keys[i] = rand() % keySpace;
      }
      if (std_map.size() < mapSize) {
        for (int i = 0; i < batchSize; ++i) {
          vals[i] = rand();
        }
        for (int i = 0; i < batchSize && std_map.size() < mapSize; ++i) {
          map.insert(keys[i], vals[i]);
          std_map[keys[i]] = vals[i];
        }
      }
      std::vector<uint8_t> findExistFlag = map.findParBatch(
          keys, vals, std::vector<bool>(batchSize, false), numThreads);
      for (size_t i = 0; i < keys.size(); ++i) {
        auto it = std_map.find(keys[i]);
        if (it != std_map.end()) {
          ASSERT_EQ(findExistFlag[i], true);
          ASSERT_EQ(vals[i], it->second);
        } else {
          ASSERT_EQ(findExistFlag[i], false);
        }
      }
    }
  }
}

template <bool isOblivious>
void testCuckooHashMapFindBatch() {
  for (int r = 0; r < 10; ++r) {
    int mapSize = 23456;
    int keySpace = mapSize;
    int batchSize = 1000;
    int numBatch = 10;
    using CHMap = CuckooHashMap<int, int, isOblivious, uint64_t, true>;
    CHMap map(mapSize,
                                                             1UL << 62);
    std::unordered_map<int, int> std_map;
    for (int r = 0; r < mapSize * 2 / 3; ++r) {
      int key = rand() % keySpace;
      int value = rand();
      std_map[key] = value;
    }
    auto it = std_map.begin();
    EM::VirtualVector::VirtualReader<std::pair<int, int>> reader(
        std_map.size(), [&it](uint64_t) { return *it++; });
    map.InitFromReaderInPlace(reader);
    for (int r = 0; r < numBatch; ++r) {
      std::vector<int> keys(batchSize);
      using ValResult = CHMap::ValResult;
      std::vector<ValResult> vals(batchSize);
      for (int i = 0; i < batchSize; ++i) {
        keys[i] = rand() % keySpace;
      }
      if (std_map.size() < mapSize) {
        // for (int i = 0; i < batchSize; ++i) {
        //   vals[i] = rand();
        // }
        for (int i = 0; i < batchSize && std_map.size() < mapSize; ++i) {
          map.find(keys[i], vals[i].value);
          // std_map[keys[i]] = vals[i];
        }
      }
      std::vector<uint8_t> findExistFlag;
      map.findBatchDeferWriteBack(
          keys.begin(), keys.end(), vals.begin());
#pragma omp parallel for num_threads(2)
      for (int i = 0; i < 2; ++i) {
        map.writeBackTable(i);
      }
      for (size_t i = 0; i < keys.size(); ++i) {
        auto it = std_map.find(keys[i]);
        if (it != std_map.end()) {
          ASSERT_EQ(vals[i].found, true);
          ASSERT_EQ(vals[i].value, it->second);
        } else {
          ASSERT_EQ(vals[i].found, false);
        }
      }
    }
  }
}

template <bool isOblivious>
void testCuckooHashMapFindParBatchPerf() {
  int mapSize = 100000;
  int keySpace = 1000;
  int batchSize = 1000;
  int numBatch = 1000;
  int numThreads = 8;
  CuckooHashMap<int64_t, int64_t, isOblivious, uint64_t, true, true> map(
      mapSize, 1UL << 62, numThreads);

  map.InitDefault();
  for (int r = 0; r < numBatch; ++r) {
    std::vector<int64_t> keys(batchSize);
    std::vector<int64_t> vals(batchSize);
    for (int i = 0; i < batchSize; ++i) {
      keys[i] = rand() % keySpace;
    }
    map.findParBatch(keys, vals, std::vector<bool>(batchSize, false),
                     numThreads);
  }
}

TEST(Cuckoo, CuckooHashMapNonOblivious) { testCuckooHashMap<false>(); }

TEST(Cuckoo, CuckooHashMapOblivious) { testCuckooHashMap<true>(); }

TEST(Cuckoo, CuckooHashMapInitFromReaderNonOblivious) {
  testCuckooHashMapInitFromReader<false>();
}

TEST(Cuckoo, CuckooHashMapInitFromReaderOblivious) {
  testCuckooHashMapInitFromReader<true>();
}

TEST(Cuckoo, CuckooHashMapFindBatchNonOblivious) {
  testCuckooHashMapFindParBatch<false>();
}

TEST(Cuckoo, CuckooHashMapFindParBatchOblivious) {
  testCuckooHashMapFindParBatch<true>();
}

TEST(Cuckoo, CuckooHashMapFindBatchOblivious) {
  testCuckooHashMapFindBatch<true>();
}

TEST(Cuckoo, CuckooHashMapFindParBatchPerfOblivious) {
  testCuckooHashMapFindParBatchPerf<true>();
}