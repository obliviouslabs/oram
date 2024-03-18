#include <gtest/gtest.h>

#include "oram/omap.hpp"
#include "unordered_map"

using namespace ODSL;

template <bool isOblivious>
void testOHashMap() {
  for (int r = 0; r < 1e2; ++r) {
    int mapSize = 50;
    int keySpace = 100;
    OHashMap<int, int, isOblivious> map(mapSize);
    map.InitDefault();
    std::unordered_map<int, int> std_map;
    for (int r = 0; r < 2 * keySpace; ++r) {
      if (std_map.size() < mapSize) {
        int key = rand() % keySpace;
        int value = rand();
        if constexpr (isOblivious) {
          map.insertOblivious(key, value);
        } else {
          map.insert(key, value);
        }
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
void testOHashMapInitFromReader() {
  for (int r = 0; r < 1e2; ++r) {
    int mapSize = 5000;
    int keySpace = 10000;
    OHashMap<int, int, isOblivious> map(mapSize);
    std::unordered_map<int, int> std_map;
    for (int r = 0; r < mapSize; ++r) {
      int key = rand() % keySpace;
      int value = rand();
      std_map[key] = value;
    }
    auto it = std_map.begin();
    EM::VirtualVector::VirtualReader<std::pair<int, int>> reader(
        std_map.size(), [&it](uint64_t) { return *it++; });
    map.InitFromReader(reader, 1UL << 26);
    for (int r = 0; r < 2 * keySpace; ++r) {
      if (std_map.size() < mapSize) {
        int key = rand() % keySpace;
        int value = rand();
        if constexpr (isOblivious) {
          map.insertOblivious(key, value);
        } else {
          map.insert(key, value);
        }
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

void testOHashMapInitFromNonObliviousWithDummy() {
  for (int r = 0; r < 1e1; ++r) {
    int mapSize = 5000;
    int keySpace = 10000;
    OHashMap<int, int, true> map(mapSize);
    OHashMap<int, int, false> nonOMap(mapSize);
    std::unordered_map<int, int> std_map;
    for (int r = 0; r < mapSize; ++r) {
      int key = rand() % keySpace;
      int value = rand();
      std_map[key] = value;
    }
    for (auto it = std_map.begin(); it != std_map.end(); ++it) {
      if (rand() % 5 == 0 && nonOMap.GetLoad() < mapSize) {
        nonOMap.insert<true>(it->first, it->second, true);
      }
      nonOMap.insert<true>(it->first, it->second);
    }
    map.InitFromNonOblivious(nonOMap);
    // printf("init done\n");
    for (int r = 0; r < 2 * keySpace; ++r) {
      if (std_map.size() < mapSize) {
        int key = rand() % keySpace;
        int value = rand();
        map.insertOblivious(key, value);
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
void testOHashMapFindBatch() {
  for (int r = 0; r < 10; ++r) {
    int mapSize = 23456;
    int keySpace = mapSize;
    int batchSize = 1000;
    int numBatch = 10;
    using CHMap = OHashMap<int, int, isOblivious, uint64_t, true>;
    CHMap map(mapSize, 1UL << 62);
    std::unordered_map<int, int> std_map;
    for (int r = 0; r < mapSize * 2 / 3; ++r) {
      int key = rand() % keySpace;
      int value = rand();
      std_map[key] = value;
    }
    auto it = std_map.begin();
    EM::VirtualVector::VirtualReader<std::pair<int, int>> reader(
        std_map.size(), [&it](uint64_t) { return *it++; });
    map.InitFromReader(reader);
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
      map.findBatchDeferWriteBack(keys.begin(), keys.end(), vals.begin());
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

TEST(Cuckoo, OHashMapNonOblivious) { testOHashMap<false>(); }

TEST(Cuckoo, OHashMapOblivious) { testOHashMap<true>(); }

TEST(Cuckoo, OHashMapInitFromReaderNonOblivious) {
  testOHashMapInitFromReader<false>();
}

TEST(Cuckoo, OHashMapInitFromReaderOblivious) {
  testOHashMapInitFromReader<true>();
}

TEST(Cuckoo, OHashMapFindBatchOblivious) { testOHashMapFindBatch<true>(); }

TEST(Cuckoo, OHashMapInitWithDummy) {
  testOHashMapInitFromNonObliviousWithDummy();
}