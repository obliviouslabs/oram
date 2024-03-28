#include <gtest/gtest.h>

#include "odsl/omap.hpp"
#include "unordered_map"

using namespace ODSL;

template <bool isOblivious>
void testOHashMap() {
  for (int r = 0; r < 1e2; ++r) {
    int mapSize = 50;
    int keySpace = 100;
    OHashMap<int, int, isOblivious> map(mapSize);
    map.Init();
    std::unordered_map<int, int> std_map;
    for (int r = 0; r < 2 * keySpace; ++r) {
      if (std_map.size() < mapSize) {
        int key = rand() % keySpace;
        int value = rand();
        if constexpr (isOblivious) {
          map.OInsert(key, value);
        } else {
          map.Insert(key, value);
        }
        std_map[key] = value;
      }

      int key = rand() % keySpace;
      int value;
      bool foundFlag = map.Find(key, value);
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
  for (int rr = 0; rr < 1e2; ++rr) {
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
    map.InitFromReader(reader);
    for (int r = 0; r < 2 * keySpace; ++r) {
      if (std_map.size() < mapSize) {
        int key = rand() % keySpace;
        int value = rand();
        if constexpr (isOblivious) {
          map.OInsert(key, value);
        } else {
          map.Insert(key, value);
        }
        std_map[key] = value;
      }

      int key = rand() % keySpace;
      int value;
      bool foundFlag = map.Find(key, value);
      auto itFind = std_map.find(key);
      if (itFind != std_map.end()) {
        ASSERT_TRUE(foundFlag);
        ASSERT_EQ(value, itFind->second);
      } else {
        ASSERT_FALSE(foundFlag);
      }
    }
  }
}

void testOHashMapPushInit() {
  for (int rr = 0; rr < 1e2; ++rr) {
    int mapSize = 5000;
    int keySpace = 10000;
    OHashMap<int, int> map(mapSize);
    std::unordered_map<int, int> std_map;
    for (int r = 0; r < mapSize; ++r) {
      int key = rand() % keySpace;
      int value = rand();
      std_map[key] = value;
    }
    auto* initContext = map.NewInitContext(1UL << 25);
    for (auto it = std_map.begin(); it != std_map.end();) {
      if (rand() % 2) {
        int batchSize = (rand() % 100) + 1;
        std::vector<std::pair<int, int>> vec;
        for (int i = 0; it != std_map.end() && i < batchSize; ++it, ++i) {
          vec.push_back(*it);
        }
        initContext->InsertBatch(vec.begin(), vec.end());
      } else {
        initContext->Insert(it->first, it->second);
        ++it;
      }
    }
    initContext->Finalize();
    delete initContext;
    initContext = nullptr;
    ASSERT_EQ(map.GetLoad(), std_map.size());
    for (int r = 0; r < 2 * keySpace; ++r) {
      if (std_map.size() < mapSize) {
        int key = rand() % keySpace;
        int value = rand();

        map.OInsert(key, value);
        std_map[key] = value;
      }

      int key = rand() % keySpace;
      int value;
      bool foundFlag = map.Find(key, value);
      auto it = std_map.find(key);
      if (it != std_map.end()) {
        ASSERT_TRUE(foundFlag);
        ASSERT_EQ(value, it->second);
      } else {
        ASSERT_FALSE(foundFlag);
      }
    }
    ASSERT_EQ(map.GetLoad(), std_map.size());
  }
}

void testOHashMapInitFromNonObliviousWithDummy() {
  for (int rr = 0; rr < 1e1; ++rr) {
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
        nonOMap.Insert<true>(it->first, it->second, true);
      }
      nonOMap.Insert<true>(it->first, it->second);
    }
    map.InitFromNonOblivious(nonOMap);
    ASSERT_EQ(map.GetLoad(), std_map.size());
    for (int r = 0; r < 2 * keySpace; ++r) {
      if (std_map.size() < mapSize) {
        int key = rand() % keySpace;
        int value = rand();
        map.OInsert(key, value);
        std_map[key] = value;
      }

      int key = rand() % keySpace;
      int value;
      bool foundFlag = map.Find(key, value);
      auto it = std_map.find(key);
      if (it != std_map.end()) {
        ASSERT_TRUE(foundFlag);
        ASSERT_EQ(value, it->second);
      } else {
        ASSERT_FALSE(foundFlag);
      }
    }
    ASSERT_EQ(map.GetLoad(), std_map.size());
  }
}

template <bool isOblivious>
void testOHashMapFindBatch() {
  for (int rr = 0; rr < 10; ++rr) {
    int mapSize = 23456;
    int keySpace = mapSize;
    int batchSize = 1000;
    int numBatch = 10;
    using CHMap = OHashMap<int, int, isOblivious, uint64_t, true>;
    CHMap map(mapSize, MAX_CACHE_SIZE);
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
    ASSERT_EQ(map.GetLoad(), std_map.size());
    for (int r = 0; r < numBatch; ++r) {
      std::vector<int> keys(batchSize);
      using ValResult = CHMap::ValResult;
      std::vector<ValResult> vals(batchSize);
      for (int i = 0; i < batchSize; ++i) {
        keys[i] = rand() % keySpace;
      }
      std::vector<uint8_t> findExistFlag;
      map.FindBatchDeferWriteBack(keys.begin(), keys.end(), vals.begin());
#pragma omp parallel for num_threads(2)
      for (int i = 0; i < 2; ++i) {
        map.WriteBackTable(i);
      }
      for (size_t i = 0; i < keys.size(); ++i) {
        auto itFind = std_map.find(keys[i]);
        if (itFind != std_map.end()) {
          ASSERT_EQ(vals[i].found, true);
          ASSERT_EQ(vals[i].value, itFind->second);
        } else {
          ASSERT_EQ(vals[i].found, false);
        }
      }
    }
  }
}

template <bool isOblivious>
void testReplaceCount() {
  // test replace count distribution
  int mapSize = 100000;
  int round = 10000;
  int outerRound = 500;

  int windowSize = 5;

  std::vector<uint64_t> stashLoads(30, 0);
  for (int rr = 0; rr < outerRound; ++rr) {
    OHashMap<int, int, false> map(mapSize, MAX_CACHE_SIZE);
    map.Init();
    const auto& stash = map.GetStash();
    for (int i = 0; i < mapSize; ++i) {
      map.Insert(rand(), 0);
    }
    for (int r = 0; r < round; ++r) {
      int key = rand();
      if constexpr (isOblivious) {
        map.OInsert(key, 0);
        map.OErase(key);
      } else {
        map.Insert(key, 0);
        map.Erase(key);
      }
      if (r % windowSize == 0) {
        int load = 0;
        for (int k = 0; k < stash.size(); ++k) {
          if (stash[k].valid) {
            ++load;
          }
        }
        stashLoads[load]++;
      }
    }
  }
  for (int i = 0; i < stashLoads.size(); ++i) {
    printf("%d %lu\n", i, stashLoads[i]);
  }
  for (int i = 10; i < stashLoads.size(); ++i) {
    // stash load should be less than 10 with high probability
    ASSERT_EQ(stashLoads[i], 0);
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

TEST(Cuckoo, ReplaceCountDistriNonOblivious) { testReplaceCount<false>(); }

TEST(Cuckoo, ReplaceCountDistriOblivious) { testReplaceCount<true>(); }

TEST(Cuckoo, OHashMapPushInit) { testOHashMapPushInit(); }