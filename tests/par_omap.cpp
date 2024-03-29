#include "odsl/par_omap.hpp"

#include <unordered_set>

#include "testutils.hpp"

using namespace ODSL;

TEST(ParOMap, InitInsertFind) {
  uint64_t mapSize = UniformRandom(100, 5000000);
  uint64_t shardCount = 27;
  uint64_t round = 1000;
  uint64_t batchSize = 100;
  ParOMap<uint64_t, uint64_t> parOMap(mapSize, shardCount);
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
  size_t BackendSize = 2e9;
  EM::Backend::g_DefaultBackend =
      new EM::Backend::MemServerBackend(BackendSize);
  std::map<uint64_t, uint64_t> kvMap;
  for (uint64_t i = 0; i < mapSize; ++i) {
    if (UniformRandomBit()) {
      kvMap[UniformRandom()] = UniformRandom();
    }
  }
  auto kvIt = kvMap.begin();
  EM::VirtualVector::VirtualReader<std::pair<uint64_t, uint64_t>> reader(
      kvMap.size(), [&](uint64_t i) {
        const auto& pr = *kvIt;
        ++kvIt;
        return pr;
      });
  parOMap.InitFromReader(reader, 128UL << 20);
  std::cout << "omp max threads: " << omp_get_max_threads() << std::endl;
  for (uint64_t r = 0; r < round; ++r) {
    std::vector<uint64_t> keys(batchSize);
    std::vector<uint64_t> vals(batchSize);
    if (r % 10 == 0) {  // perform fewer inserts than finds
      std::unordered_set<uint64_t> keySet;
      for (uint64_t i = 0; i < batchSize; ++i) {
        do {
          keys[i] = UniformRandom() % mapSize;
        } while (keySet.count(keys[i]) > 0);
        keySet.insert(keys[i]);
        vals[i] = UniformRandom();
      }
      std::vector<uint8_t> insertExistFlag =
          parOMap.InsertBatch(keys.begin(), keys.end(), vals.begin());
      for (size_t i = 0; i < keys.size(); ++i) {
        if (insertExistFlag[i] != (kvMap.count(keys[i]) > 0)) {
          std::cout << "key = " << keys[i] << std::endl;
        }
        ASSERT_EQ(insertExistFlag[i], kvMap.count(keys[i]) > 0);
        kvMap[keys[i]] = vals[i];
      }
    }
    for (uint64_t i = 0; i < batchSize; ++i) {
      keys[i] = UniformRandom() % mapSize;
    }
    std::vector<uint64_t> foundVals(keys.size());
    std::vector<uint8_t> findExistFlag =
        parOMap.FindBatch(keys.begin(), keys.end(), foundVals.begin());
    for (size_t i = 0; i < keys.size(); ++i) {
      if (findExistFlag[i] != (kvMap.count(keys[i]) > 0)) {
        std::cout << "key = " << keys[i] << std::endl;
        std::cout << "expected value " << kvMap[keys[i]] << std::endl;
        std::cout << "found value " << foundVals[i] << std::endl;
      }
      ASSERT_EQ(findExistFlag[i], kvMap.count(keys[i]) > 0);
      if (findExistFlag[i]) {
        ASSERT_EQ(foundVals[i], kvMap[keys[i]]);
      }
    }
  }
}

TEST(ParOMap, InitCatchDup) {
  uint64_t mapSize = UniformRandom(100, 500000);
  uint64_t shardCount = 5;
  ParOMap<uint64_t, uint64_t> parOMap(mapSize, shardCount);
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
  size_t BackendSize = 2e9;
  EM::Backend::g_DefaultBackend =
      new EM::Backend::MemServerBackend(BackendSize);

  size_t initSize = mapSize * 2 / 3;
  EM::VirtualVector::VirtualReader<std::pair<uint64_t, uint64_t>> reader(
      initSize,
      [&](uint64_t i) { return std::make_pair(i % (initSize / 2), 0UL); });
  try {
    parOMap.InitFromReader(reader, 128UL << 20);
    ASSERT_TRUE(false);  // should not reach here
  } catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}

TEST(ParOMap, PushInitInsertFind) {
  uint64_t mapSize = UniformRandom(100, 500000);
  uint64_t shardCount = 32;
  uint64_t round = 1000;
  ParOMap<uint64_t, uint64_t> parOMap(mapSize, shardCount);
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
  size_t BackendSize = 2e9;
  EM::Backend::g_DefaultBackend =
      new EM::Backend::MemServerBackend(BackendSize);
  std::map<uint64_t, uint64_t> kvMap;
  for (uint64_t i = 0; i < mapSize; ++i) {
    kvMap[UniformRandom()] = UniformRandom();
  }
  auto initContext = parOMap.NewInitContext(kvMap.size(), 1UL << 28);
  for (auto it = kvMap.begin(); it != kvMap.end();) {
    if (rand() % 2) {
      int batchSize = (rand() % 100) + 1;
      std::vector<std::pair<uint64_t, uint64_t>> vec;
      for (int i = 0; it != kvMap.end() && i < batchSize; ++it, ++i) {
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
  std::cout << "omp max threads: " << omp_get_max_threads() << std::endl;
  for (uint64_t r = 0; r < round; ++r) {
    int batchSize = (rand() % 100) + 1;
    std::vector<uint64_t> keys(batchSize);
    std::vector<uint64_t> vals(batchSize);
    if (r % 10 == 0) {  // perform fewer inserts than finds
      std::unordered_set<uint64_t> keySet;
      for (uint64_t i = 0; i < batchSize; ++i) {
        do {
          keys[i] = UniformRandom() % mapSize;
        } while (keySet.count(keys[i]) > 0);
        keySet.insert(keys[i]);
        vals[i] = UniformRandom();
      }
      std::vector<uint8_t> insertExistFlag =
          parOMap.InsertBatch(keys.begin(), keys.end(), vals.begin());
      for (size_t i = 0; i < keys.size(); ++i) {
        if (insertExistFlag[i] != (kvMap.count(keys[i]) > 0)) {
          std::cout << "key = " << keys[i] << std::endl;
        }
        ASSERT_EQ(insertExistFlag[i], kvMap.count(keys[i]) > 0);
        kvMap[keys[i]] = vals[i];
      }
    }
    for (uint64_t i = 0; i < batchSize; ++i) {
      keys[i] = UniformRandom() % mapSize;
    }
    std::vector<uint64_t> foundVals(keys.size());
    std::vector<uint8_t> findExistFlag =
        parOMap.FindBatch(keys.begin(), keys.end(), foundVals.begin());
    for (size_t i = 0; i < keys.size(); ++i) {
      if (findExistFlag[i] != (kvMap.count(keys[i]) > 0)) {
        std::cout << "key = " << keys[i] << std::endl;
        std::cout << "expected value " << kvMap[keys[i]] << std::endl;
        std::cout << "found value " << foundVals[i] << std::endl;
      }
      ASSERT_EQ(findExistFlag[i], kvMap.count(keys[i]) > 0);
      if (findExistFlag[i]) {
        ASSERT_EQ(foundVals[i], kvMap[keys[i]]);
      }
    }
  }
}

TEST(ParOMap, InsertAndFindSimple) {
  uint64_t mapSize = 1234;
  uint64_t shardCount = 2;
  ParOMap<uint64_t, uint64_t> parOMap(mapSize, shardCount);
  parOMap.Init();
  std::vector<uint64_t> keys = {1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<uint64_t> vals = {10, 20, 30, 40, 50, 60, 70, 80};
  std::vector<uint8_t> insertExistFlag =
      parOMap.InsertBatch(keys.begin(), keys.end(), vals.begin());
  for (size_t i = 0; i < keys.size(); ++i) {
    ASSERT_FALSE(insertExistFlag[i]);
  }
  std::vector<uint64_t> foundVals(keys.size());
  std::vector<uint8_t> findExistFlag =
      parOMap.FindBatch(keys.begin(), keys.end(), foundVals.begin());
  for (size_t i = 0; i < keys.size(); ++i) {
    ASSERT_TRUE(findExistFlag[i]);
    ASSERT_EQ(foundVals[i], vals[i]);
  }
}

TEST(ParOMap, InsertAndFind) {
  uint64_t mapSize = 123456;
  uint64_t shardCount = 8;
  uint64_t round = 1000;
  uint64_t batchSize = 1000;
  ParOMap<uint64_t, uint64_t> parOMap(mapSize, shardCount);
  parOMap.Init();
  std::unordered_map<uint64_t, uint64_t> kvMap;
  std::cout << "omp max threads: " << omp_get_max_threads() << std::endl;
  for (uint64_t r = 0; r < round; ++r) {
    std::vector<uint64_t> keys(batchSize);
    std::vector<uint64_t> vals(batchSize);
    if (r % 10 == 0) {  // perform fewer inserts than finds
      std::unordered_set<uint64_t> keySet;
      for (uint64_t i = 0; i < batchSize; ++i) {
        do {
          keys[i] = UniformRandom() % mapSize;
        } while (keySet.count(keys[i]) > 0);
        keySet.insert(keys[i]);
        vals[i] = UniformRandom();
      }
      std::vector<uint8_t> insertExistFlag =
          parOMap.InsertBatch(keys.begin(), keys.end(), vals.begin());
      for (size_t i = 0; i < keys.size(); ++i) {
        if (insertExistFlag[i] != (kvMap.count(keys[i]) > 0)) {
          std::cout << "key = " << keys[i] << std::endl;
        }
        ASSERT_EQ(insertExistFlag[i], kvMap.count(keys[i]) > 0);
        kvMap[keys[i]] = vals[i];
      }
    }
    for (uint64_t i = 0; i < batchSize; ++i) {
      keys[i] = UniformRandom() % mapSize;
    }
    std::vector<uint64_t> foundVals(keys.size());
    std::vector<uint8_t> findExistFlag =
        parOMap.FindBatch(keys.begin(), keys.end(), foundVals.begin());
    for (size_t i = 0; i < keys.size(); ++i) {
      if (findExistFlag[i] != (kvMap.count(keys[i]) > 0)) {
        std::cout << "key = " << keys[i] << std::endl;
        std::cout << "expected value " << kvMap[keys[i]] << std::endl;
        std::cout << "found value " << foundVals[i] << std::endl;
      }
      ASSERT_EQ(findExistFlag[i], kvMap.count(keys[i]) > 0);
      if (findExistFlag[i]) {
        ASSERT_EQ(foundVals[i], kvMap[keys[i]]);
      }
    }
  }
}

TEST(ParOMap, InsertAndDuplicateFind) {
  uint64_t mapSize = 123456;
  uint64_t shardCount = 8;
  uint64_t round = 1000;
  uint64_t batchSize = 1000;
  ParOMap<uint64_t, uint64_t> parOMap(mapSize, shardCount);
  parOMap.Init();
  std::unordered_map<uint64_t, uint64_t> kvMap;
  std::cout << "omp max threads: " << omp_get_max_threads() << std::endl;
  for (uint64_t r = 0; r < round; ++r) {
    std::vector<uint64_t> keys(batchSize);
    std::vector<uint64_t> vals(batchSize);
    if (r % 10 == 0) {  // perform fewer inserts than finds
      std::unordered_set<uint64_t> keySet;
      for (uint64_t i = 0; i < batchSize; ++i) {
        do {
          keys[i] = UniformRandom() % mapSize;
        } while (keySet.count(keys[i]) > 0);
        keySet.insert(keys[i]);
        vals[i] = UniformRandom();
      }
      std::vector<uint8_t> insertExistFlag =
          parOMap.InsertBatch(keys.begin(), keys.end(), vals.begin());
      for (size_t i = 0; i < keys.size(); ++i) {
        if (insertExistFlag[i] != (kvMap.count(keys[i]) > 0)) {
          std::cout << "key = " << keys[i] << std::endl;
        }
        ASSERT_EQ(insertExistFlag[i], kvMap.count(keys[i]) > 0);
        kvMap[keys[i]] = vals[i];
      }
    }
    uint64_t base = UniformRandom() % mapSize;
    for (uint64_t i = 0; i < batchSize; ++i) {
      keys[i] = base + UniformRandom() % 10;
    }
    std::vector<uint64_t> foundVals(keys.size());
    std::vector<uint8_t> findExistFlag =
        parOMap.FindBatch(keys.begin(), keys.end(), foundVals.begin());
    for (size_t i = 0; i < keys.size(); ++i) {
      if (findExistFlag[i] != (kvMap.count(keys[i]) > 0)) {
        std::cout << "key = " << keys[i] << std::endl;
        std::cout << "expected value " << kvMap[keys[i]] << std::endl;
        std::cout << "found value " << foundVals[i] << std::endl;
      }
      ASSERT_EQ(findExistFlag[i], kvMap.count(keys[i]) > 0);
      if (findExistFlag[i]) {
        ASSERT_EQ(foundVals[i], kvMap[keys[i]]);
      }
    }
  }
}

TEST(ParOMap, InsertFindErase) {
  uint64_t mapSize = 12345;
  uint64_t shardCount = 8;
  uint64_t round = 1000;
  uint64_t batchSize = 1000;
  ParOMap<uint64_t, uint64_t> parOMap(mapSize, shardCount);
  parOMap.Init();
  std::unordered_map<uint64_t, uint64_t> kvMap;
  std::cout << "omp max threads: " << omp_get_max_threads() << std::endl;
  for (uint64_t r = 0; r < round; ++r) {
    std::vector<uint64_t> keys(batchSize);
    std::vector<uint64_t> vals(batchSize);
    if (UniformRandom32() % 10 == 0) {  // perform fewer inserts than finds
      std::unordered_set<uint64_t> keySet;
      for (uint64_t i = 0; i < batchSize; ++i) {
        do {
          keys[i] = UniformRandom() % mapSize;
        } while (keySet.count(keys[i]) > 0);
        keySet.insert(keys[i]);
        vals[i] = UniformRandom();
      }
      std::vector<uint8_t> insertExistFlag =
          parOMap.InsertBatch(keys.begin(), keys.end(), vals.begin());
      for (size_t i = 0; i < keys.size(); ++i) {
        if (insertExistFlag[i] != (kvMap.count(keys[i]) > 0)) {
          std::cout << "key = " << keys[i] << std::endl;
        }
        ASSERT_EQ(insertExistFlag[i], kvMap.count(keys[i]) > 0);
        kvMap[keys[i]] = vals[i];
      }
    }
    if (UniformRandom32() % 10 == 0) {
      for (uint64_t i = 0; i < batchSize; ++i) {
        keys[i] = UniformRandom() % mapSize;
      }
      std::vector<uint8_t> eraseExistFlag =
          parOMap.EraseBatch(keys.begin(), keys.end());
      for (size_t i = 0; i < keys.size(); ++i) {
        if (eraseExistFlag[i] != (kvMap.count(keys[i]) > 0)) {
          std::cout << "key = " << keys[i] << std::endl;
        }
        ASSERT_EQ(eraseExistFlag[i], kvMap.count(keys[i]) > 0);
      }
      for (size_t i = 0; i < keys.size(); ++i) {
        kvMap.erase(keys[i]);
      }
    }
    for (uint64_t i = 0; i < batchSize; ++i) {
      keys[i] = UniformRandom() % mapSize;
    }
    std::vector<uint64_t> foundVals(keys.size());
    std::vector<uint8_t> findExistFlag =
        parOMap.FindBatch(keys.begin(), keys.end(), foundVals.begin());
    for (size_t i = 0; i < keys.size(); ++i) {
      if (findExistFlag[i] != (kvMap.count(keys[i]) > 0)) {
        std::cout << "key = " << keys[i] << std::endl;
        std::cout << "expected value " << kvMap[keys[i]] << std::endl;
        std::cout << "found value " << foundVals[i] << std::endl;
      }
      ASSERT_EQ(findExistFlag[i], kvMap.count(keys[i]) > 0);
      if (findExistFlag[i]) {
        ASSERT_EQ(foundVals[i], kvMap[keys[i]]);
      }
    }
  }
}

TEST(ParOMap, MixedLarge) {
  int round = 10;
  size_t BackendSize = 4e9;
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
  EM::Backend::g_DefaultBackend =
      new EM::Backend::MemServerBackend(BackendSize);
  for (int rr = 0; rr < round; ++rr) {
    std::cout << "Test round " << rr << std::endl;
    uint64_t mapSize = UniformRandom(100, 5000000);
    uint64_t shardCount = 1UL << UniformRandom(1, 5);
    uint64_t accessRound = 100;

    ParOMap<uint64_t, uint64_t> parOMap(mapSize, shardCount);

    std::map<uint64_t, uint64_t> kvMap;
    size_t mapInitSize = UniformRandom(mapSize);
    for (uint64_t i = 0; i < mapInitSize; ++i) {
      if (UniformRandomBit()) {
        kvMap[i] = UniformRandom();
      }
    }
    auto kvIt = kvMap.begin();
    EM::VirtualVector::VirtualReader<std::pair<uint64_t, uint64_t>> reader(
        kvMap.size(), [&](uint64_t i) {
          const auto& pr = *kvIt;
          ++kvIt;
          return pr;
        });
    size_t cacheSize = UniformRandom(1UL << 26, 1UL << 30);
    parOMap.InitFromReader(reader, cacheSize);
    for (uint64_t r = 0; r < accessRound; ++r) {
      uint64_t batchSize = UniformRandom(1, 1000);
      std::vector<uint64_t> keys(batchSize);
      std::vector<uint64_t> vals(batchSize);
      // perform fewer inserts than finds
      if (rand() % 10 == 0 && kvMap.size() < mapSize) {
        std::unordered_set<uint64_t> keySet;
        for (uint64_t i = 0; i < batchSize; ++i) {
          do {
            keys[i] = UniformRandom() % mapSize;
          } while (keySet.count(keys[i]) > 0);
          keySet.insert(keys[i]);
          vals[i] = UniformRandom();
        }
        std::vector<uint8_t> insertExistFlag =
            parOMap.InsertBatch(keys.begin(), keys.end(), vals.begin());
        for (size_t i = 0; i < keys.size(); ++i) {
          if (insertExistFlag[i] != (kvMap.count(keys[i]) > 0)) {
            std::cout << "key = " << keys[i] << std::endl;
          }
          ASSERT_EQ(insertExistFlag[i], kvMap.count(keys[i]) > 0);
          kvMap[keys[i]] = vals[i];
        }
      }
      if (rand() % 10 == 0) {
        for (uint64_t i = 0; i < batchSize; ++i) {
          keys[i] = UniformRandom() % mapSize;
        }
        std::vector<uint8_t> eraseExistFlag =
            parOMap.EraseBatch(keys.begin(), keys.end());
        for (size_t i = 0; i < keys.size(); ++i) {
          if (eraseExistFlag[i] != (kvMap.count(keys[i]) > 0)) {
            std::cout << "key = " << keys[i] << std::endl;
          }
          ASSERT_EQ(eraseExistFlag[i], kvMap.count(keys[i]) > 0);
        }
        for (size_t i = 0; i < keys.size(); ++i) {
          kvMap.erase(keys[i]);
        }
      }
      for (uint64_t i = 0; i < batchSize; ++i) {
        keys[i] = UniformRandom() % mapSize;
      }
      std::vector<uint64_t> foundVals(keys.size());
      std::vector<uint8_t> findExistFlag =
          parOMap.FindBatch(keys.begin(), keys.end(), foundVals.begin());
      for (size_t i = 0; i < keys.size(); ++i) {
        if (findExistFlag[i] != (kvMap.count(keys[i]) > 0)) {
          std::cout << "key = " << keys[i] << std::endl;
          std::cout << "expected value " << kvMap[keys[i]] << std::endl;
          std::cout << "found value " << foundVals[i] << std::endl;
        }
        ASSERT_EQ(findExistFlag[i], kvMap.count(keys[i]) > 0);
        if (findExistFlag[i]) {
          ASSERT_EQ(foundVals[i], kvMap[keys[i]]);
        }
      }
    }
  }
}