#include "oram/par_omap.hpp"

#include "testutils.hpp"

using namespace ODSL;
TEST(ParOMap, maxQueryPerShard) {
  uint64_t batchSize = 1000;
  uint64_t shardCount = 8;
  double logFailProb = -40;
  uint64_t maxQuery = ParOMap<uint64_t, uint64_t>::maxQueryPerShard(
      batchSize, shardCount, logFailProb);
  std::cout << "maxQuery: " << maxQuery << std::endl;
}

TEST(ParOMap, Init) {
  uint64_t mapSize = 1234;
  uint64_t shardCount = 8;
  ParOMap<uint64_t, uint64_t> parOMap(mapSize, shardCount);
  parOMap.Init();
}

TEST(ParOMap, InsertAndFindSimple) {
  uint64_t mapSize = 1234;
  uint64_t shardCount = 2;
  ParOMap<uint64_t, uint64_t> parOMap(mapSize, shardCount);
  parOMap.Init();
  std::vector<uint64_t> keys = {1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<uint64_t> vals = {10, 20, 30, 40, 50, 60, 70, 80};
  std::vector<bool> insertExistFlag =
      parOMap.insertBatch(keys.begin(), keys.end(), vals.begin());
  for (size_t i = 0; i < keys.size(); ++i) {
    ASSERT_FALSE(insertExistFlag[i]);
  }
  std::vector<uint64_t> foundVals(keys.size());
  std::vector<bool> findExistFlag =
      parOMap.findBatch(keys.begin(), keys.end(), foundVals.begin());
  for (size_t i = 0; i < keys.size(); ++i) {
    ASSERT_TRUE(findExistFlag[i]);
    ASSERT_EQ(foundVals[i], vals[i]);
  }
}

TEST(ParOMap, InsertAndFind) {
  uint64_t mapSize = 1234;
  uint64_t shardCount = 2;
  uint64_t round = 1000;
  uint64_t batchSize = 10;
  ParOMap<uint64_t, uint64_t> parOMap(mapSize, shardCount);
  parOMap.Init();
  std::unordered_map<uint64_t, uint64_t> kvMap;
  for (uint64_t r = 0; r < round; ++r) {
    std::vector<uint64_t> keys(batchSize);
    std::vector<uint64_t> vals(batchSize);
    for (uint64_t i = 0; i < batchSize; ++i) {
      keys[i] = UniformRandom() % mapSize;
      vals[i] = UniformRandom();
    }
    std::vector<bool> insertExistFlag =
        parOMap.insertBatch(keys.begin(), keys.end(), vals.begin());
    for (size_t i = 0; i < keys.size(); ++i) {
      if (insertExistFlag[i] != (kvMap.count(keys[i]) > 0)) {
        std::cout << "key = " << keys[i] << std::endl;
      }
      ASSERT_EQ(insertExistFlag[i], kvMap.count(keys[i]) > 0);
      kvMap[keys[i]] = vals[i];
    }
    for (uint64_t i = 0; i < batchSize; ++i) {
      keys[i] = UniformRandom() % mapSize;
    }
    std::vector<uint64_t> foundVals(keys.size());
    std::vector<bool> findExistFlag =
        parOMap.findBatch(keys.begin(), keys.end(), foundVals.begin());
    for (size_t i = 0; i < keys.size(); ++i) {
      if (findExistFlag[i] != (kvMap.count(keys[i]) > 0)) {
        std::cout << "key = " << keys[i] << std::endl;
      }
      ASSERT_EQ(findExistFlag[i], kvMap.count(keys[i]) > 0);
      if (findExistFlag[i]) {
        ASSERT_EQ(foundVals[i], kvMap[keys[i]]);
      }
    }
  }
}