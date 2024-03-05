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

TEST(Cuckoo, CuckooHashMapNonOblivious) { testCuckooHashMap<false>(); }

TEST(Cuckoo, CuckooHashMapOblivious) { testCuckooHashMap<true>(); }