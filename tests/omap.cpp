#include "oram/omap.hpp"

#include <gtest/gtest.h>

#include <unordered_map>

#include "external_memory/algorithm/sort_def.hpp"
#include "testutils.hpp"

using namespace ORAM;

TEST(OMap, Init) {
  OMap<uint64_t, SortElement> omap(512);
  StdVector<std::pair<uint64_t, SortElement>> vec(16);
  for (int i = 0; i < 16; i++) {
    vec[i].first = i * 10;
    vec[i].second.key = i * 3;
  }
  StdVector<std::pair<uint64_t, SortElement>>::Reader reader(vec.begin(),
                                                             vec.end());

  omap.InitFromReader(reader);
}

TEST(OMap, Find) {
  size_t mapSize = 4321;
  size_t initSize = 1234;
  OMap<uint64_t, SortElement> omap(mapSize);
  StdVector<std::pair<uint64_t, SortElement>> vec(initSize);
  for (int i = 0; i < initSize; i++) {
    vec[i].first = i * 10;
    vec[i].second.key = i * 3;
  }
  StdVector<std::pair<uint64_t, SortElement>>::Reader reader(vec.begin(),
                                                             vec.end());

  omap.InitFromReader(reader);
  SortElement val;
  for (int r = 0; r < 1000; ++r) {
    uint64_t i = UniformRandom(initSize - 1);
    if (UniformRandomBit()) {
      ASSERT_TRUE(omap.find(10 * i, val));
      // printf("find %lu %lu\n", 10 * i, val.key);
      ASSERT_EQ(val.key, 3 * i);
    } else {
      ASSERT_FALSE(omap.find(10 * i + 1, val));
    }
  }
}

TEST(OMap, Insert) {
  size_t mapSize = 1e5;
  size_t initSize = 5e4;
  OMap<uint64_t, int64_t> omap(mapSize);
  StdVector<std::pair<uint64_t, int64_t>> vec(initSize);
  std::unordered_map<uint64_t, int64_t> map;
  for (int i = 0; i < initSize; i++) {
    vec[i].first = i * 10;
    vec[i].second = i * 3;
    map[i * 10] = i * 3;
  }

  StdVector<std::pair<uint64_t, int64_t>>::Reader reader(vec.begin(),
                                                         vec.end());

  omap.InitFromReader(reader);
  for (size_t r = 0; r < mapSize - initSize; ++r) {
    uint64_t i = UniformRandom(mapSize * 10);
    int64_t val = UniformRandom(mapSize * 3);
    bool res = omap.insert(i, val);
    if (map.find(i) != map.end()) {
      ASSERT_TRUE(res);
    } else {
      ASSERT_FALSE(res);
    }
    map[i] = val;
  }
  for (auto& p : map) {
    int64_t val;
    ASSERT_TRUE(omap.find(p.first, val));
    ASSERT_EQ(val, p.second);
  }
  for (size_t r = 0; r < 1000; ++r) {
    uint64_t i = UniformRandom(mapSize * 10);
    int64_t val;
    bool res = omap.find(i, val);
    if (map.find(i) != map.end()) {
      ASSERT_TRUE(res);
      ASSERT_EQ(val, map[i]);
    } else {
      ASSERT_FALSE(res);
    }
  }
}

TEST(OMap, FindPerf) {
  size_t mapSize = 1e7;
  size_t initSize = 1e7;
  OMap<uint64_t, int64_t, 9, uint32_t, uint32_t> omap(mapSize);
  StdVector<std::pair<uint64_t, int64_t>> vec(initSize);
  for (int i = 0; i < initSize; i++) {
    vec[i].first = i * 10;
    vec[i].second = i * 3;
  }
  StdVector<std::pair<uint64_t, int64_t>>::Reader reader(vec.begin(),
                                                         vec.end());
  // time initialization
  auto start = std::chrono::system_clock::now();
  omap.InitFromReader(reader);
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "omap init time: " << elapsed_seconds.count() << "s\n";
  auto start2 = std::chrono::system_clock::now();
  int64_t val;
  int round = 1e5;
  for (int r = 0; r < round; ++r) {
    uint64_t i = UniformRandom(initSize - 1);
    if (UniformRandomBit()) {
      ASSERT_TRUE(omap.find(10 * i, val));
      // printf("find %lu %lu\n", 10 * i, val.key);
      ASSERT_EQ(val, 3 * i);
    } else {
      ASSERT_FALSE(omap.find(10 * i + 1, val));
    }
  }
  auto end2 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds2 = end2 - start2;
  std::cout << "omap find time per access: "
            << elapsed_seconds2.count() / round * 1e6 << "us\n";
}

TEST(OMap, InsertPerf) {
  size_t mapSize = 1e7;
  size_t initSize = 1e5;
  OMap<uint64_t, int64_t> omap(mapSize);
  StdVector<std::pair<uint64_t, int64_t>> vec(initSize);
  std::unordered_map<uint64_t, int64_t> map;
  for (int i = 0; i < initSize; i++) {
    vec[i].first = i * 10;
    vec[i].second = i * 3;
    map[i * 10] = i * 3;
  }

  StdVector<std::pair<uint64_t, int64_t>>::Reader reader(vec.begin(),
                                                         vec.end());

  omap.InitFromReader(reader);
  int round = 1e5;
  auto start = std::chrono::system_clock::now();
  for (size_t r = 0; r < round; ++r) {
    uint64_t i = UniformRandom(mapSize * 10);
    int64_t val = UniformRandom(mapSize * 3);
    bool res = omap.insert(i, val);
    if (map.find(i) != map.end()) {
      ASSERT_TRUE(res);
    } else {
      ASSERT_FALSE(res);
    }
    map[i] = val;
  }
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "omap insert time per access: "
            << elapsed_seconds.count() / round * 1e6 << "us\n";
}

TEST(OMap, InsertAndFind) {
  printf("test omap\n");
  size_t mapSize = 1e5;
  size_t initSize = 5e4;
  OMap<uint64_t, int64_t> omap(mapSize);
  std::unordered_map<uint64_t, int64_t> map;
  for (int i = 0; i < initSize; i++) {
    map[i * 10] = i * 3;
  }

  std::function<std::pair<uint64_t, int64_t>(uint64_t)> readerFunc =
      [](uint64_t i) { return std::pair<uint64_t, int64_t>(i * 10, i * 3); };

  EM::VirtualVector::VirtualReader<std::pair<uint64_t, int64_t>> reader(
      initSize, readerFunc);
  uint64_t start, end;
  printf("init omap\n");
  omap.InitFromReader(reader);

  int round = 1e4;

  for (size_t r = 0; r < round; ++r) {
    uint64_t i = UniformRandom(mapSize * 10);
    int64_t val = UniformRandom(mapSize * 3);
    bool res = omap.insert(i, val);
    if (map.find(i) != map.end()) {
      if (!res) {
        printf("insert failed at round %lu, does not replace element\n", r);
        abort();
      }
    } else {
      if (res) {
        printf("insert failed at round %lu, found element that doesn't exist\n",
               r);
        abort();
      }
    }
    map[i] = val;
  }

  for (size_t r = 0; r < round; ++r) {
    uint64_t i = UniformRandom(mapSize * 10);
    int64_t val;
    bool res = omap.find(i, val);
    if (map.find(i) != map.end()) {
      if (!res) {
        printf("find failed at round %lu, does not find element\n", r);
        abort();
      }
      if (val != map[i]) {
        printf("find failed at round %lu, value not match\n", r);
        abort();
      }
    } else {
      if (res) {
        printf("find failed at round %lu, found element that doesn't exist\n",
               r);
        abort();
      }
    }
  }
}