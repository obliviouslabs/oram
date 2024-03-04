#include "oram/omap.hpp"

#include <gtest/gtest.h>

#include <unordered_map>

#include "external_memory/algorithm/sort_def.hpp"
#include "testutils.hpp"

using namespace ODSL;

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

TEST(OMap, InitInPlace) {
  OMap<uint64_t, SortElement> omap(512);
  StdVector<std::pair<uint64_t, SortElement>> vec(16);
  for (int i = 0; i < 16; i++) {
    vec[i].first = i * 10;
    vec[i].second.key = i * 3;
  }
  StdVector<std::pair<uint64_t, SortElement>>::Reader reader(vec.begin(),
                                                             vec.end());

  omap.InitFromReaderInPlace(reader);
}

TEST(OMap, Find) {
  for (int rr = 0; rr < 100; ++rr) {
    size_t mapSize = UniformRandom(2, 100000);
    size_t initSize = UniformRandom(1, mapSize);
    // std::cout << "mapSize: " << mapSize << " initSize: " << initSize
    //           << std::endl;
    OMap<uint64_t, int64_t> omap(mapSize);
    StdVector<std::pair<uint64_t, int64_t>> vec(initSize);
    for (int i = 0; i < initSize; i++) {
      vec[i].first = i * 10;
      vec[i].second = i * 3;
    }
    StdVector<std::pair<uint64_t, int64_t>>::Reader reader(vec.begin(),
                                                           vec.end());

    omap.InitFromReaderInPlace(reader);
    int64_t val;
    for (int r = 0; r < 1000; ++r) {
      uint64_t i = UniformRandom(initSize - 1);
      if (UniformRandomBit()) {
        ASSERT_TRUE(omap.find(10 * i, val));

        ASSERT_EQ(val, 3 * i);
      } else {
        ASSERT_FALSE(omap.find(10 * i + 1, val));
      }
    }
  }
}

TEST(OMap, Erase) {
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
  for (int i = 0; i < initSize; i++) {
    val.key = i * 3 + 2;
    omap.insert(i * 10 + 6, val);
  }
  for (int r = 0; r < 1000; ++r) {
    uint64_t i = UniformRandom(initSize - 1);
    if (vec[i].second.key != -1) {
      ASSERT_TRUE(omap.erase(10 * i));
      // printf("find %lu %lu\n", 10 * i, val.key);
      vec[i].second.key = -1;
    } else {
      ASSERT_FALSE(omap.erase(10 * i));
    }
  }
}

TEST(OMap, EraseAll) {
  size_t mapSize = 4321;
  size_t initSize = 4321;
  OMap<uint64_t, int> omap(mapSize);
  StdVector<std::pair<uint64_t, int>> vec(initSize);
  for (int i = 0; i < initSize; i++) {
    vec[i].first = i * 10;
    vec[i].second = i * 3;
  }
  StdVector<std::pair<uint64_t, int>>::Reader reader(vec.begin(), vec.end());

  omap.InitFromReader(reader);
  fisherYatesShuffle(vec.begin(), vec.end());
  for (int r = 0; r < initSize; ++r) {
    uint64_t i = vec[r].first;
    ASSERT_TRUE(omap.erase(i));
  }
  std::unordered_map<uint64_t, int> map;
  for (size_t r = 0; r < mapSize; ++r) {
    uint64_t i = UniformRandom(mapSize * 2);
    int val = UniformRandom(mapSize * 3);
    bool res = omap.insert(i, val);
    if (map.find(i) != map.end()) {
      ASSERT_TRUE(res);
    } else {
      ASSERT_FALSE(res);
    }
    map[i] = val;
  }
}

TEST(OMap, Update) {
  size_t mapSize = 54321;
  size_t initSize = 12345;
  OMap<uint64_t, SortElement> omap(mapSize);
  std::map<uint64_t, uint64_t> keyMap;
  StdVector<std::pair<uint64_t, SortElement>> vec(initSize);
  for (int i = 0; i < initSize; i++) {
    vec[i].first = i * 10;
    vec[i].second.key = i * 3;
    keyMap[i * 10] = i * 3;
  }
  StdVector<std::pair<uint64_t, SortElement>>::Reader reader(vec.begin(),
                                                             vec.end());

  omap.InitFromReader(reader);
  SortElement val;
  auto valUpdateFunc = [](SortElement& val) { ++val.key; };
  for (int r = 0; r < 10000; ++r) {
    uint64_t i = UniformRandom(initSize - 1);
    if (UniformRandomBit()) {
      ASSERT_TRUE(omap.update(10 * i, valUpdateFunc, val));
      keyMap[10 * i] += 1;
      ASSERT_EQ(val.key, keyMap[10 * i]);
    } else {
      ASSERT_FALSE(omap.update(10 * i + 1, valUpdateFunc));
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

TEST(OMap, InsertFromEmpty) {
  size_t mapSize = 1e5;

  OMap<uint64_t, int64_t> omap(mapSize);
  std::unordered_map<uint64_t, int64_t> map;
  omap.Init();
  for (size_t r = 0; r < mapSize; ++r) {
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

template <typename Node>
void setNodeHelper(Node& node, short numChildren, short firstKey) {
  for (short i = 0; i < numChildren; ++i) {
    node.children[i] = {(uint64_t)(10 * (firstKey + i)),
                        (uint64_t)(10 * (firstKey + i))};
    if (i > 0) node.keys[i - 1] = firstKey + i;
  }
  node.numChildren = numChildren;
}

void testRedistributeHelper(int leftNum, int rightNum) {
  using OMap_ = OMap<int64_t, int64_t, 9>;
  typename OMap_::BPlusNode_ nodeLeft, nodeRight;
  int64_t parentKey = leftNum;
  setNodeHelper(nodeLeft, leftNum, 0);
  setNodeHelper(nodeRight, rightNum, leftNum);
  std::cout << "left: " << nodeLeft << std::endl;
  std::cout << "right: " << nodeRight << std::endl;
  std::cout << "parentKey: " << parentKey << std::endl;
  OMap_::redistributeOrCoalesceNode(true, nodeLeft, nodeRight, parentKey);
  std::cout << "after redistribute: " << std::endl;
  std::cout << "left: " << nodeLeft << std::endl;
  std::cout << "right: " << nodeRight << std::endl;
  std::cout << "parentKey: " << parentKey << std::endl << std::endl;
}

TEST(OMap, TestRedistribute) {
  testRedistributeHelper(5, 5);
  testRedistributeHelper(4, 6);
  testRedistributeHelper(6, 4);
  testRedistributeHelper(4, 9);
  testRedistributeHelper(9, 4);
  testRedistributeHelper(4, 5);
  testRedistributeHelper(5, 4);
}

TEST(OMap, FindPerf) {
  size_t mapSize = 1UL << 26;
  size_t initSize = 1UL << 26;
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
  size_t BackendSize = 1e8;
  EM::Backend::g_DefaultBackend =
      new EM::Backend::MemServerBackend(BackendSize);
  OMap<uint64_t, int64_t, 9, uint64_t, uint64_t> omap(mapSize, 1UL << 35);
  std::function<std::pair<uint64_t, int64_t>(uint64_t)> readerFunc =
      [](uint64_t i) { return std::pair<uint64_t, int64_t>(i * 10, i * 3); };

  EM::VirtualVector::VirtualReader<std::pair<uint64_t, int64_t>> reader(
      initSize, readerFunc);
  // time initialization
  auto start = std::chrono::system_clock::now();
  omap.InitFromReaderInPlace(reader);
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "omap init time: " << elapsed_seconds.count() << "s\n";
  auto start2 = std::chrono::system_clock::now();
  int64_t val;
  int round = 1e5;
  for (int r = 0; r < round; ++r) {
    uint64_t i = UniformRandom(initSize - 1);
    if (UniformRandomBit()) {
      if (!omap.find(10 * i, val)) {
        printf("find failed at round %d\n", r);
        abort();
      }

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
  size_t initSize = 1e4;
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
  size_t BackendSize = 1e10;
  EM::Backend::g_DefaultBackend =
      new EM::Backend::MemServerBackend(BackendSize);
  OMap<uint64_t, int64_t> omap(mapSize);
  StdVector<std::pair<uint64_t, int64_t>> vec(initSize);
  std::unordered_map<uint64_t, int64_t> map;
  for (int i = 0; i < initSize; i++) {
    vec[i].first = i;
    vec[i].second = i * 3;
    map[i] = i * 3;
  }

  StdVector<std::pair<uint64_t, int64_t>>::Reader reader(vec.begin(),
                                                         vec.end());

  omap.InitFromReader(reader);
  int round = 1e5;
  auto start = std::chrono::system_clock::now();
  for (size_t r = 0; r < round; ++r) {
    uint64_t i = UniformRandom(mapSize);
    int64_t val = UniformRandom(mapSize * 3);
    bool res = omap.insert(i, val);
    if (map.find(i) != map.end()) {
      if (!res) {
        printf("insert failed at round %lu, does not replace element\n", r);
      }
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

TEST(OMap, ErasePerf) {
  size_t mapSize = 1e7;
  size_t initSize = 1e5;
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
  size_t BackendSize = 1e10;
  EM::Backend::g_DefaultBackend =
      new EM::Backend::MemServerBackend(BackendSize);
  OMap<uint64_t, int64_t> omap(mapSize);
  StdVector<std::pair<uint64_t, int64_t>> vec(initSize);
  for (int i = 0; i < initSize; i++) {
    vec[i].first = i;
    vec[i].second = i * 3;
  }

  StdVector<std::pair<uint64_t, int64_t>>::Reader reader(vec.begin(),
                                                         vec.end());

  omap.InitFromReader(reader);
  int round = 1e5;
  auto start = std::chrono::system_clock::now();
  for (size_t r = 0; r < round; ++r) {
    uint64_t i = UniformRandom(mapSize);
    bool res = omap.erase(i);
  }
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "omap erase time per access: "
            << elapsed_seconds.count() / round * 1e6 << "us\n";
}

TEST(OMap, InsertAndFind) {
  printf("test omap\n");
  size_t mapSize = 1e5;
  size_t initSize = 5e4;
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
  size_t BackendSize = 1e10;
  EM::Backend::g_DefaultBackend =
      new EM::Backend::MemServerBackend(BackendSize);
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

TEST(OMap, Mixed) {
  printf("test omap\n");
  size_t mapSize = 1e5;
  size_t initSize = 1e4;
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
  size_t BackendSize = 1e10;
  EM::Backend::g_DefaultBackend =
      new EM::Backend::MemServerBackend(BackendSize);
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
  omap.InitFromReader(reader);

  int round = 1e6;

  for (size_t r = 0; r < round; ++r) {
    uint64_t op = UniformRandom(2);
    bool isDummy = !(UniformRandom() % 5);
    if (op == 0) {  // insert
      uint64_t i = UniformRandom(mapSize * 10);
      int64_t val = UniformRandom(mapSize * 3);
      bool res = omap.insert(i, val, isDummy);
      if (isDummy) {
        continue;
      }
      if (map.find(i) != map.end()) {
        if (!res) {
          printf("insert failed at round %lu, does not replace element\n", r);
          abort();
        }
      } else {
        if (res) {
          printf(
              "insert failed at round %lu, found element that doesn't exist\n",
              r);
          abort();
        }
      }
      map[i] = val;
    } else if (op == 1) {  // find
      uint64_t i = UniformRandom(mapSize * 10);
      int64_t val;
      bool res = omap.find(i, val, isDummy);
      if (isDummy) {
        continue;
      }
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
    } else if (op == 2) {  // erase
      uint64_t i = UniformRandom(mapSize * 10);
      bool res = omap.erase(i, isDummy);
      if (isDummy) {
        continue;
      }
      if (map.find(i) != map.end()) {
        ASSERT_TRUE(res);
        map.erase(i);
      } else {
        ASSERT_FALSE(res);
      }
    }
  }
}