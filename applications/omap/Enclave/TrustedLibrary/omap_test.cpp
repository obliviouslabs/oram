#include "../Enclave.h"
#include "Enclave_t.h"
////
#include <functional>
#include <unordered_map>

#include "oram/omap.hpp"

using namespace ORAM;
EM::Backend::MemServerBackend* EM::Backend::g_DefaultBackend = nullptr;

#define ASSERT_EQ(a, b)                             \
  if ((a) != (b)) {                                 \
    printf("assert failed at line %d\n", __LINE__); \
    printf("%lu != %lu\n", (a), (b));               \
    abort();                                        \
  }

void testORAMReadWrite() {
  for (int memSize : {1, 3, 5, 7, 9, 33, 40, 55, 127, 129, 543, 678, 1023, 1025,
                      2000, 10000, 50000}) {
    PathORAM::PathORAM<uint64_t> oram(memSize, 6);
    std::vector<uint64_t> posMap(memSize);
    std::vector<uint64_t> valMap(memSize);
    for (uint64_t i = 0; i < memSize; i++) {
      uint64_t val = UniformRandom();
      uint64_t pos = oram.Write(i, val);
      posMap[i] = pos;
      valMap[i] = val;

      // printf("write %lu %lu to pos %lu\n", i, val, pos);
    }
    for (int r = 0; r < 7; ++r) {
      for (uint64_t i = 0; i < memSize; i++) {
        uint64_t val;
        // printf("read %lu %lu at pos %lu\n", i, val, posMap[i]);
        uint64_t pos = oram.Read(posMap[i], i, val);
        posMap[i] = pos;
        ASSERT_EQ(val, valMap[i]);
      }
    }
  }
}

void testBuildBottomUp() {
  for (size_t leafCount = 2; leafCount < 1000; ++leafCount) {
    int totalLevel = GetLogBaseTwo(leafCount - 1) + 2;
    for (int cacheLevel = 1; cacheLevel <= totalLevel; ++cacheLevel) {
      for (int packLevel = 1; packLevel <= totalLevel; ++packLevel) {
        // printf("leafCount: %lu, cacheLevel: %d, packLevel: %d\n", leafCount,
        //        cacheLevel, packLevel);
        size_t size = leafCount * 2 - 1;
        std::vector<uint64_t> heap(size);
        std::function<uint64_t(uint64_t&, const uint64_t&, const uint64_t&)>
            reduceFunc =
                [](uint64_t& val, const uint64_t& i,
                   const uint64_t& j) -> uint64_t { return val = i + j; };
        std::function<uint64_t(uint64_t&, size_t)> leafFunc =
            [](uint64_t& i, size_t idx) { return i = idx; };
        ASSERT_EQ(HeapTree<uint64_t>::BuildBottomUp<uint64_t>(
                      heap.begin(), heap.end(), reduceFunc, leafFunc,
                      cacheLevel, packLevel),
                  leafCount * (leafCount - 1) / 2);
      }
    }
  }
}

void testORAMInit() {
  uint64_t memSize = 4321;
  uint64_t size = 1234;

  PathORAM::PathORAM<SortElement> oram(memSize);
  std::vector<uint64_t> valMap(size);
  StdVector<SortElement> vec(size);
  for (int i = 0; i < size; ++i) {
    valMap[i] = UniformRandom();
    vec[i].key = valMap[i];
  }

  StdVector<UidBlock<uint64_t>> posMap(size);
  StdVector<SortElement>::Reader reader(vec.begin(), vec.end());
  StdVector<UidBlock<uint64_t>>::Writer posMapWriter(posMap.begin(),
                                                     posMap.end());

  oram.InitFromReader(reader, posMapWriter);
  // for (size_t i = 0; i < size; ++i) {
  //   printf("posMap[%lu] = %lu\n", posMap[i].uid, posMap[i].data);
  // }

  for (int r = 0; r < 2; ++r) {
    for (uint64_t i = 0; i < size; i++) {
      SortElement val;
      uint64_t pos = oram.Read(posMap[i].data, i, val);
      posMap[i].data = pos;
      ASSERT_EQ(val.key, valMap[i]);
    }
  }
}

void testOMap() {
  printf("test omap\n");
  size_t mapSize = 1e6;
  size_t initSize = 1e6;
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
  ocall_measure_time(&start);
  omap.InitFromReader(reader);
  ocall_measure_time(&end);
  uint64_t timediff = end - start;
  printf("oram init time %d.%d s\n", timediff / 1'000'000'000,
         timediff % 1'000'000'000);
  int round = 1e4;
  ocall_measure_time(&start);
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
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram insert time %d.%d us\n", timediff / round / 1'000,
         timediff / round % 1'000);

  ocall_measure_time(&start);
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
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram find time %d.%d us\n", timediff / round / 1'000,
         timediff / round % 1'000);
}

void ecall_omap_perf() {
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
  size_t BackendSize = 1e9;
  EM::Backend::g_DefaultBackend =
      new EM::Backend::MemServerBackend(BackendSize);
  try {
    testOMap();
    // testORAMInit();
    // testORAMReadWrite();
    // testBuildBottomUp();
  } catch (std::exception& e) {
    printf("exception: %s\n", e.what());
  }
  return;
}