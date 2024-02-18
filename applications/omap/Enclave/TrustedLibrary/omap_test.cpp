#include "../Enclave.h"
#include "Enclave_t.h"
////
#include <omp.h>

#include <functional>
#include <unordered_map>

#include "oram/omap.hpp"
#include "oram/par_omap.hpp"
#include "sgx_thread.h"
#include "sgx_trts.h"

using namespace ODSL;
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
    PathORAM::ORAM<uint64_t> oram(memSize, 6);
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

void testORAMInit() {
  uint64_t memSize = 4321;
  uint64_t size = 1234;

  PathORAM::ORAM<SortElement> oram(memSize);
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
    map[i] = i * 3;
  }

  std::function<std::pair<uint64_t, int64_t>(uint64_t)> readerFunc =
      [](uint64_t i) { return std::pair<uint64_t, int64_t>(i, i * 3); };

  EM::VirtualVector::VirtualReader<std::pair<uint64_t, int64_t>> reader(
      initSize, readerFunc);
  uint64_t start, end;
  printf("init omap\n");
  ocall_measure_time(&start);
  omap.InitFromReaderInPlace(reader);
  ocall_measure_time(&end);
  uint64_t timediff = end - start;
  printf("oram init time %d.%d s\n", timediff / 1'000'000'000,
         timediff % 1'000'000'000);
  int round = 1e5;
  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    uint64_t i = UniformRandom(mapSize);
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
#pragma omp parallel for
  for (size_t r = 0; r < round; ++r) {
    uint64_t i = UniformRandom(mapSize);
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

struct ETH_Addr {
  // 20 bytes

  uint32_t part[5];

  // define less operator
  bool operator<(const ETH_Addr& other) const {
    bool res = false;
    bool eq = true;
    for (int i = 0; i < 5; ++i) {
      res |= (eq & (part[i] < other.part[i]));
      eq &= part[i] == other.part[i];
    }
    return res;
  }

  bool operator==(const ETH_Addr& other) const {
    bool eq = true;
    for (int i = 0; i < 5; ++i) {
      eq &= (part[i] == other.part[i]);
    }
    return eq;
  }

  bool operator!=(const ETH_Addr& other) const { return !(*this == other); }
};

struct ERC20_Balance {
  uint64_t part[4];
};

void testOMapPerf() {
  printf("test omap perf with %d threads\n", TCS_NUM);
  printf("actual working thread max %d\n", omp_get_max_threads());
  size_t mapSize = 5e6;
  size_t initSize = 4e6;
  OMap<ETH_Addr, ERC20_Balance> omap(mapSize);

  std::function<std::pair<ETH_Addr, ERC20_Balance>(uint64_t)> readerFunc =
      [](uint64_t i) { return std::pair<ETH_Addr, ERC20_Balance>(); };

  EM::VirtualVector::VirtualReader<std::pair<ETH_Addr, ERC20_Balance>> reader(
      initSize, readerFunc);
  uint64_t start, end;
  printf("init omap of size %lu\n", mapSize);
  ocall_measure_time(&start);
  omap.InitFromReaderInPlace(reader);
  ocall_measure_time(&end);
  uint64_t timediff = end - start;
  printf("oram init time %d.%d s\n", timediff / 1'000'000'000,
         timediff % 1'000'000'000);
  int round = 1e6;
  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    ETH_Addr addr;
    ERC20_Balance balance;
    bool res = omap.insert(addr, balance);
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram insert time %d.%d us\n", timediff / round / 1'000,
         timediff / round % 1'000);

  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    ETH_Addr addr;
    ERC20_Balance balance;
    omap.find(addr, balance);
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram find time %d.%d us\n", timediff / round / 1'000,
         timediff / round % 1'000);

  ocall_measure_time(&start);
#pragma omp parallel for
  for (size_t r = 0; r < round; ++r) {
    ETH_Addr addr;
    ERC20_Balance balance;
    omap.find(addr, balance);
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram find with %d threads time %d.%d us\n", TCS_NUM,
         timediff / round / 1'000, timediff / round % 1'000);

  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    ETH_Addr addr;
    omap.erase(addr);
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram erase time %d.%d us\n", timediff / round / 1'000,
         timediff / round % 1'000);
}

void testParOMapPerf() {
  printf("test parallel omap perf with %d threads\n", TCS_NUM);
  printf("actual working thread max %d\n", omp_get_max_threads());
  size_t mapSize = 5e6;
  size_t initSize = 4e6;
  ParOMap<ETH_Addr, ERC20_Balance> omap(mapSize, 2);
  uint64_t start, end;
  omap.Init();
  uint64_t timediff = end - start;
  int round = 1e6;
  uint32_t batchSize = 1e3;
  ocall_measure_time(&start);
  for (size_t r = 0; r < round / batchSize; ++r) {
    std::vector<ETH_Addr> addr(batchSize);
    std::vector<ERC20_Balance> balance(batchSize);
    omap.insertBatch(addr.begin(), addr.end(), balance.begin());
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram insert time %d.%d us\n", timediff / round / 1'000,
         timediff / round % 1'000);

  ocall_measure_time(&start);
  for (size_t r = 0; r < round / batchSize; ++r) {
    std::vector<ETH_Addr> addr(batchSize);
    std::vector<ERC20_Balance> balance(batchSize);
    omap.findBatch(addr.begin(), addr.end(), balance.begin());
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram find time %d.%d us\n", timediff / round / 1'000,
         timediff / round % 1'000);
}

void ecall_omap() {
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
  size_t BackendSize = 1e10;
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

void ecall_omap_perf() {
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
  size_t BackendSize = 4e9;
  EM::Backend::g_DefaultBackend =
      new EM::Backend::MemServerBackend(BackendSize);
  try {
    testParOMapPerf();
    // testOMap();
    // testORAMInit();
    // testORAMReadWrite();
    // testBuildBottomUp();
  } catch (std::exception& e) {
    printf("exception: %s\n", e.what());
  }
  return;
}