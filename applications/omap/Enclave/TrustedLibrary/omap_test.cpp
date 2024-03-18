#include "../Enclave.h"
#include "Enclave_t.h"
////
#include <omp.h>

#include <functional>
#include <unordered_map>

#include "oram/omap.hpp"
#include "oram/par_omap.hpp"
#include "oram/recursive_oram.hpp"
#include "sgx_thread.h"
#include "sgx_trts.h"

#define MB << 20

using namespace ODSL;
EM::Backend::MemServerBackend* EM::Backend::g_DefaultBackend = nullptr;

#define ASSERT_EQ(a, b)                             \
  if ((a) != (b)) {                                 \
    printf("assert failed at line %d\n", __LINE__); \
    printf("%lu != %lu\n", (a), (b));               \
    abort();                                        \
  }

void testOmpSpeedup() {
  size_t maxThread = omp_get_max_threads();
  uint64_t start, end;
  int64_t timediffOneThread, timediffMultipleThread;
  printf("bitonic sort benchmark\n");
  ocall_measure_time(&start);
  for (int i = 0; i < maxThread; ++i) {
    std::vector<uint64_t> vec(65546);
    for (int r = 0; r < 100; ++r) {
      EM::Algorithm::BitonicSort(vec);
    }
  }
  ocall_measure_time(&end);
  timediffOneThread = end - start;
  printf("one thread %f s\n", (double)timediffOneThread * 1e-9);

  ocall_measure_time(&start);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < maxThread; ++i) {
    std::vector<uint64_t> vec(65546);
    for (int r = 0; r < 100; ++r) {
      EM::Algorithm::BitonicSort(vec);
    }
  }
  ocall_measure_time(&end);
  timediffMultipleThread = end - start;
  printf("%d thread %f s\n", omp_get_max_threads(),
         (double)timediffMultipleThread * 1e-9);
  printf("speedup for bitonic sort benchmark (memory intensive task) = %f\n",
         (double)timediffOneThread / timediffMultipleThread);

  printf("\nReverse benchmark\n");
  ocall_measure_time(&start);
  for (int i = 0; i < maxThread; ++i) {
    std::vector<uint64_t> vec(65546);
    for (int t = 0; t < 50000; ++t) {
      std::reverse(vec.begin(), vec.end());
    }
  }
  ocall_measure_time(&end);
  timediffOneThread = end - start;
  printf("one thread %f s\n", (double)timediffOneThread * 1e-9);

  ocall_measure_time(&start);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < maxThread; ++i) {
    std::vector<uint64_t> vec(65546);
    for (int t = 0; t < 50000; ++t) {
      std::reverse(vec.begin(), vec.end());
    }
  }
  ocall_measure_time(&end);
  timediffMultipleThread = end - start;
  printf("%d thread %f s\n", omp_get_max_threads(),
         (double)timediffMultipleThread * 1e-9);
  printf("speedup for reverse benchmark (memory intensive task) = %f\n",
         (double)timediffOneThread / timediffMultipleThread);

  printf("\nFloating point benchmark\n");
  ocall_measure_time(&start);
  double totalSum;
  for (int i = 0; i < maxThread; ++i) {
    double sum = 1.0 + UniformRandom() % 2;
    for (uint64_t i = 0; i < 1e8; ++i) {
      sum += 1.0 / sum;
    }
    totalSum += sum;
  }
  ocall_measure_time(&end);
  printf("total sum = %lu\n", totalSum);
  timediffOneThread = end - start;
  printf("one thread %f s\n", (double)timediffOneThread * 1e-9);

  ocall_measure_time(&start);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < maxThread; ++i) {
    double sum = 1.0 + UniformRandom() % 2;
    for (uint64_t i = 0; i < 1e8; ++i) {
      sum += 1.0 / sum;
    }
    totalSum += sum;
  }
  ocall_measure_time(&end);
  printf("total sum = %lu\n", totalSum);
  timediffMultipleThread = end - start;
  printf("%d thread %f s\n", omp_get_max_threads(),
         (double)timediffMultipleThread * 1e-9);
  printf("speedup for floating point benchmark (cpu intensive task) = %f\n",
         (double)timediffOneThread / timediffMultipleThread);

  printf("\nHash benchmark\n");
  ocall_measure_time(&start);
  struct Key {
    uint8_t key[20];
  };
  Key k;
  uint8_t salt[16];
  for (int i = 0; i < maxThread; ++i) {
    for (int i = 0; i < 1e5; ++i) {
      uint64_t hash = secure_hash_with_salt(k, salt);
    }
  }
  ocall_measure_time(&end);
  timediffOneThread = end - start;
  printf("one thread %f s\n", (double)timediffOneThread * 1e-9);

  ocall_measure_time(&start);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < maxThread; ++i) {
    for (int i = 0; i < 1e5; ++i) {
      uint64_t hash = secure_hash_with_salt(k, salt);
    }
  }
  ocall_measure_time(&end);
  timediffMultipleThread = end - start;
  printf("%d thread %f s\n", omp_get_max_threads(),
         (double)timediffMultipleThread * 1e-9);
  printf("speedup for hash benchmark (cpu intensive task) = %f\n",
         (double)timediffOneThread / timediffMultipleThread);

  printf("\nSgx PRNG benchmark\n");
  ocall_measure_time(&start);
  uint64_t sumRand = 0;
  for (int i = 0; i < maxThread; ++i) {
    uint64_t localRand = 0;
    for (uint64_t i = 0; i < (uint64_t)1e6; ++i) {
      uint64_t rd;
      read_rand((uint8_t*)&rd, sizeof(uint64_t));
      localRand += rd;
    }
    sumRand += localRand;
  }

  ocall_measure_time(&end);
  timediffOneThread = end - start;
  printf("one thread %f s\n", (double)timediffOneThread * 1e-9);

  ocall_measure_time(&start);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < maxThread; ++i) {
    uint64_t localRand = 0;
    for (uint64_t i = 0; i < (uint64_t)1e6; ++i) {
      uint64_t rd;
      read_rand((uint8_t*)&rd, sizeof(uint64_t));
      localRand += rd;
    }
    sumRand += localRand;
  }
  ocall_measure_time(&end);
  printf("total sum = %lu\n", sumRand);
  timediffMultipleThread = end - start;
  printf("%d thread %f s\n", omp_get_max_threads(),
         (double)timediffMultipleThread * 1e-9);
  printf("speedup for Sgx PRNG benchmark = %f\n",
         (double)timediffOneThread / timediffMultipleThread);

  printf("\nHeapTree benchmark\n");
  ocall_measure_time(&start);

  for (int i = 0; i < maxThread; ++i) {
    HeapTree<uint64_t> tree;
    uint64_t size = 2048;
    tree.Init(size);
    std::vector<uint64_t> path(tree.GetDepth());
    for (int r = 0; r < 1e7; ++r) {
      tree.ReadPath(r % size, path.begin());
      tree.WritePath(r % size, path.begin());
    }
  }

  ocall_measure_time(&end);
  timediffOneThread = end - start;
  printf("one thread %f s\n", (double)timediffOneThread * 1e-9);

  ocall_measure_time(&start);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < maxThread; ++i) {
    HeapTree<uint64_t> tree;
    uint64_t size = 2048;
    tree.Init(size);
    std::vector<uint64_t> path(tree.GetDepth());
    for (int r = 0; r < 1e7; ++r) {
      tree.ReadPath(r % size, path.begin());
      tree.WritePath(r % size, path.begin());
    }
  }
  ocall_measure_time(&end);
  timediffMultipleThread = end - start;
  printf("%d thread %f s\n", omp_get_max_threads(),
         (double)timediffMultipleThread * 1e-9);
  printf("one thread %f s\n", (double)timediffOneThread * 1e-9);
  printf("speedup for heap tree benchmark = %f\n",
         (double)timediffOneThread / timediffMultipleThread);

  printf("\nLambda benchmark\n");
  ocall_measure_time(&start);
  for (int i = 0; i < maxThread; ++i) {
    uint64_t localRand = UniformRandom();
    auto inc = [&]() { localRand = localRand * 2 - 1; };
    for (uint64_t i = 0; i < (uint64_t)1e9; ++i) {
      inc();
    }
    sumRand += localRand;
  }

  ocall_measure_time(&end);
  timediffOneThread = end - start;
  printf("one thread %f s\n", (double)timediffOneThread * 1e-9);

  ocall_measure_time(&start);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < maxThread; ++i) {
    uint64_t localRand = UniformRandom();
    auto inc = [&]() { localRand = localRand * 2 - 1; };
    for (uint64_t i = 0; i < (uint64_t)1e9; ++i) {
      inc();
    }
    sumRand += localRand;
  }
  ocall_measure_time(&end);
  printf("total sum = %lu\n", sumRand);
  timediffMultipleThread = end - start;
  printf("%d thread %f s\n", omp_get_max_threads(),
         (double)timediffMultipleThread * 1e-9);
  printf("speedup for lambda benchmark = %f\n",
         (double)timediffOneThread / timediffMultipleThread);

  printf("\nCircuit ORAM update benchmark\n");
  ocall_measure_time(&start);

  for (int i = 0; i < maxThread; ++i) {
    CircuitORAM::ORAM<TestElement> oram;
    uint64_t size = 65536;
    oram.SetSize(size);
    uint64_t out;
    for (int r = 0; r < 1e5; ++r) {
      oram.Update(UniformRandom() % size, 0, [](TestElement& val) {
        val.key++;
        return true;
      });
    }
  }

  ocall_measure_time(&end);
  timediffOneThread = end - start;
  printf("one thread %f s\n", (double)timediffOneThread * 1e-9);

  ocall_measure_time(&start);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < maxThread; ++i) {
    CircuitORAM::ORAM<TestElement> oram;
    uint64_t size = 65536;
    oram.SetSize(size);
    uint64_t out;
    for (int r = 0; r < 1e5; ++r) {
      oram.Update(UniformRandom() % size, 0, [](TestElement& val) {
        val.key++;
        return true;
      });
    }
  }
  ocall_measure_time(&end);
  timediffMultipleThread = end - start;
  printf("%d thread %f s\n", omp_get_max_threads(),
         (double)timediffMultipleThread * 1e-9);
  printf("one thread %f s\n", (double)timediffOneThread * 1e-9);
  printf("speedup for circuit oram benchmark = %f\n",
         (double)timediffOneThread / timediffMultipleThread);

  printf("\nMulti Circuit ORAM update benchmark\n");
  ocall_measure_time(&start);

  for (int i = 0; i < maxThread; ++i) {
    std::vector<CircuitORAM::ORAM<TestElement>> orams(8);
    uint64_t size = 8192;
    for (int i = 0; i < 8; ++i) {
      orams[i].SetSize(size);
    }
    uint64_t out;
    for (int r = 0; r < 1e5; ++r) {
      orams[r % 8].Update(UniformRandom() % size, 0, [](TestElement& val) {
        val.key++;
        return true;
      });
    }
  }

  ocall_measure_time(&end);
  timediffOneThread = end - start;
  printf("one thread %f s\n", (double)timediffOneThread * 1e-9);

  ocall_measure_time(&start);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < maxThread; ++i) {
    std::vector<CircuitORAM::ORAM<TestElement>> orams(8);
    uint64_t size = 8192;
    for (int i = 0; i < 8; ++i) {
      orams[i].SetSize(size);
    }
    uint64_t out;
    for (int r = 0; r < 1e5; ++r) {
      orams[r % 8].Update(UniformRandom() % size, 0, [](TestElement& val) {
        val.key++;
        return true;
      });
    }
  }
  ocall_measure_time(&end);
  timediffMultipleThread = end - start;
  printf("%d thread %f s\n", omp_get_max_threads(),
         (double)timediffMultipleThread * 1e-9);

  printf("speedup for multiple (round robin) circuit oram benchmark = %f\n",
         (double)timediffOneThread / timediffMultipleThread);

  printf("\nLinear ORAM update benchmark\n");
  ocall_measure_time(&start);

  for (int i = 0; i < maxThread; ++i) {
    LinearORAM::ORAM<TestElement> oram;
    uint64_t size = 256;
    oram.SetSize(size);
    uint64_t out;
    for (int r = 0; r < 1e5; ++r) {
      oram.Update(0, [](TestElement& val) {
        val.key++;
        return true;
      });
    }
  }

  ocall_measure_time(&end);
  timediffOneThread = end - start;
  printf("one thread %f s\n", (double)timediffOneThread * 1e-9);

  ocall_measure_time(&start);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < maxThread; ++i) {
    LinearORAM::ORAM<TestElement> oram;
    uint64_t size = 256;
    oram.SetSize(size);
    uint64_t out;
    for (int r = 0; r < 1e5; ++r) {
      oram.Update(0, [](TestElement& val) {
        val.key++;
        return true;
      });
    }
  }
  ocall_measure_time(&end);
  timediffMultipleThread = end - start;
  printf("%d thread %f s\n", omp_get_max_threads(),
         (double)timediffMultipleThread * 1e-9);

  printf("speedup for linear oram benchmark = %f\n",
         (double)timediffOneThread / timediffMultipleThread);

  printf("\nRecursive ORAM find benchmark\n");

  std::vector<RecursiveORAM<uint64_t>> orams(maxThread);

  for (int i = 0; i < maxThread; ++i) {
    orams[i].SetSize(4096);
    orams[i].InitDefault(0);
  }
  ocall_measure_time(&start);
  for (int i = 0; i < maxThread; ++i) {
    RecursiveORAM<uint64_t>& oram = orams[i];
    uint64_t size = 4096;
    uint64_t out;
    for (int r = 0; r < 1e5; ++r) {
      oram.Read(UniformRandom() % size, out);
    }
  }

  ocall_measure_time(&end);
  timediffOneThread = end - start;
  printf("one thread %f s\n", (double)timediffOneThread * 1e-9);

  ocall_measure_time(&start);

#pragma omp parallel for schedule(static)
  for (int i = 0; i < maxThread; ++i) {
    RecursiveORAM<uint64_t>& oram = orams[i];
    uint64_t size = 4096;
    uint64_t out;
    for (int r = 0; r < 1e5; ++r) {
      oram.Read(UniformRandom() % size, out);
    }
  }
  ocall_measure_time(&end);

  timediffMultipleThread = end - start;
  printf("%d thread %f s\n", omp_get_max_threads(),
         (double)timediffMultipleThread * 1e-9);

  printf("speedup for recursive oram benchmark = %f\n",
         (double)timediffOneThread / timediffMultipleThread);
}

void testORAMReadWrite() {
  for (int memSize : {1, 3, 5, 7, 9, 33, 40, 55, 127, 129, 543, 678, 1023, 1025,
                      2000, 10000, 50000}) {
    CircuitORAM::ORAM<uint64_t> oram(memSize, 10 MB);
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

  CircuitORAM::ORAM<TestElement> oram(memSize);
  std::vector<uint64_t> valMap(size);
  StdVector<TestElement> vec(size);
  for (int i = 0; i < size; ++i) {
    valMap[i] = UniformRandom();
    vec[i].key = valMap[i];
  }

  StdVector<UidBlock<uint64_t>> posMap(size);
  StdVector<TestElement>::Reader reader(vec.begin(), vec.end());
  StdVector<UidBlock<uint64_t>>::Writer posMapWriter(posMap.begin(),
                                                     posMap.end());

  oram.InitFromReader(reader, posMapWriter);
  // for (size_t i = 0; i < size; ++i) {
  //   printf("posMap[%lu] = %lu\n", posMap[i].uid, posMap[i].data);
  // }

  for (int r = 0; r < 2; ++r) {
    for (uint64_t i = 0; i < size; i++) {
      TestElement val;
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
  OHashMap<uint64_t, int64_t, true> omap(mapSize);
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
  omap.InitFromReader(reader);
  ocall_measure_time(&end);
  uint64_t timediff = end - start;
  printf("oram init time %f s\n", timediff * 1e-9);
  int round = 1e5;
  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    uint64_t i = UniformRandom(mapSize);
    int64_t val = UniformRandom(mapSize * 3);
    bool res = omap.InsertOblivious(i, val);
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
  printf("oram insert time %f us\n", timediff * 1e-3 / round);

  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    uint64_t i = UniformRandom(mapSize);
    int64_t val;
    bool res = omap.Find(i, val);
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
  printf("oram find time %f us\n", timediff * 1e-3 / round);
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
  OHashMap<ETH_Addr, ERC20_Balance, true> omap(mapSize);

  std::function<std::pair<ETH_Addr, ERC20_Balance>(uint64_t)> readerFunc =
      [](uint64_t i) { return std::pair<ETH_Addr, ERC20_Balance>(); };

  EM::VirtualVector::VirtualReader<std::pair<ETH_Addr, ERC20_Balance>> reader(
      initSize, readerFunc);
  uint64_t start, end;
  printf("init omap of size %lu\n", mapSize);
  ocall_measure_time(&start);
  omap.InitFromReader(reader);
  ocall_measure_time(&end);
  uint64_t timediff = end - start;
  printf("oram init time %f s\n", timediff * 1e-9);
  int round = 1e5;
  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    ETH_Addr addr;
    ERC20_Balance balance;
    bool res = omap.Insert(addr, balance);
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram insert time %f us\n", timediff * 1e-3 / round);

  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    ETH_Addr addr;
    ERC20_Balance balance;
    omap.Find(addr, balance);
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram find time %f us\n", timediff * 1e-3 / round);

  ocall_measure_time(&start);
#pragma omp parallel for
  for (size_t r = 0; r < round; ++r) {
    ETH_Addr addr;
    ERC20_Balance balance;
    omap.Find(addr, balance);
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram find with %d threads time %d.%d us\n", TCS_NUM,
         timediff / round / 1'000, timediff / round % 1'000);

  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    ETH_Addr addr;
    omap.Erase(addr);
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram erase time %f us\n", timediff * 1e-3 / round);
}

void testOHashMapPerf(size_t mapSize = 5e6) {
  int round = 1e5;
  size_t initSize = mapSize - round;
  OHashMap<ETH_Addr, ERC20_Balance, true, uint32_t> omap(mapSize);

  std::function<std::pair<ETH_Addr, ERC20_Balance>(uint64_t)> readerFunc =
      [](uint64_t i) { return std::pair<ETH_Addr, ERC20_Balance>(); };

  EM::VirtualVector::VirtualReader<std::pair<ETH_Addr, ERC20_Balance>> reader(
      initSize, readerFunc);
  uint64_t start, end;
  printf("mapSize = %u, threadCount = %d, batchSize = %u\n", mapSize, 1, 1);
  ocall_measure_time(&start);
  omap.InitFromReader(reader);
  ocall_measure_time(&end);
  uint64_t timediff = end - start;
  printf("oram init time %f s\n", timediff * 1e-9);

  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    ETH_Addr addr;
    addr.part[0] = r;
    ERC20_Balance balance;
    bool res = omap.Insert(addr, balance);
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram insert time %f us\n", timediff * 1e-3 / round);

  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    ETH_Addr addr;
    addr.part[0] = r;
    ERC20_Balance balance;
    omap.Find(addr, balance);
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram find time %f us\n", timediff * 1e-3 / round);

  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    ETH_Addr addr;
    omap.Erase(addr);
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram erase time %f us\n", timediff * 1e-3 / round);
}

template <const size_t size>
struct Bytes {
  uint8_t data[size];
};

void testOHashMapPerfSignal(size_t mapSize = 5e6) {
  int round = 1e5;
  size_t initSize = mapSize - round;
  OHashMap<uint64_t, Bytes<240>, true, uint32_t> omap(mapSize);

  std::function<std::pair<uint64_t, Bytes<240>>(uint64_t)> readerFunc =
      [](uint64_t i) { return std::pair<uint64_t, Bytes<240>>(); };

  EM::VirtualVector::VirtualReader<std::pair<uint64_t, Bytes<240>>> reader(
      initSize, readerFunc);
  uint64_t start, end;
  printf("mapSize = %u, threadCount = %d, batchSize = %u\n", mapSize, 1, 1);
  ocall_measure_time(&start);
  omap.InitFromReader(reader);
  ocall_measure_time(&end);
  uint64_t timediff = end - start;
  printf("oram init time %f s\n", timediff * 1e-9);

  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    uint64_t addr;
    addr = r;
    Bytes<240> balance;
    bool res = omap.Insert(addr, balance);
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram insert time %f us\n", timediff * 1e-3 / round);

  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    uint64_t addr;
    addr = r;
    Bytes<240> balance;
    omap.Find(addr, balance);
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram find time %f us\n", timediff * 1e-3 / round);

  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    uint64_t addr;
    addr = r;
    omap.Erase(addr);
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram erase time %f us\n", timediff * 1e-3 / round);
}

void testRecursiveORAMPerf() {
  printf("test recursive oram perf with %d threads\n", TCS_NUM);
  printf("actual working thread max %d\n", omp_get_max_threads());
  size_t mapSize = 5e6;
  struct AddrBalance {
    ETH_Addr addr;
    ERC20_Balance balance;
  };
  RecursiveORAM<AddrBalance, uint32_t> roram(mapSize);

  std::function<AddrBalance(uint64_t)> readerFunc = [](uint64_t i) {
    return AddrBalance();
  };

  EM::VirtualVector::VirtualReader<AddrBalance> reader(mapSize, readerFunc);
  uint64_t start, end;
  printf("init recursive oram of size %lu\n", mapSize);
  ocall_measure_time(&start);
  roram.InitFromReader(reader);
  ocall_measure_time(&end);
  uint64_t timediff = end - start;
  printf("oram init time %f s\n", timediff * 1e-9);
  int round = 1e5;
  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    AddrBalance pr;
    uint64_t addr = UniformRandom(mapSize - 1);
    roram.Read(addr, pr);
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("recursive oram access time %d.%d us\n", timediff / round / 1'000,
         timediff / round % 1'000);
}

void testParOMapPerf(size_t mapSize = 5e6,
                     int threadCount = omp_get_max_threads()) {
  size_t initSize = mapSize - 1e5;
  ParOMap<ETH_Addr, ERC20_Balance, uint32_t> omap(mapSize, threadCount / 2);
  std::function<std::pair<ETH_Addr, ERC20_Balance>(uint64_t)> readerFunc =
      [](uint64_t i) {
        std::pair<ETH_Addr, ERC20_Balance> pr;
        pr.first.part[0] = i;
        return pr;
      };

  EM::VirtualVector::VirtualReader<std::pair<ETH_Addr, ERC20_Balance>> reader(
      initSize, readerFunc);
  uint64_t start, end;
  // printf("init omap of size %lu\n", mapSize);
  ocall_measure_time(&start);
  // omap.InitFromReader(reader);
  omap.Init();
  ocall_measure_time(&end);
  uint64_t initTimediff = end - start;
  for (uint32_t batchSize : {100, 200, 500, 1000, 2000, 5000, 10000}) {
    printf("mapSize = %u, threadCount = %d, batchSize = %u\n", mapSize,
           threadCount, batchSize);
    printf("oram init time %f s\n", initTimediff * 1e-9);
    int round = 5e5;
    // ocall_measure_time(&start);
    // for (size_t r = 0; r < round / batchSize; ++r) {
    //   std::vector<ETH_Addr> addr(batchSize);
    //   std::vector<ERC20_Balance> balance(batchSize);
    //   for (int i = 0; i < batchSize; ++i) {
    //     addr[i].part[0] = initSize + r * batchSize + i;
    //   }
    //   omap.InsertParBatch(addr.begin(), addr.end(), balance.begin());
    // }
    // ocall_measure_time(&end);
    // uint64_t timediff = end - start;
    // printf("oram insert time %f us\n", timediff * 1e-3 / round);
    ocall_measure_time(&start);
    for (size_t r = 0; r < round / batchSize; ++r) {
      std::vector<ETH_Addr> addr(batchSize);
      std::vector<ERC20_Balance> balance(batchSize);
      for (int i = 0; i < batchSize; ++i) {
        addr[i].part[0] = initSize + r * batchSize + i;
      }
      omap.FindParBatch(addr.begin(), addr.end(), balance.begin());
    }
    ocall_measure_time(&end);
    uint64_t timediff = end - start;
    printf("oram find time %f us\n", timediff * 1e-3 / round);
  }
}

void testParOMapPerfDeferWriteBack(size_t mapSize = 5e6,
                                   int threadCount = omp_get_max_threads()) {
  size_t initSize = mapSize - 1e5;
  ParOMap<ETH_Addr, ERC20_Balance, uint32_t> omap(mapSize, threadCount / 2);
  std::function<std::pair<ETH_Addr, ERC20_Balance>(uint64_t)> readerFunc =
      [](uint64_t i) {
        std::pair<ETH_Addr, ERC20_Balance> pr;
        pr.first.part[0] = i;
        return pr;
      };

  EM::VirtualVector::VirtualReader<std::pair<ETH_Addr, ERC20_Balance>> reader(
      initSize, readerFunc);
  uint64_t start, end;
  // printf("init omap of size %lu\n", mapSize);
  ocall_measure_time(&start);
  // omap.InitFromReader(reader);
  omap.Init();
  ocall_measure_time(&end);
  uint64_t initTimediff = end - start;
  for (uint32_t batchSize :
       {100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000}) {
    printf("mapSize = %u, threadCount = %d, batchSize = %u\n", mapSize,
           threadCount, batchSize);
    printf("oram init time %f s\n", initTimediff * 1e-9);
    int round = 1e6;
    uint64_t queryTimediff = 0;

    ocall_measure_time(&start);

    for (size_t r = 0; r < round / batchSize; ++r) {
      uint64_t queryStart, queryEnd;
      ocall_measure_time(&queryStart);
      std::vector<ETH_Addr> addr(batchSize);
      std::vector<ERC20_Balance> balance(batchSize);
      for (int i = 0; i < batchSize; ++i) {
        addr[i].part[0] = initSize + r * batchSize + i;
      }
      omap.FindParBatchDeferWriteBack(addr.begin(), addr.end(),
                                      balance.begin());
      ocall_measure_time(&queryEnd);
      queryTimediff += queryEnd - queryStart;
      omap.ParWriteBack();
    }
    ocall_measure_time(&end);
    uint64_t timediff = end - start;
    printf("oram find time %f us\n", queryTimediff * 1e-3 / round);
    printf("oram find and evict time %f us\n", timediff * 1e-3 / round);
  }
}

void testParOMapPerfSignal(size_t mapSize = 5e6,
                           int threadCount = omp_get_max_threads()) {
  size_t initSize = mapSize;

  ParOMap<uint64_t, Bytes<240>, uint32_t> omap(mapSize, threadCount / 2);
  std::function<std::pair<uint64_t, Bytes<240>>(uint64_t)> readerFunc =
      [](uint64_t i) {
        std::pair<uint64_t, Bytes<240>> pr;
        pr.first = i;
        return pr;
      };

  EM::VirtualVector::VirtualReader<std::pair<uint64_t, Bytes<240>>> reader(
      initSize, readerFunc);
  uint64_t start, end;
  // printf("init omap of size %lu\n", mapSize);
  ocall_measure_time(&start);
  // omap.InitFromReader(reader);
  omap.Init();
  ocall_measure_time(&end);
  uint64_t initTimediff = end - start;
  for (uint32_t batchSize : {100, 200, 500, 1000, 2000, 5000, 10000}) {
    printf("mapSize = %u, threadCount = %d, batchSize = %u\n", mapSize,
           threadCount, batchSize);
    printf("oram init time %f s\n", initTimediff * 1e-9);
    int round = 2e5;
    ocall_measure_time(&start);
    for (size_t r = 0; r < round / batchSize; ++r) {
      std::vector<uint64_t> addr(batchSize);
      std::vector<Bytes<240>> balance(batchSize);
      for (int i = 0; i < batchSize; ++i) {
        addr[i] = initSize + batchSize * (100000UL + r) + i;
      }
      omap.FindParBatch(addr.begin(), addr.end(), balance.begin());
    }
    ocall_measure_time(&end);
    uint64_t timediff = end - start;
    printf("oram find time %f us\n", timediff * 1e-3 / round);
  }
}

void testParOMapPerfDiffCond() {
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
  size_t BackendSize = 1e10;
  EM::Backend::g_DefaultBackend =
      new EM::Backend::MemServerBackend(BackendSize);
  for (uint32_t mapSize :
       {1e5, 2e5, 5e5, 1e6, 2e6, 5e6, 1e7, 2e7, 5e7, 1e8, 2e8, 5e8, 1e9}) {
    for (int threadCount : {2, 4, 8, 16, 32}) {
      try {
        // testParOMapPerf(mapSize, threadCount);
        testParOMapPerfDeferWriteBack(mapSize, threadCount);
        // testParOMapPerfSignal(mapSize, threadCount);
      } catch (const std::runtime_error& e) {
        printf("Caught a runtime_error: %s\n", e.what());
      }
    }
  }
}

void testOHashMapPerfDiffCond() {
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
  size_t BackendSize = 1e11;
  EM::Backend::g_DefaultBackend =
      new EM::Backend::MemServerBackend(BackendSize);
  for (uint32_t mapSize :
       {1e5, 2e5, 5e5, 1e6, 2e6, 5e6, 1e7, 2e7, 5e7, 1e8, 2e8, 5e8, 1e9}) {
    try {
      testOHashMapPerf(mapSize);
      // testOHashMapPerfSignal(mapSize);
    } catch (const std::runtime_error& e) {
      printf("Caught a runtime_error: %s\n", e.what());
    }
  }
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
    // testOmpSpeedup();
    testParOMapPerfDiffCond();
    // testParOMapPerf(5e6, 32);
    // testParOMapPerfDeferWriteBack(5e6, 32);
    // testOHashMapPerf();
    // testOHashMapPerfDiffCond();
    // testRecursiveORAMPerf();
    // testOMapPerf();
    // testOMap();
    // testORAMInit();
    // testORAMReadWrite();
  } catch (std::exception& e) {
    printf("exception: %s\n", e.what());
  }
  return;
}