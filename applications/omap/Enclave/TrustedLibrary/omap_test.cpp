#include "../Enclave.h"
#include "Enclave_t.h"
////
#include <omp.h>

#include <functional>
#include <unordered_map>

#include "odsl/omap.hpp"
#include "odsl/omap_short_kv.hpp"
#include "odsl/page_oram.hpp"
#include "odsl/par_omap.hpp"
#include "odsl/recursive_oram.hpp"
#include "sgx_thread.h"
#include "sgx_trts.h"

#define MB << 20

#define ASSERT_TRUE(expr)                           \
  if (!expr) {                                      \
    printf("assert failed at line %d\n", __LINE__); \
    abort();                                        \
  }

using namespace ODSL;
EM::Backend::MemServerBackend* EM::Backend::g_DefaultBackend = nullptr;

#define ASSERT_EQ(a, b)                             \
  if ((a) != (b)) {                                 \
    printf("assert failed at line %d\n", __LINE__); \
    printf("%lu != %lu\n", (a), (b));               \
    abort();                                        \
  }

void testOmpSpeedup() {
  int maxThread = omp_get_max_threads();
  uint64_t start, end;
  int64_t timediffOneThread, timediffMultipleThread;
  printf("bitonic sort benchmark\n");
  ocall_measure_time(&start);
  for (int i = 0; i < maxThread; ++i) {
    std::vector<uint64_t> vec(65546);
    for (int r = 0; r < 100; ++r) {
      Algorithm::BitonicSort(vec);
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
      Algorithm::BitonicSort(vec);
    }
  }
  ocall_measure_time(&end);
  timediffMultipleThread = end - start;
  printf("%d thread %f s\n", omp_get_max_threads(),
         (double)timediffMultipleThread * 1e-9);
  printf("speedup for bitonic sort benchmark (memory intensive task) = %f\n",
         (double)timediffOneThread / (double)timediffMultipleThread);

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
         (double)timediffOneThread / (double)timediffMultipleThread);

  printf("\nFloating point benchmark\n");
  ocall_measure_time(&start);
  double totalSum;
  for (int i = 0; i < maxThread; ++i) {
    double sum = 1.0 + (double)(UniformRandom() % 2);
    for (uint64_t i = 0; i < (uint64_t)1e8; ++i) {
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
    double sum = 1.0 + (double)(UniformRandom() % 2);
    for (uint64_t i = 0; i < (uint64_t)1e8; ++i) {
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
         (double)timediffOneThread / (double)timediffMultipleThread);

  printf("\nHash benchmark\n");
  ocall_measure_time(&start);
  struct Key {
    uint8_t key[20];
  };
  Key k;
  uint8_t salt[16];
  for (int i = 0; i < maxThread; ++i) {
    for (int i = 0; i < (int)1e5; ++i) {
      secure_hash_with_salt(k, salt);
    }
  }
  ocall_measure_time(&end);
  timediffOneThread = end - start;
  printf("one thread %f s\n", (double)timediffOneThread * 1e-9);

  ocall_measure_time(&start);
#pragma omp parallel for schedule(static)
  for (int i = 0; i < maxThread; ++i) {
    for (int i = 0; i < (int)1e5; ++i) {
      secure_hash_with_salt(k, salt);
    }
  }
  ocall_measure_time(&end);
  timediffMultipleThread = end - start;
  printf("%d thread %f s\n", omp_get_max_threads(),
         (double)timediffMultipleThread * 1e-9);
  printf("speedup for hash benchmark (cpu intensive task) = %f\n",
         (double)timediffOneThread / (double)timediffMultipleThread);

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
         (double)timediffOneThread / (double)timediffMultipleThread);

  printf("\nHeapTree benchmark\n");
  ocall_measure_time(&start);

  for (int i = 0; i < maxThread; ++i) {
    HeapTree<uint64_t> tree;
    uint64_t size = 2048;
    tree.Init(size);
    std::vector<uint64_t> path(tree.GetDepth());
    for (int r = 0; r < 1e7; ++r) {
      uint64_t nodeIdxArr[64];
      int depth = tree.GetNodeIdxArr(nodeIdxArr, r % size);
      for (int i = 0; i < depth; ++i) {
        ++tree.GetNodeByIdx(nodeIdxArr[i]);
      }
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
    for (int r = 0; r < 1e7; ++r) {
      uint64_t nodeIdxArr[64];
      int depth = tree.GetNodeIdxArr(nodeIdxArr, r % size);
      for (int i = 0; i < depth; ++i) {
        ++tree.GetNodeByIdx(nodeIdxArr[i]);
      }
    }
  }
  ocall_measure_time(&end);
  timediffMultipleThread = end - start;
  printf("%d thread %f s\n", omp_get_max_threads(),
         (double)timediffMultipleThread * 1e-9);
  printf("one thread %f s\n", (double)timediffOneThread * 1e-9);
  printf("speedup for heap tree benchmark = %f\n",
         (double)timediffOneThread / (double)timediffMultipleThread);

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
         (double)timediffOneThread / (double)timediffMultipleThread);

  printf("\nCircuit ORAM update benchmark\n");
  ocall_measure_time(&start);

  for (int i = 0; i < maxThread; ++i) {
    CircuitORAM::ORAM<TestElement> oram;
    uint64_t size = 65536;
    oram.SetSize(size);
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
    try {
      CircuitORAM::ORAM<TestElement> oram;
      uint64_t size = 8192;
      oram.SetSize(size);
      for (int r = 0; r < 1e5; ++r) {
        oram.Update(UniformRandom() % size, 0, [](TestElement& val) {
          val.key++;
          return true;
        });
      }
    } catch (std::runtime_error& e) {
      printf("%s\n", e.what());
    }
  }
  ocall_measure_time(&end);
  timediffMultipleThread = end - start;
  printf("%d thread %f s\n", omp_get_max_threads(),
         (double)timediffMultipleThread * 1e-9);
  printf("one thread %f s\n", (double)timediffOneThread * 1e-9);
  printf("speedup for circuit oram benchmark = %f\n",
         (double)timediffOneThread / (double)timediffMultipleThread);

  printf("\nLinear ORAM update benchmark\n");
  ocall_measure_time(&start);

  for (int i = 0; i < maxThread; ++i) {
    LinearORAM::ORAM<TestElement> oram;
    uint64_t size = 256;
    oram.SetSize(size);
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
         (double)timediffOneThread / (double)timediffMultipleThread);

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
         (double)timediffOneThread / (double)timediffMultipleThread);
}

template <const uint64_t size>
void testEncrypted() {
  struct TestBlock {
    uint8_t data[size];
  };
  TestBlock b;
  for (uint64_t i = 0; i < size; i++) {
    b.data[i] = 7 * i + 3;
  }
  // check decryption correctness
  Encrypted<TestBlock> e;
  e.Encrypt(b);
  TestBlock b2;
  e.Decrypt(b2);
  for (size_t j = 0; j < size; j++) {
    ASSERT_EQ(b.data[j], b2.data[j]);
  }
  size_t i = UniformRandom(size - 1);
  // Check changing the encrypted data gets different decrypted data
  e.data[i]++;
  TestBlock b3;
  e.Decrypt(b3);
  bool isDifferent = false;
  for (size_t j = 0; j < size; j++) {
    if (b.data[j] != b3.data[j]) {
      isDifferent = true;
      break;
    }
  }
  ASSERT_TRUE(isDifferent);
  e.data[i]--;
  // Check encryption is not deterministic
  Encrypted<TestBlock> e2;
  e2.Encrypt(b);
  bool isDifferent2 = false;
  for (size_t j = 0; j < size; j++) {
    if (e.data[j] != e2.data[j]) {
      isDifferent2 = true;
      break;
    }
  }
  ASSERT_TRUE(isDifferent2);

  // Check FreshEncrypted
  FreshEncrypted<TestBlock> fe;
  uint8_t iv[IV_SIZE];
  GetRandIV(iv);
  fe.Encrypt(b, iv);
  TestBlock b4;
  fe.Decrypt(b4, iv);
  for (size_t j = 0; j < size; j++) {
    ASSERT_EQ(b.data[j], b4.data[j]);
  }
  FreshEncrypted<TestBlock> fe2;
  fe2.Encrypt(b, iv);
  for (size_t j = 0; j < size; j++) {
    ASSERT_EQ(fe.data[j], fe2.data[j]);
  }
  for (size_t j = 0; j < MAC_SIZE; j++) {
    ASSERT_EQ(fe.tag[j], fe2.tag[j]);
  }
  // Check changing the encrypted data can be detected
  i = UniformRandom(size - 1);
  fe.data[i]++;
  TestBlock b5;
  try {
    fe.Decrypt(b5, iv);
    ASSERT_TRUE(false);
  } catch (std::exception& err) {
    ASSERT_TRUE(true);
  }
  fe.data[i]--;
  // Check changing the tag can be detected
  i = UniformRandom(MAC_SIZE - 1);
  fe.tag[i]++;
  try {
    fe.Decrypt(b5, iv);
    ASSERT_TRUE(false);
  } catch (std::exception& err) {
    ASSERT_TRUE(true);
  }
  fe.tag[i]--;
  // Check changing the iv can be detected
  i = UniformRandom(IV_SIZE - 1);
  iv[i]++;
  try {
    fe.Decrypt(b5, iv);
    ASSERT_TRUE(false);
  } catch (std::exception& err) {
    ASSERT_TRUE(true);
  }
}

template <const uint64_t size>
void testEncryptPerf() {
  uint64_t start, end;
  uint64_t round = 1000000UL;
  ocall_measure_time(&start);
  for (uint64_t r = 0; r < round; ++r) {
    struct TestBlock {
      uint8_t data[size];
    };
    TestBlock b;
    for (uint64_t i = 0; i < size; i++) {
      b.data[i] = 7 * i + r;
    }
    uint8_t iv[IV_SIZE];
    FreshEncrypted<TestBlock> e2;
    e2.Encrypt(b, iv);
    e2.Decrypt(b, iv);
  }
  ocall_measure_time(&end);
  uint64_t timediff = end - start;
  printf("encrypt + decrypt %f us\n", (double)timediff * 1e-3 / (double)round);
}

void testORAMReadWrite() {
  for (uint64_t memSize : {1, 3, 5, 7, 9, 33, 40, 55, 127, 129, 543, 678, 1023,
                           1025, 2000, 10000, 50000}) {
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
        uint64_t val = 0;
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
  for (uint64_t i = 0; i < size; ++i) {
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
  OHashMap<uint64_t, int64_t, FULL_OBLIVIOUS> omap(mapSize);
  std::unordered_map<uint64_t, int64_t> map;
  for (size_t i = 0; i < initSize; i++) {
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
  printf("oram init time %f s\n", (double)timediff * 1e-9);
  size_t round = 1e5;
  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    uint64_t i = UniformRandom(mapSize);
    int64_t val = UniformRandom(mapSize * 3);
    bool res = omap.OInsert(i, val);
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
  printf("oram insert time %f us\n", (double)timediff * 1e-3 / (double)round);

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
  printf("oram find time %f us\n", (double)timediff * 1e-3 / (double)round);
}

using ETH_Addr = Bytes<32>;

using ERC20_Balance = Bytes<104>;

void testOMapBatchAccess(size_t mapSize = 1e5) {
  size_t initSize = mapSize;
  printf("default heap size %lu\n", DEFAULT_HEAP_SIZE);
  using OMap =
      OHashMap<ETH_Addr, ERC20_Balance, FULL_OBLIVIOUS, uint32_t, false>;
  OMap omap((uint32_t)mapSize, 10 MB);
  std::unordered_map<ETH_Addr, ERC20_Balance> map;
  std::vector<ETH_Addr> addrs(initSize);
  for (size_t i = 0; i < initSize; i++) {
    ETH_Addr addr;
    addr.SetRand();
    ERC20_Balance balance;
    balance.SetRand();
    map[addr] = balance;
    addrs[i] = addr;
  }

  auto initializer = omap.NewInitContext();
  for (auto& kv : map) {
    initializer->Insert(kv.first, kv.second);
  }
  initializer->Finalize();
  delete initializer;
  printf("omap init done\n");
  size_t round = 1e5;
  uint32_t maxBatchSize = 100;
  uint64_t start, end;
  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    uint32_t batchSize = UniformRandom(1u, maxBatchSize);
    std::vector<ETH_Addr> batchAddrs(batchSize);
    for (size_t i = 0; i < batchSize; i++) {
      if (UniformRandomBit()) {
        batchAddrs[i] = addrs[UniformRandom(initSize - 1)];
      } else {
        batchAddrs[i].SetRand();
      }
    }
    using ValResult = typename OMap::ValResult;
    std::vector<ValResult> balances(batchSize);
    omap.FindBatchDeferWriteBack(batchAddrs.begin(), batchAddrs.end(),
                                 balances.begin());
    for (size_t i = 0; i < batchSize; i++) {
      auto it = map.find(batchAddrs[i]);
      if (it != map.end()) {
        if (balances[i].value != it->second) {
          printf("find failed at round %lu, value not match\n", r);
          abort();
        }
      } else {
        if (balances[i].found) {
          printf("find failed at round %lu, found element that doesn't exist\n",
                 r);
          abort();
        }
      }
    }
    omap.WriteBackTable(0);
    omap.WriteBackTable(1);
  }
  ocall_measure_time(&end);
  uint64_t timediff = end - start;
  printf("oram find batch time %f us\n",
         (double)timediff * 1e-3 / (double)round);
}

void testOHashMapPerf(size_t mapSize = 5e6) {
  size_t round = 1e5;
  size_t initSize = mapSize;
  printf("default heap size %lu\n", DEFAULT_HEAP_SIZE);
  OHashMap<ETH_Addr, ERC20_Balance, FULL_OBLIVIOUS, uint32_t, false> omap(
      (uint32_t)mapSize, DEFAULT_HEAP_SIZE);

  std::function<std::pair<ETH_Addr, ERC20_Balance>(uint64_t)> readerFunc =
      [](uint64_t) { return std::pair<ETH_Addr, ERC20_Balance>(); };

  EM::VirtualVector::VirtualReader<std::pair<ETH_Addr, ERC20_Balance>> reader(
      initSize, readerFunc);
  uint64_t start, end;
  printf("mapSize = %u, threadCount = %d, batchSize = %u\n", mapSize, 1, 1);
  ocall_measure_time(&start);
  omap.InitFromReader(reader);
  ocall_measure_time(&end);
  uint64_t timediff = end - start;
  printf("oram init time %f s\n", (double)timediff * 1e-9);

  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    ETH_Addr addr;
    ERC20_Balance balance;
    omap.Insert(addr, balance);
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram insert time %f us\n", (double)timediff * 1e-3 / (double)round);

  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    ETH_Addr addr;
    addr.SetRand();
    ERC20_Balance balance;
    omap.Find(addr, balance);
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram find time %f us\n", (double)timediff * 1e-3 / (double)round);

  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    ETH_Addr addr = {};
    omap.Erase(addr);
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram erase time %f us\n", (double)timediff * 1e-3 / (double)round);
}

void testOHashMapImproved(size_t mapSize = 5e6) {
  size_t round = 1e5;
  size_t initSize = mapSize;
  printf("default heap size %lu\n", DEFAULT_HEAP_SIZE);
  OMap<ETH_Addr, ERC20_Balance, uint32_t> omap((uint32_t)mapSize,
                                               DEFAULT_HEAP_SIZE / 3 * 2);

  // std::function<std::pair<ETH_Addr, ERC20_Balance>(uint64_t)> readerFunc =
  //     [](uint64_t) { return std::pair<ETH_Addr, ERC20_Balance>(); };

  // EM::VirtualVector::VirtualReader<std::pair<ETH_Addr, ERC20_Balance>>
  // reader(
  //     initSize, readerFunc);
  uint64_t start, end;
  printf("mapSize = %u, threadCount = %d, batchSize = %u\n", mapSize, 1, 1);
  ocall_measure_time(&start);
  omap.Init();
  ocall_measure_time(&end);
  uint64_t timediff = end - start;
  printf("oram init time %f s\n", (double)timediff * 1e-9);
  // std::vector<uint64_t> cacheEvictor(1UL << 24);
  // for (size_t i = 0; i < cacheEvictor.size(); i++) {
  //   cacheEvictor[i] = i;
  // }
  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    ETH_Addr addr;
    ERC20_Balance balance;
    omap.Insert(addr, balance);
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram insert time %f us\n", (double)timediff * 1e-3 / (double)round);
  // for (size_t i = 0; i < cacheEvictor.size(); i++) {
  //   cacheEvictor[i]++;
  // }
  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    ETH_Addr addr;
    addr.SetRand();
    ERC20_Balance balance;
    omap.Find(addr, balance);
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram find time %f us\n", (double)timediff * 1e-3 / (double)round);
}

void testPageOMap(size_t mapSize = 5e6) {
  size_t round = 1e5;
  size_t initSize = mapSize;
  printf("default heap size %lu\n", DEFAULT_HEAP_SIZE);
  OHashMap<ETH_Addr, ERC20_Balance, PAGE_OBLIVIOUS, uint32_t> omap(
      (uint32_t)mapSize, DEFAULT_HEAP_SIZE / 3 * 2);

  // std::function<std::pair<ETH_Addr, ERC20_Balance>(uint64_t)> readerFunc =
  //     [](uint64_t) { return std::pair<ETH_Addr, ERC20_Balance>(); };

  // EM::VirtualVector::VirtualReader<std::pair<ETH_Addr, ERC20_Balance>>
  // reader(
  //     initSize, readerFunc);
  uint64_t start, end;
  printf("mapSize = %u, threadCount = %d, batchSize = %u\n", mapSize, 1, 1);
  ocall_measure_time(&start);
  omap.Init();
  ocall_measure_time(&end);
  uint64_t timediff = end - start;
  printf("oram init time %f s\n", (double)timediff * 1e-9);

  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    ETH_Addr addr;
    ERC20_Balance balance;
    omap.Insert(addr, balance);
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram insert time %f us\n", (double)timediff * 1e-3 / (double)round);

  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    ETH_Addr addr;
    addr.SetRand();
    ERC20_Balance balance;
    omap.Find(addr, balance);
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram find time %f us\n", (double)timediff * 1e-3 / (double)round);
}

void testOHashMapPerfSignal(size_t mapSize = 5e6) {
  size_t round = 1e5;
  size_t initSize = mapSize - round;
  OHashMap<uint64_t, Bytes<240>, FULL_OBLIVIOUS, uint32_t> omap(
      (uint32_t)mapSize);

  std::function<std::pair<uint64_t, Bytes<240>>(uint64_t)> readerFunc =
      [](uint64_t) { return std::pair<uint64_t, Bytes<240>>(); };

  EM::VirtualVector::VirtualReader<std::pair<uint64_t, Bytes<240>>> reader(
      initSize, readerFunc);
  uint64_t start, end;
  printf("mapSize = %u, threadCount = %d, batchSize = %u\n", mapSize, 1, 1);
  ocall_measure_time(&start);
  omap.InitFromReader(reader);
  ocall_measure_time(&end);
  uint64_t timediff = end - start;
  printf("oram init time %f s\n", (double)timediff * 1e-9);

  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    uint64_t addr;
    addr = r;
    Bytes<240> balance;
    omap.Insert(addr, balance);
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram insert time %f us\n", (double)timediff * 1e-3 / (double)round);

  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    uint64_t addr;
    addr = r;
    Bytes<240> balance;
    omap.Find(addr, balance);
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram find time %f us\n", (double)timediff * 1e-3 / (double)round);

  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    uint64_t addr;
    addr = r;
    omap.Erase(addr);
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("oram erase time %f us\n", (double)timediff * 1e-3 / (double)round);
}

void testRecursiveORAMPerf() {
  printf("test recursive oram perf with %d threads\n", TCS_NUM);
  printf("actual working thread max %d\n", omp_get_max_threads());
  size_t mapSize = 5e6;
  struct AddrBalance {
    ETH_Addr addr;
    ERC20_Balance balance;
  };
  RecursiveORAM<AddrBalance, uint32_t> roram((uint32_t)mapSize);

  std::function<AddrBalance(uint64_t)> readerFunc = [](uint64_t) {
    return AddrBalance();
  };

  EM::VirtualVector::VirtualReader<AddrBalance> reader(mapSize, readerFunc);
  uint64_t start, end;
  printf("init recursive oram of size %lu\n", mapSize);
  ocall_measure_time(&start);
  roram.InitFromReader(reader);
  ocall_measure_time(&end);
  uint64_t timediff = end - start;
  printf("oram init time %f s\n", (double)timediff * 1e-9);
  size_t round = 1e5;
  ocall_measure_time(&start);
  for (size_t r = 0; r < round; ++r) {
    AddrBalance pr;
    uint32_t addr = UniformRandom32((uint32_t)mapSize - 1);
    roram.Read(addr, pr);
  }
  ocall_measure_time(&end);
  timediff = end - start;
  printf("recursive oram access time %d.%d us\n", timediff / round / 1'000,
         timediff / round % 1'000);
}

void testParOMapPerf(size_t mapSize = 5e6,
                     int threadCount = omp_get_max_threads()) {
  size_t initSize = mapSize - 100000;
  ParOMap<ETH_Addr, ERC20_Balance, uint32_t> omap(mapSize, threadCount / 2);
  std::function<std::pair<ETH_Addr, ERC20_Balance>(uint64_t)> readerFunc =
      [](uint64_t i) {
        std::pair<ETH_Addr, ERC20_Balance> pr;
        pr.first.SetRand();
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
    printf("oram init time %f s\n", (double)initTimediff * 1e-9);
    size_t round = std::min(500000UL, mapSize);
    ocall_measure_time(&start);
    for (size_t r = 0; r < round / batchSize; ++r) {
      std::vector<ETH_Addr> addr(batchSize);
      std::vector<ERC20_Balance> balance(batchSize);
      omap.InsertBatch(addr.begin(), addr.end(), balance.begin());
    }
    ocall_measure_time(&end);
    uint64_t timediff = end - start;
    printf("oram insert time %f us\n", (double)timediff * 1e-3 / (double)round);
    ocall_measure_time(&start);
    for (size_t r = 0; r < round / batchSize; ++r) {
      std::vector<ETH_Addr> addr(batchSize);
      std::vector<ERC20_Balance> balance(batchSize);
      for (size_t i = 0; i < batchSize; ++i) {
        addr[i].SetRand();
      }
      omap.FindBatch(addr.begin(), addr.end(), balance.begin());
    }
    ocall_measure_time(&end);
    timediff = end - start;
    printf("oram find time %f us\n", (double)timediff * 1e-3 / (double)round);
  }
}

void testParOMapPerfDeferWriteBack(size_t mapSize = 5e6,
                                   int threadCount = omp_get_max_threads()) {
  size_t initSize = mapSize - 100000;
  ParOMap<ETH_Addr, ERC20_Balance, uint32_t> omap(mapSize, threadCount / 2);
  std::function<std::pair<ETH_Addr, ERC20_Balance>(uint64_t)> readerFunc =
      [](uint64_t i) {
        std::pair<ETH_Addr, ERC20_Balance> pr;
        pr.first.SetRand();
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
    printf("oram init time %f s\n", (double)initTimediff * 1e-9);
    int round = 1e6;
    uint64_t queryTimediff = 0;

    ocall_measure_time(&start);

    for (size_t r = 0; r < round / batchSize; ++r) {
      uint64_t queryStart, queryEnd;
      ocall_measure_time(&queryStart);
      std::vector<ETH_Addr> addr(batchSize);
      std::vector<ERC20_Balance> balance(batchSize);
      omap.FindBatchDeferWriteBack(addr.begin(), addr.end(), balance.begin());
      ocall_measure_time(&queryEnd);
      queryTimediff += queryEnd - queryStart;
      omap.WriteBack();
    }
    ocall_measure_time(&end);
    uint64_t timediff = end - start;
    printf("oram find time %f us\n",
           (double)queryTimediff * 1e-3 / (double)round);
    printf("oram find and evict time %f us\n",
           (double)timediff * 1e-3 / (double)round);
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
    printf("oram init time %f s\n", (double)initTimediff * 1e-9);
    int round = 2e5;
    ocall_measure_time(&start);
    for (size_t r = 0; r < round / batchSize; ++r) {
      std::vector<uint64_t> addr(batchSize);
      std::vector<Bytes<240>> balance(batchSize);
      for (uint32_t i = 0; i < batchSize; ++i) {
        addr[i] = initSize + batchSize * (100000UL + r) + i;
      }
      omap.FindBatch(addr.begin(), addr.end(), balance.begin());
    }
    ocall_measure_time(&end);
    uint64_t timediff = end - start;
    printf("oram find time %f us\n", (double)timediff * 1e-3 / (double)round);
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
        testParOMapPerf(mapSize, threadCount);
        // testParOMapPerfDeferWriteBack(mapSize, threadCount);
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
    // testParOMapPerfDiffCond();
    // testParOMapPerf(5e6, 2);
    // testParOMapPerfDeferWriteBack(5e6, 32);
    // testOHashMapPerf();
    // testEncrypted<4096>();
    // testEncryptPerf<4096>();
    // testOHashMapPerfDiffCond();
    // testRecursiveORAMPerf();

    testOHashMapImproved();
    // testPageOMap();
    // printf("heap used %lu\n", g_peak_heap_used);
    // testOMapBatchAccess();
    // testORAMInit();
    // testORAMReadWrite();
  } catch (std::exception& e) {
    printf("exception: %s\n", e.what());
  }
  return;
}