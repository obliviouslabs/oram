#include "oram/oram.hpp"

#include <gtest/gtest.h>

#include <unordered_map>

#include "external_memory/algorithm/element.hpp"
#include "oram/circuit_oram.hpp"
#include "oram/recursive_oram.hpp"
#include "testutils.hpp"

using namespace ODSL::CircuitORAM;

TEST(ORAM, Init) { ORAM<int> oram(1024); }

TEST(ORAM, CommonSuffixLen) {
  uint64_t a = 3;
  uint64_t b = 5;
  ASSERT_EQ(ODSL::commonSuffixLength(a, b), 1);
}

TEST(ORAM, WithoutPositionMap1) {
  int memSize = 1000;
  ODSL::CircuitORAM::ORAM<uint64_t, 2, 20> oram(memSize);
  std::vector<uint64_t> posMap(memSize);
  std::vector<uint64_t> valMap(memSize);
  for (uint64_t i = 0; i < memSize; i++) {
    uint64_t val = UniformRandom();
    uint64_t pos = oram.Write(i, val);
    posMap[i] = pos;
    valMap[i] = val;
  }
  for (int r = 0; r < 7; ++r) {
    for (uint64_t i = 0; i < memSize; i++) {
      uint64_t val;
      uint64_t pos = oram.Read(posMap[i], i, val);
      posMap[i] = pos;
      ASSERT_EQ(val, valMap[i]);
      // printf("read %lu %lu\n", i, val);
    }
  }
  for (uint64_t i = 0; i < memSize; i++) {
    uint64_t val;
    uint64_t pos = oram.Update(
        posMap[i], i,
        [](uint64_t& x) {
          ++x;
          return true;
        },
        val);
    posMap[i] = pos;
    ++valMap[i];
    ASSERT_EQ(val, valMap[i]);
    // printf("write %lu %lu\n", i, val);
  }
  for (int r = 0; r < 7; ++r) {
    for (uint64_t i = 0; i < memSize; i++) {
      uint64_t val;
      uint64_t pos = oram.Read(posMap[i], i, val);
      posMap[i] = pos;
      ASSERT_EQ(val, valMap[i]);
      // printf("read %lu %lu\n", i, val);
    }
  }
}

TEST(ORAM, BatchUpdate) {
  int memSize = 400;
  ODSL::ORAM<uint64_t> oram(memSize, 1UL << 62);
  std::vector<uint64_t> posMap(memSize);
  std::vector<uint64_t> valMap(memSize);
  for (uint64_t i = 0; i < memSize; i++) {
    uint64_t val = UniformRandom();
    uint64_t pos = oram.Write(i, val);
    posMap[i] = pos;
    valMap[i] = val;
  }
  for (int r = 0; r < 7; ++r) {
    for (uint64_t i = 0; i < memSize; i++) {
      uint64_t val;
      uint64_t pos = oram.Read(posMap[i], i, val);
      posMap[i] = pos;
      ASSERT_EQ(val, valMap[i]);
      // printf("read %lu %lu\n", i, val);
    }
  }
  printf("test batch update with duplicate\n");
  for (uint64_t i = 0; i < memSize; i++) {
    uint64_t peer = UniformRandom(memSize - 1);
    if (!(UniformRandom() % 10)) {
      peer = i;  // test self swap
    }
    auto swapFunc = [&](uint64_t batchSize,
                        uint64_t* pair) -> std::vector<bool> {
      if (pair[0] != valMap[i]) {
        printf("Read wrong result for uid %lu, expected %lu, got %lu\n", i,
               valMap[i], pair[0]);
      }
      if (pair[1] != valMap[peer]) {
        printf("Read wrong result for uid %lu, expected %lu, got %lu\n", peer,
               valMap[peer], pair[1]);
      }
      std::swap(pair[0], pair[1]);
      return {true, true};
    };
    uint64_t posPair[2] = {posMap[i], posMap[peer]};
    uint64_t uidPair[2] = {i, peer};
    const auto& newPoses = oram.BatchUpdate(2, posPair, uidPair, swapFunc);
    posMap[i] = newPoses[0];
    posMap[peer] = newPoses[1];
    std::swap(valMap[i], valMap[peer]);
    // printf("write %lu %lu\n", i, val);
  }
  for (int r = 0; r < 7; ++r) {
    for (uint64_t i = 0; i < memSize; i++) {
      uint64_t val;
      uint64_t pos = oram.Read(posMap[i], i, val);
      posMap[i] = pos;
      ASSERT_EQ(val, valMap[i]);
      // printf("read %lu %lu\n", i, val);
    }
  }
}

TEST(ORAM, BatchUpdateLarge) {
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
  size_t BackendSize = 1e9;
  EM::Backend::g_DefaultBackend =
      new EM::Backend::MemServerBackend(BackendSize);
  int memSize = 2e3;
  ODSL::CircuitORAM::ORAM<uint64_t> oram(memSize, 1UL << 62);
  std::vector<uint64_t> posMap(memSize);
  std::vector<uint64_t> valMap(memSize);
  for (uint64_t i = 0; i < memSize; i++) {
    uint64_t val = UniformRandom();
    uint64_t pos = oram.Write(i, val);
    posMap[i] = pos;
    valMap[i] = val;
  }
  for (uint64_t i = 0; i < memSize; i++) {
    for (int j = 0; j < 100; ++j) {
      uint64_t val;
      uint64_t uid = UniformRandom(memSize - 1);
      uint64_t pos = oram.Read(posMap[uid], uid, val);
      posMap[uid] = pos;
      ASSERT_EQ(val, valMap[uid]);
    }
    uint64_t batchSize = 1000;
    std::vector<uint64_t> batchUid(batchSize);
    std::vector<uint64_t> batchPos(batchSize);
    for (uint64_t j = 0; j < batchSize; j++) {
      batchUid[j] = UniformRandom(memSize - 1);
    }
    std::sort(batchUid.begin(), batchUid.end());
    for (uint64_t j = 0; j < batchSize; j++) {
      batchPos[j] = posMap[batchUid[j]];
    }
    auto updateFunc = [&](uint64_t batchSize,
                          uint64_t* vals) -> std::vector<bool> {
      for (uint64_t j = 0; j < batchSize; j++) {
        if (vals[j] != valMap[batchUid[j]]) {
          printf("Read wrong result for uid %lu, expected %lu, got %lu\n",
                 batchUid[j], valMap[batchUid[j]], vals[j]);
        }
        vals[j]++;
      }
      return std::vector<bool>(batchSize, true);
    };
    const auto& newPoses =
        oram.BatchUpdate(batchSize, &batchPos[0], &batchUid[0], updateFunc);
    for (uint64_t j = 0; j < batchSize; j++) {
      if (j == 0 || batchUid[j] != batchUid[j - 1]) {
        posMap[batchUid[j]] = newPoses[j];
      }
    }
    for (uint64_t j = 0; j < batchSize; j++) {
      if (j == 0 || batchUid[j] != batchUid[j - 1]) {
        valMap[batchUid[j]]++;
      }
    }
  }

  for (uint64_t i = 0; i < memSize; i++) {
    uint64_t val;
    uint64_t pos = oram.Read(posMap[i], i, val);
    posMap[i] = pos;
    ASSERT_EQ(val, valMap[i]);
    // printf("read %lu %lu\n", i, val);
  }
}

TEST(ORAM, WithoutPositionMapMixed) {
  for (int memSize : {2, 3, 5, 7, 9, 33, 40, 55, 127, 129, 543, 678, 1023, 1025,
                      2000, 71429}) {
    ODSL::CircuitORAM::ORAM<uint64_t> oram(memSize);
    std::vector<uint64_t> posMap(memSize, -1);
    std::vector<uint64_t> valMap(memSize);
    int opCount = 1e5;
    for (int i = 0; i < opCount; ++i) {
      int rd = UniformRandom() % 3;
      uint64_t uid = UniformRandom() % memSize;
      switch (rd) {
        case 0: {  // write
          if (posMap[uid] != -1) {
            --i;
            continue;
          }
          uint64_t val = UniformRandom();
          uint64_t pos = oram.Write(uid, val);
          posMap[uid] = pos;
          valMap[uid] = val;
        } break;
        case 1: {  // read
          if (posMap[uid] == -1) {
            --i;
            continue;
          }
          uint64_t val;
          uint64_t pos = oram.Read(posMap[uid], uid, val);
          posMap[uid] = pos;
          ASSERT_EQ(val, valMap[uid]);
        } break;
        case 2: {  // update
          if (posMap[uid] == -1) {
            --i;
            continue;
          }
          uint64_t val;
          uint64_t pos = oram.Update(
              posMap[uid], uid,
              [](uint64_t& x) {
                ++x;
                return true;
              },
              val);
          posMap[uid] = pos;
          ++valMap[uid];
          ASSERT_EQ(val, valMap[uid]);
        } break;

        default:
          break;
      }
    }
  }
}

TEST(ORAM, WithoutPositionMapMixedOverflowHandling) {
  for (int memSize : {2, 3, 5, 7, 9, 33, 40, 55, 127, 129, 543, 678, 1023, 1025,
                      2000, 71429}) {
    ODSL::CircuitORAM::ORAM<uint64_t, 2, 5, uint64_t, uint64_t, 4096, 1> oram(
        memSize);
    std::vector<uint64_t> posMap(memSize, -1);
    std::vector<uint64_t> valMap(memSize);
    int opCount = 1e5;
    for (int i = 0; i < opCount; ++i) {
      int rd = UniformRandom() % 3;
      uint64_t uid = UniformRandom() % memSize;
      switch (rd) {
        case 0: {  // write
          if (posMap[uid] != -1) {
            --i;
            continue;
          }
          uint64_t val = UniformRandom();
          uint64_t pos = oram.Write(uid, val);
          posMap[uid] = pos;
          valMap[uid] = val;
        } break;
        case 1: {  // read
          if (posMap[uid] == -1) {
            --i;
            continue;
          }
          uint64_t val;
          uint64_t pos = oram.Read(posMap[uid], uid, val);
          posMap[uid] = pos;
          ASSERT_EQ(val, valMap[uid]);
        } break;
        case 2: {  // update
          if (posMap[uid] == -1) {
            --i;
            continue;
          }
          uint64_t val;
          uint64_t pos = oram.Update(
              posMap[uid], uid,
              [](uint64_t& x) {
                ++x;
                return true;
              },
              val);
          posMap[uid] = pos;
          ++valMap[uid];
          ASSERT_EQ(val, valMap[uid]);
        } break;

        default:
          break;
      }
    }
  }
}

TEST(ORAM, StashLoad) {
  GTEST_SKIP();
  size_t memSize = 1UL << 16;
  static constexpr int Z = 2;
  static constexpr int stashSize = 50;
  ODSL::CircuitORAM::ORAM<int, Z, stashSize, uint64_t, uint64_t> oram(memSize);
  size_t warmupWindowCount = 1e5;
  size_t windowCount = 1e8;
  size_t windowSize = 10;
  std::vector<uint64_t> posMap(memSize);
  for (int i = 0; i < memSize; ++i) {
    posMap[i] = oram.Write(i, 0);
  }
  std::vector<uint64_t> elementDistribute(60);
  for (int i = 0; i < warmupWindowCount + windowCount; ++i) {
    int windowMaxStashLoad = 0;
    for (int j = 0; j < windowSize; ++j) {
      int val;
      uint64_t idx = UniformRandom(memSize - 1);
      uint64_t oldpos = posMap[idx];
      uint64_t pos = oram.Read(oldpos, idx, val);
      posMap[idx] = pos;

      int stashLoad = 0;
      for (int j = 0; j < stashSize; ++j) {
        if (!oram.getStash()->blocks[j].isDummy()) {
          ++stashLoad;
        }
      }
      windowMaxStashLoad = std::max(windowMaxStashLoad, stashLoad);
    }
    if (i >= warmupWindowCount) {
      ++elementDistribute[windowMaxStashLoad];
    }
  }
  for (int i = 0; i < 60; ++i) {
    printf("%d %lu\n", i, elementDistribute[i]);
  }
}

TEST(ORAM, NonPowerOfTwo) {
  for (int memSize :
       {2, 3, 5, 7, 9, 33, 40, 55, 127, 129, 543, 678, 1023, 1025, 2000}) {
    ORAM<uint64_t> oram(memSize);
    std::vector<uint64_t> posMap(memSize);
    std::vector<uint64_t> valMap(memSize);
    for (uint64_t i = 0; i < memSize; i++) {
      uint64_t val = UniformRandom();
      uint64_t pos = oram.Write(i, val);
      posMap[i] = pos;
      valMap[i] = val;
    }
    for (int r = 0; r < 7; ++r) {
      for (uint64_t i = 0; i < memSize; i++) {
        uint64_t val;
        uint64_t pos = oram.Read(posMap[i], i, val);
        posMap[i] = pos;
        ASSERT_EQ(val, valMap[i]);
      }
    }
  }
}

TEST(ORAM, DiffCacheSize) {
  for (int memSize : {2, 3, 5, 7, 9, 33, 40, 55, 127, 129, 543, 678, 1023, 1025,
                      2000, 3000}) {
    for (size_t cacheSize = 1UL << 17; cacheSize > 1; cacheSize *= 0.95) {
      ODSL::CircuitORAM::ORAM<uint64_t, 2, 20, uint64_t, uint64_t, 288> oram;
      try {
        oram.SetSize(memSize, cacheSize);
      } catch (const std::runtime_error& e) {
        if (strcmp(e.what(), "Circuit ORAM cache size too small") == 0) {
          break;
        } else {
          throw;
        }
      }

      std::vector<uint64_t> posMap(memSize);
      std::vector<uint64_t> valMap(memSize);
      for (uint64_t i = 0; i < memSize; i++) {
        uint64_t val = UniformRandom();
        uint64_t pos = oram.Write(i, val);
        posMap[i] = pos;
        valMap[i] = val;
      }
      for (int r = 0; r < 7; ++r) {
        for (uint64_t i = 0; i < memSize; i++) {
          uint64_t val = 0;
          uint64_t pos = oram.Read(posMap[i], i, val);
          posMap[i] = pos;
          ASSERT_EQ(val, valMap[i]);
        }
      }
    }
  }
}

TEST(ORAM, WithoutPositionMapLargePerf) {
  int memSize = 1 << 20;
  ODSL::CircuitORAM::ORAM<TestElement> oram(memSize);

  uint64_t numAccesses = 1e6;
  auto start = std::chrono::system_clock::now();
  for (uint64_t i = 0; i < numAccesses; i++) {
    TestElement val;
    uint64_t pos = oram.Read(i % memSize, 0, val);
  }
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> diff = end - start;
  printf("Time per access: %f us\n", diff.count() * 1e6 / numAccesses);
}

TEST(ORAM, testInitWithReader) {
  for (int memSize :
       {2,   3,    5,    7,    9,    33,    40,    55,    127,    129,   543,
        678, 1023, 1025, 2000, 4567, 12345, 54321, 71429, 100000, 200000}) {
    uint64_t size = memSize * 1111112UL / 2000000UL;
    if (EM::Backend::g_DefaultBackend) {
      delete EM::Backend::g_DefaultBackend;
    }
    size_t BackendSize = 1e9;
    EM::Backend::g_DefaultBackend =
        new EM::Backend::MemServerBackend(BackendSize);
    ODSL::CircuitORAM::ORAM<TestElement, 2, 50, uint64_t, uint64_t> oram(
        memSize);
    std::vector<uint64_t> valMap(size);
    StdVector<TestElement> vec(size);
    for (int i = 0; i < size; ++i) {
      valMap[i] = UniformRandom();
      vec[i].key = valMap[i];
    }
    using ODSL::UidBlock;
    StdVector<UidBlock<uint64_t>> posMap(size);
    StdVector<TestElement>::Reader reader(vec.begin(), vec.end());
    StdVector<UidBlock<uint64_t>>::Writer posMapWriter(posMap.begin(),
                                                       posMap.end());
    auto start = std::chrono::system_clock::now();
    oram.InitFromReader(reader, posMapWriter);
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start;
    printf("Time for initializing %lu elements of %lu bytes: %f s\n", size,
           sizeof(TestElement), diff.count());
    for (int r = 0; r < 2; ++r) {
      for (uint64_t i = 0; i < size; i++) {
        TestElement val;
        uint64_t pos = oram.Read(posMap[i].data, i, val);
        posMap[i].data = pos;
        ASSERT_EQ(val.key, valMap[i]);
      }
    }
  }
}

TEST(RecursiveORAM, testInitDefault) {
  uint64_t size = 12345;
  ODSL::RecursiveORAM<TestElement> oram(size);
  oram.InitDefault(TestElement());
}

TEST(RecursiveORAM, testReadAfterWrite) {
  uint64_t size = 1234;
  ODSL::RecursiveORAM<int64_t> oram(size);
  oram.InitDefault(0);
  oram.Write(0, 1235);
  int64_t val;
  oram.Read(0, val);
  ASSERT_EQ(val, 1235);
  for (uint64_t i = 0; i < size; i++) {
    oram.Write(i, i * 3);
  }
  for (uint64_t i = 0; i < size; i++) {
    int64_t val;
    oram.Read(i, val);
    ASSERT_EQ(val, i * 3);
  }
}

TEST(RecursiveORAM, testReadAfterInit) {
  uint64_t size = 1234;
  StdVector<int64_t> ref(size);
  ODSL::RecursiveORAM<int64_t> oram(size);
  for (uint64_t i = 0; i < size; i++) {
    ref[i] = UniformRandom();
  }
  StdVector<int64_t>::Reader reader(ref.begin(), ref.end());
  oram.InitFromReaderInPlace(reader);
  for (int round = 0; round < 1e5; ++round) {
    uint64_t addr = UniformRandom(size - 1);
    int64_t val;
    oram.Read(addr, val);
    ASSERT_EQ(val, ref[addr]);
  }
}

TEST(RecursiveORAM, testMixed) {
  uint64_t size = 4312;
  using Vec = EM::CacheFrontVector::Vector<int64_t>;
  Vec ref(size, 0);
  ODSL::RecursiveORAM<int64_t> oram(size);
  for (uint64_t i = 0; i < size; i++) {
    ref[i] = UniformRandom();
  }
  oram.InitFromVector(ref);
  for (int round = 0; round < 1e6; ++round) {
    int op = UniformRandom(2);
    uint64_t addr = UniformRandom(size - 1);
    switch (op) {
      case 0: {
        int64_t val;
        oram.Read(addr, val);
        ASSERT_EQ(val, ref[addr]);
      } break;
      case 1: {
        int64_t val = UniformRandom();
        oram.Write(addr, val);
        ref[addr] = val;
      } break;
      case 2: {
        auto incFunc = [](int64_t& x) {
          ++x;
          return true;
        };
        oram.Access(addr, incFunc);
        ++ref[addr];
      } break;
    }
  }
}

TEST(RecursiveORAM, testBatchAccessDefer) {
  uint64_t size = 1234;
  StdVector<uint64_t> ref(size);
  ODSL::RecursiveORAM<uint64_t, uint64_t> oram(size, 1UL << 62);
  for (uint64_t i = 0; i < size; i++) {
    ref[i] = UniformRandom();
  }
  StdVector<uint64_t>::Reader reader(ref.begin(), ref.end());
  oram.InitFromReaderInPlace(reader);
  for (int round = 0; round < 1e4; ++round) {
    // uint64_t addr = UniformRandom(size - 1);
    // int64_t val;
    // oram.Read(addr, val);
    // ASSERT_EQ(val, ref[addr]);
    std::vector<uint64_t> addrs(100);
    for (uint64_t& addr : addrs) {
      addr = UniformRandom(size - 1);
    }
    std::sort(addrs.begin(), addrs.end());
    auto incFunc = [&](std::vector<uint64_t>& xs) {
      for (size_t i = 0; i < addrs.size(); ++i) {
        if (xs[i] != ref[addrs[i]]) {
          printf("Read wrong result for uid %lu, expected %lu, got %lu\n",
                 addrs[i], ref[addrs[i]], xs[i]);
          abort();
        }
        // else {
        //   printf("Read correct result for uid %lu, got %lu\n", addrs[i],
        //          ref[addrs[i]]);
        // }
        ++xs[i];
      }
    };

    oram.BatchAccessDeferWriteBack(addrs, incFunc);
    for (size_t i = 0; i < addrs.size(); ++i) {
      if (i == 0 || addrs[i] != addrs[i - 1]) {
        ++ref[addrs[i]];
      }
    }
    oram.WriteBack();
    for (uint64_t i = 0; i < addrs.size(); i++) {
      uint64_t val;
      oram.Read(addrs[i], val);
      ASSERT_EQ(val, ref[addrs[i]]);
    }
  }
}