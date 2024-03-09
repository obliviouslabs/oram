#include <gtest/gtest.h>

#include <unordered_map>

#include "external_memory/algorithm/sort_def.hpp"
#include "oram/circuit_oram.hpp"
#include "oram/path_oram.hpp"
#include "oram/recursive_oram.hpp"
#include "testutils.hpp"

using namespace ODSL::PathORAM;

TEST(ORAM, Init) { ORAM<int> oram(1024); }

TEST(ORAM, ReadPath) {
  ORAM<int> oram(32);
  std::vector<ORAM<int>::Block_> path = oram.ReadPath(UniformRandom(31));
}

TEST(ORAM, ReadPath2) {
  ORAM<int> oram(1e7);
  std::vector<ORAM<int>::Block_> path = oram.ReadPath(1801922);
}

TEST(ORAM, ExtractFromPath) {
  ORAM<int> oram(32);
  std::vector<ORAM<int>::Block_> path = oram.ReadPath(UniformRandom(31));
  int out;
  ODSL::ReadElementAndRemoveFromPath(path, UniformRandom(31), out);
}

TEST(ORAM, EvictPath) {
  constexpr int Z = 5;
  constexpr int stashSize = 10;
  int depth = 4;
  int pathSize = stashSize + Z * depth;
  uint64_t maxPos = (1UL << (depth - 1)) - 1;
  using ORAM_ = ORAM<int, Z, stashSize>;
  for (int r = 0; r < 100; ++r) {
    std::vector<ORAM_::Block_> blocks(pathSize);
    uint64_t pos = UniformRandom(maxPos);
    for (int i = 0; i < pathSize; i++) {
      while (true) {
        blocks[i].uid = UniformRandom() % 5 ? 1 : DUMMY<uint64_t>();
        blocks[i].position = UniformRandom(maxPos);
        int suffixLen = ODSL::commonSuffixLength(blocks[i].position, pos);
        int level = (i - stashSize) / Z;
        if (blocks[i].isDummy() || suffixLen == 64 || suffixLen >= level) {
          break;
        }
      }
    }
    ORAM_::EvictPath(blocks, pos);
    int maxCommonSuffix = -1;
    for (int i = 0; i < pathSize; i++) {
      int suffixLen = blocks[i].isDummy()
                          ? -1
                          : ODSL::commonSuffixLength(blocks[i].position, pos);
      maxCommonSuffix = std::max(maxCommonSuffix, suffixLen);
      int level = std::max(0, (i - stashSize) / Z);
      if (!blocks[i].isDummy()) {
        // filled slot must have suffixLen > level
        ASSERT_LE(level, suffixLen);
      } else {
        // empty slot must be filled
        ASSERT_LT(maxCommonSuffix, level);
      }
    }
  }
}

TEST(ORAM, EvictPathPerf) {
  constexpr int Z = 5;
  constexpr int stashSize = 10;
  int depth = 4;
  int pathSize = stashSize + Z * depth;
  using ORAM_ = ORAM<int, Z, stashSize>;
  for (int r = 0; r < 1000000; ++r) {
    std::vector<ORAM_::Block_> blocks(pathSize);
    ORAM_::EvictPath(blocks, 0);
  }
}

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
    auto swapFunc = [](std::vector<uint64_t>& pair) -> std::vector<bool> {
      std::swap(pair[0], pair[1]);
      return {true, true};
    };
    uint64_t peer = UniformRandom(memSize - 1);
    if (peer == i) {
      peer = (peer + 1) % memSize;
    }
    const auto& newPoses =
        oram.BatchUpdate({posMap[i], posMap[peer]}, {i, peer}, swapFunc);
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

TEST(ORAM, ParBatchUpdate) {
  int memSize = 6543;
  int numThreads = 4;
  ODSL::CircuitORAM::ParORAM<uint64_t> oram(memSize, 1UL << 62, numThreads);
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
  printf("test par batch update\n");
  for (uint64_t i = 0; i < memSize; i++) {
    uint64_t peer = UniformRandom(memSize - 1);
    if (!(UniformRandom() % 10)) {
      peer = i;  // test self swap
    }
    auto swapFunc = [&](std::vector<uint64_t>& pair) -> std::vector<bool> {
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
    std::vector<uint64_t> posPair = {posMap[i], posMap[peer]};
    const auto& newPoses =
        oram.ParBatchUpdate(posPair, {i, peer}, swapFunc, numThreads);
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

TEST(ORAM, ParBatchUpdateLarge) {
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
  size_t BackendSize = 1e9;
  EM::Backend::g_DefaultBackend =
      new EM::Backend::MemServerBackend(BackendSize);
  int memSize = 1e5;
  int numThreads = 32;
  ODSL::CircuitORAM::ParORAM<uint64_t> oram(memSize, 1UL << 62, numThreads);
  std::vector<uint64_t> posMap(memSize);
  std::vector<uint64_t> valMap(memSize);
  for (uint64_t i = 0; i < memSize; i++) {
    uint64_t val = UniformRandom();
    uint64_t pos = oram.Write(i, val);
    posMap[i] = pos;
    valMap[i] = val;
  }
  for (uint64_t i = 0; i < memSize; i++) {
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
    auto updateFunc = [&](std::vector<uint64_t>& vals) -> std::vector<bool> {
      for (uint64_t j = 0; j < batchSize; j++) {
        if (vals[j] != valMap[batchUid[j]]) {
          printf("Read wrong result for uid %lu, expected %lu, got %lu\n",
                 batchUid[j], valMap[batchUid[j]], vals[j]);
        }
        vals[j] = i * batchUid[j];
      }
      return std::vector<bool>(batchSize, true);
    };
    const auto& newPoses =
        oram.ParBatchUpdate(batchPos, batchUid, updateFunc, numThreads);
    for (uint64_t j = 0; j < batchSize; j++) {
      posMap[batchUid[j]] = newPoses[j];
    }
    for (uint64_t j = 0; j < batchSize; j++) {
      // valMap[batchUid[j]] = valMap[batchUid[(j + 1) % batchSize]] + 1;
      valMap[batchUid[j]] = i * batchUid[j];
    }
    // printf("write %lu %lu\n", i, val);
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
        if (!oram.stash->blocks[j].isDummy()) {
          ++stashLoad;
        }
      }
      for (int j = 0; j < Z; ++j) {
        if (!oram.tree.GetByPathAndLevel(0, 0).blocks[j].isDummy()) {
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

TEST(ORAM, StashLoadOPRAM) {
  size_t memSize = 1UL << 16;
  static constexpr int Z = 2;
  static constexpr int stashSize = 50;
  ODSL::CircuitORAM::ORAM<int, Z, stashSize, uint64_t, uint64_t, 4096, 2> oram(
      memSize);
  size_t warmupWindowCount = 1e5;
  size_t windowCount = 1e8;
  size_t windowSize = 10;
  double dummyProb = 0;
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
      if (UniformRandomDouble() < dummyProb) {
        idx = DUMMY<uint64_t>();
        oldpos = UniformRandom(memSize - 1);
      }
      uint64_t pos = oram.Read(oldpos, idx, val);
      if (idx != DUMMY<uint64_t>()) {
        posMap[idx] = pos;
      }
      int stashLoad = 0;
      for (int j = 0; j < stashSize; ++j) {
        if (!oram.stash->blocks[j].isDummy()) {
          ++stashLoad;
        }
      }
      for (int j = 0; j < Z; ++j) {
        if (!oram.tree.GetByPathAndLevel(0, 0).blocks[j].isDummy()) {
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

// TEST(ORAM, DiffCacheSize) {
//   for (int memSize :
//        {2, 3, 5, 7, 9, 33, 40, 55, 127, 129, 543, 678, 1023, 1025, 2000}) {
//     for (size_t cacheSize = 6910886 / 512; cacheSize > 8192;
//          cacheSize *= 0.95) {
//       ODSL::CircuitORAM::ORAM<uint64_t, 2, 50, uint64_t, uint64_t, 288> oram(
//           memSize, cacheSize);

//       std::vector<uint64_t> posMap(memSize);
//       std::vector<uint64_t> valMap(memSize);
//       for (uint64_t i = 0; i < memSize; i++) {
//         uint64_t val = UniformRandom();
//         uint64_t pos = oram.Write(i, val);
//         posMap[i] = pos;
//         valMap[i] = val;

//         // printf("write %lu %lu to pos %lu\n", i, val, pos);
//       }
//       for (int r = 0; r < 7; ++r) {
//         for (uint64_t i = 0; i < memSize; i++) {
//           uint64_t val = 0;
//           uint64_t pos = oram.Read(posMap[i], i, val);
//           posMap[i] = pos;
//           if (val != valMap[i]) {
//             std::vector<uint64_t> pathIdxs(GetLogBaseTwo(memSize) + 2);
//             int pathLen = oram.tree.GetPathIdx(pathIdxs.begin(),
//             pathIdxs.end(),
//                                                posMap[i], memSize,
//                                                oram.tree.GetCacheLevel());
//             for (int i = 0; i < pathLen; i++) {
//               printf("%lu ", pathIdxs[i]);
//             }
//             printf("\n");
//           }
//           ASSERT_EQ(val, valMap[i]);
//         }
//       }
//     }
//   }
// }

TEST(ORAM, WithoutPositionMapLargePerf) {
  int memSize = 1 << 20;
  ODSL::CircuitORAM::ORAM<SortElement> oram(memSize);

  uint64_t numAccesses = 1e6;
  auto start = std::chrono::system_clock::now();
  for (uint64_t i = 0; i < numAccesses; i++) {
    SortElement val;
    uint64_t pos = oram.Read(i % memSize, 0, val);
  }
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> diff = end - start;
  printf("Time per access: %f us\n", diff.count() * 1e6 / numAccesses);
}

// TEST(ORAM, testInit) {
//   uint64_t size = 1024;
//   ORAM<SortElement> oram(size);
//   StdVector<SortElement> vec(size);
//   for (int i = 0; i < size; ++i) {
//     vec[i] = SortElement();
//     vec[i].key = i;
//   }
//   StdVector<uint64_t> posMap(size);
//   StdVector<uint64_t>::Writer posMapWriter(posMap.begin(), posMap.end());
//   oram.InitFromVector(vec, posMapWriter);
//   for (uint64_t i = 0; i < size; i++) {
//     SortElement val;
//     printf("read %lu at pos %lu\n", i, posMap[i]);
//     uint64_t pos = oram.Read(posMap[i], i, val);
//     posMap[i] = pos;
//     ASSERT_EQ(val.key, i);
//   }
// }

TEST(ORAM, testInitNaive) {
  uint64_t memSize = 123456;
  uint64_t size = 100000;
  ORAM<SortElement> oram(memSize);
  for (uint64_t i = 0; i < size; i++) {
    oram.Write(i, SortElement());
  }
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
    ODSL::CircuitORAM::ORAM<SortElement, 2, 50, uint64_t, uint64_t> oram(
        memSize);
    std::vector<uint64_t> valMap(size);
    StdVector<SortElement> vec(size);
    for (int i = 0; i < size; ++i) {
      valMap[i] = UniformRandom();
      vec[i].key = valMap[i];
    }
    using ODSL::UidBlock;
    StdVector<UidBlock<uint64_t>> posMap(size);
    StdVector<SortElement>::Reader reader(vec.begin(), vec.end());
    StdVector<UidBlock<uint64_t>>::Writer posMapWriter(posMap.begin(),
                                                       posMap.end());
    auto start = std::chrono::system_clock::now();
    oram.InitFromReader(reader, posMapWriter);
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start;
    printf("Time for initializing %lu elements of %lu bytes: %f s\n", size,
           sizeof(SortElement), diff.count());
    // for (auto& p : posMap) {
    //   printf("(%lu, %lu)\n", p.uid, p.data);
    // }
    // for (uint64_t i = 0; i < size * 2 - 1; i++) {
    //   for (int j = 0; j < 5; ++j) {
    //     printf("%ld, ",
    //            (int64_t)oram.tree.GetByInternalIdx(i).blocks[j].position);
    //   }
    //   printf("\n");
    // }
    for (int r = 0; r < 2; ++r) {
      for (uint64_t i = 0; i < size; i++) {
        // printf("read %lu at pos %lu\n", i, posMap[i].data);
        SortElement val;
        // printf("read %lu %lu at pos %lu\n", i, val, posMap[i]);
        uint64_t pos = oram.Read(posMap[i].data, i, val);
        posMap[i].data = pos;
        ASSERT_EQ(val.key, valMap[i]);
      }
    }
  }
}

TEST(RecursiveORAM, testInitDefault) {
  uint64_t size = 12345;
  ODSL::RecursiveORAM<SortElement> oram(size);
  oram.InitDefault(SortElement());
}

TEST(RecuriveORAM, testReadAfterWrite) {
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

TEST(RecuriveORAM, testReadAfterInit) {
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

TEST(RecuriveORAM, testMixed) {
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

TEST(RecuriveORAM, testBatchAccessUnique) {
  uint64_t size = 123456;
  StdVector<uint64_t> ref(size);
  int numThreads = 8;
  ODSL::RecursiveORAM<uint64_t, uint64_t, true> oram(size, 1UL << 62,
                                                    numThreads);
  for (uint64_t i = 0; i < size; i++) {
    ref[i] = UniformRandom();
  }
  StdVector<uint64_t>::Reader reader(ref.begin(), ref.end());
  oram.InitFromReaderInPlace(reader);
  for (int round = 0; round < 1e3; ++round) {
    // uint64_t addr = UniformRandom(size - 1);
    // int64_t val;
    // oram.Read(addr, val);
    // ASSERT_EQ(val, ref[addr]);
    std::vector<uint64_t> addrs(1000);
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

    oram.ParBatchAccess(addrs, incFunc, numThreads);
    for (size_t i = 0; i < addrs.size(); ++i) {
      if (i == 0 || addrs[i] != addrs[i - 1]) {
        ++ref[addrs[i]];
      }
    }
  }
}