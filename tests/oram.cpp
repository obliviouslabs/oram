#include <gtest/gtest.h>

#include <unordered_map>

#include "external_memory/algorithm/sort_def.hpp"
#include "oram/path_oram.hpp"
#include "testutils.hpp"

using namespace ORAM::PathORAM;

TEST(PathORAM, Init) {
  PathORAM<int> oram;
  oram.Init(1024);
  oram.Destroy();
}

TEST(PathORAM, ReadPath) {
  PathORAM<int> oram;
  oram.Init(32);
  std::vector<PathORAM<int>::Block_> path = oram.ReadPath(UniformRandom(31));
}

TEST(PathORAM, ExtractFromPath) {
  PathORAM<int> oram;
  oram.Init(32);
  std::vector<PathORAM<int>::Block_> path = oram.ReadPath(UniformRandom(31));
  int out;
  PathORAM<int>::ReadElementAndRemoveFromPath(path, UniformRandom(31), out);
}

TEST(PathORAM, EvictPath) {
  constexpr int Z = 5;
  constexpr int stashSize = 10;
  int depth = 4;
  int pathSize = stashSize + Z * depth;
  uint64_t maxPos = (1UL << (depth - 1)) - 1;
  using PathORAM_ = PathORAM<int, Z, stashSize>;
  for (int r = 0; r < 100; ++r) {
    std::vector<PathORAM_::Block_> blocks(pathSize);
    uint64_t pos = UniformRandom(maxPos);
    for (int i = 0; i < pathSize; i++) {
      while (true) {
        blocks[i].uid = UniformRandom() % 5 ? 1 : DUMMY<uint64_t>();
        blocks[i].position = UniformRandom(maxPos);
        int suffixLen =
            PathORAM<int>::commonSuffixLength(blocks[i].position, pos);
        int level = (i - stashSize) / Z;
        if (blocks[i].isDummy() || suffixLen == 64 || suffixLen >= level) {
          break;
        }
      }
    }
    PathORAM_::EvictPath(blocks, pos);
    int maxCommonSuffix = -1;
    for (int i = 0; i < pathSize; i++) {
      int suffixLen = blocks[i].isDummy() ? -1
                                          : PathORAM<int>::commonSuffixLength(
                                                blocks[i].position, pos);
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

TEST(PathORAM, EvictPathPerf) {
  constexpr int Z = 5;
  constexpr int stashSize = 10;
  int depth = 4;
  int pathSize = stashSize + Z * depth;
  using PathORAM_ = PathORAM<int, Z, stashSize>;
  for (int r = 0; r < 1000000; ++r) {
    std::vector<PathORAM_::Block_> blocks(pathSize);
    PathORAM_::EvictPath(blocks, 0);
  }
}

TEST(PathORAM, CommonSuffixLen) {
  uint64_t a = 3;
  uint64_t b = 5;
  ASSERT_EQ(PathORAM<int>::commonSuffixLength(a, b), 1);
}

TEST(PathORAM, WithoutPositionMap1) {
  int memSize = 1024;
  PathORAM<uint64_t, 5, 63> oram;
  oram.Init(memSize);
  std::vector<uint64_t> posMap(memSize);
  std::vector<uint64_t> valMap(memSize);
  for (uint64_t i = 0; i < memSize; i++) {
    uint64_t val = UniformRandom();
    uint64_t pos = oram.Write(i, val);
    posMap[i] = pos;
    valMap[i] = val;

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
  for (uint64_t i = 0; i < memSize; i++) {
    uint64_t val;
    uint64_t pos = oram.Update(
        posMap[i], i, [](uint64_t x) { return x + 1; }, val);
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

TEST(PathORAM, WithoutPositionMapMixed) {
  int memSize = 4000;
  PathORAM<uint64_t, 5, 63> oram;
  oram.Init(memSize);
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
            posMap[uid], uid, [](uint64_t x) { return x + 1; }, val);
        posMap[uid] = pos;
        ++valMap[uid];
        ASSERT_EQ(val, valMap[uid]);
      } break;

      default:
        break;
    }
  }
}

TEST(PathORAM, NonPowerOfTwo) {
  for (int memSize :
       {1, 3, 5, 7, 9, 33, 40, 55, 127, 129, 543, 678, 1023, 1025, 2000}) {
    PathORAM<uint64_t, 5, 20> oram;
    oram.Init(memSize);
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

TEST(PathORAM, WithoutPositionMapLargePerf) {
  int memSize = 1 << 24;
  PathORAM<SortElement, 5, 63> oram;
  oram.Init(memSize);

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