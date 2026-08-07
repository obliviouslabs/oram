#include <gtest/gtest.h>

#include <cstdlib>
#include <random>
#include <unordered_map>

#include <omp.h>

#include "algorithm/element.hpp"
#include "odsl/adaptive_oram.hpp"
#include "odsl/page_oram.hpp"
#include "odsl/recursive_oram.hpp"
#include "testutils.hpp"

using namespace ODSL::CircuitORAM;

TEST(ORAMBasic, CommonSuffixLen) {
  size_t round = 100000UL;
  for (size_t r = 0; r < round; ++r) {
    int desiredLen = UniformRandom(63);
    uint64_t a = UniformRandom();
    uint64_t b = a ^ (((UniformRandom() << 1) | 1) << desiredLen);
    ASSERT_EQ(ODSL::CommonSuffixLength(a, b), desiredLen);
  }
}

TEST(CircuitORAM, SingleAccess) {
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

TEST(CircuitORAM, BatchUpdate) {
  int memSize = 400;
  ODSL::AdaptiveORAM::ORAM<uint64_t> oram(memSize, MAX_CACHE_SIZE);
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

TEST(CircuitORAM, BatchEvictionScheduleIsBalanced) {
  using Schedule =
      ODSL::CircuitORAM::detail::ReverseLexicographicBatchSchedule<uint32_t>;
  // 63, 99, 1001, and 32769 are odd and divisible by the default eviction
  // count of 3, which forces WindowStride to search past its default
  // candidate. That search previously had no bias against even strides.
  for (uint32_t pathCount : {63, 64, 65, 99, 130, 260, 390, 1001, 32769}) {
    uint32_t counter = 0;
    std::vector<uint32_t> evictionCounts(pathCount, 0);

    for (uint32_t i = 0; i < pathCount; ++i) {
      uint32_t windowBegin = Schedule::BeginWindow(counter, pathCount, 3);
      uint32_t firstPath = Schedule::Path(windowBegin, 0, pathCount);
      uint32_t secondPath = Schedule::Path(windowBegin, 1, pathCount);
      ++evictionCounts[firstPath];
      evictionCounts[secondPath] += 2;
      if (pathCount % 2 == 0) {
        ASSERT_NE(firstPath % 2, secondPath % 2);
      }
    }

    for (uint32_t count : evictionCounts) {
      ASSERT_EQ(count, 3);
    }
  }
}

// Regression test for a bug where WindowStride's search for a stride
// coprime with pathCount could land on an even candidate (e.g. pathCount=63
// picks stride=4). Since the root-level tree split is selected by the
// lowest bit of the leaf index (see HeapTree::GetNodeIdxArr), an even
// stride leaves that bit unchanged from one window to the next, pinning
// thousands of consecutive windows to a single root subtree and starving
// the other half of deterministic evictions. An odd stride is required to
// keep the root split alternating every window.
TEST(CircuitORAM, BatchEvictionScheduleStrideIsOdd) {
  using Schedule =
      ODSL::CircuitORAM::detail::ReverseLexicographicBatchSchedule<uint32_t>;
  for (uint32_t evictionCount : {2, 3, 4, 5}) {
    for (uint32_t pathCount = 2; pathCount <= 5000; ++pathCount) {
      uint32_t stride = Schedule::WindowStride(pathCount, evictionCount);
      ASSERT_EQ(stride % 2, 1u)
          << "pathCount=" << pathCount << " evictionCount=" << evictionCount;
    }
    for (uint32_t pathCount : {32769u, 98307u, 131067u}) {
      uint32_t stride = Schedule::WindowStride(pathCount, evictionCount);
      ASSERT_EQ(stride % 2, 1u)
          << "pathCount=" << pathCount << " evictionCount=" << evictionCount;
    }
  }
}

// Directly checks that the root-level split does not stay pinned to one
// side of the tree for long stretches, which is the property that
// BatchEvictionScheduleStrideIsOdd's oddness check is meant to guarantee.
// Bounds the longest run of consecutive windows landing on the same root
// parity over one full coprime cycle.
TEST(CircuitORAM, BatchEvictionScheduleDoesNotStarveRootSubtree) {
  using Schedule =
      ODSL::CircuitORAM::detail::ReverseLexicographicBatchSchedule<uint32_t>;
  for (uint32_t evictionCount : {2, 3, 4, 5}) {
    for (uint32_t pathCount :
        {9u, 27u, 33u, 63u, 99u, 1001u, 32769u, 98307u}) {
      uint32_t counter = 0;
      uint32_t longestRun = 0;
      uint32_t currentRun = 0;
      int lastParity = -1;
      for (uint32_t i = 0; i < pathCount; ++i) {
        uint32_t windowBegin =
            Schedule::BeginWindow(counter, pathCount, evictionCount);
        int parity = static_cast<int>(windowBegin % 2);
        currentRun = (parity == lastParity) ? currentRun + 1 : 1;
        lastParity = parity;
        longestRun = std::max(longestRun, currentRun);
      }
      // An odd stride can only repeat a parity across a modular wrap, never
      // twice in a row otherwise, so runs longer than 2 indicate the root
      // subtree is being starved.
      ASSERT_LE(longestRun, 2u)
          << "pathCount=" << pathCount << " evictionCount=" << evictionCount;
    }
  }
}

template <const int stashSize = 20>
void testBatchUpdateLargeCustomStashSize() {
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
  size_t BackendSize = 1e9;
  EM::Backend::g_DefaultBackend =
      new EM::Backend::MemServerBackend(BackendSize);
  int memSize = UniformRandom(20000, 100000);
  uint64_t round = 500;
  uint64_t maxBatchSize = 1000;
  ODSL::CircuitORAM::ORAM<uint64_t> oram(memSize, 1UL << 20);
  std::vector<uint64_t> posMap(memSize);
  std::vector<uint64_t> valMap(memSize);
  for (uint64_t i = 0; i < memSize; i++) {
    uint64_t val = UniformRandom();
    uint64_t pos = oram.Write(i, val);
    posMap[i] = pos;
    valMap[i] = val;
  }
  for (uint64_t r = 0; r < round; ++r) {
    // add some random reads
    for (int j = 0; j < 10; ++j) {
      uint64_t val;
      uint64_t uid = UniformRandom(memSize - 1);
      uint64_t pos = oram.Read(posMap[uid], uid, val);
      posMap[uid] = pos;
      ASSERT_EQ(val, valMap[uid]);
    }
    uint64_t batchSize = UniformRandom(1UL, maxBatchSize);
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

TEST(CircuitORAM, BatchUpdateLarge) {
  testBatchUpdateLargeCustomStashSize<20>();
}

TEST(CircuitORAM, BatchUpdateLargeOverflows) {
  testBatchUpdateLargeCustomStashSize<3>();
}

TEST(CircuitORAM, BatchReadAndRemovePerf) {
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
  size_t BackendSize = 1e9;
  EM::Backend::g_DefaultBackend =
      new EM::Backend::MemServerBackend(BackendSize);
  int memSize = 1e5;
  uint64_t round = 5000;
  uint64_t batchSize = 1000;
  ODSL::CircuitORAM::ORAM<TestElement> oram(memSize, MAX_CACHE_SIZE);
  std::vector<uint64_t> posMap(memSize, 0);
  std::vector<uint64_t> batchUid(batchSize, 0);
  std::vector<uint64_t> batchPos(batchSize);
  std::vector<TestElement> batchVals(batchSize);
  for (uint64_t r = 0; r < round; ++r) {
    // add some random reads
    for (uint64_t j = 0; j < batchSize; j++) {
      batchPos[j] = posMap[batchUid[j]];
    }
    auto updateFunc = [&](uint64_t batchSize,
                          TestElement* vals) -> std::vector<bool> {
      for (uint64_t j = 0; j < batchSize; j++) {
        vals[j].key++;
      }
      return std::vector<bool>(batchSize, true);
    };
    oram.BatchReadAndRemove(batchSize, &batchPos[0], &batchUid[0],
                            &batchVals[0]);
  }
}

TEST(CircuitORAM, Mixed) {
  for (int memSize : {2, 3, 5, 7, 9, 33, 40, 55, 127, 129, 543, 678, 1023, 1025,
                      2000, 4567, 12345, 25432, 71429}) {
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

TEST(CircuitORAM, OverflowHandling) {
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

// Batch-access counterpart to OverflowHandling above. BatchReadAndRemove /
// BatchWriteBack go through ReverseLexicographicBatchSchedule instead of the
// scalar path's read-path eviction. WindowStride previously could select an
// even stride for pathCounts (size / Z, rounded up) that are odd and a
// multiple of the eviction count, which pins the deterministic schedule to
// one root subtree for thousands of consecutive writebacks and starves the
// other half of the tree (see BatchEvictionScheduleStrideIsOdd and
// BatchEvictionScheduleDoesNotStarveRootSubtree above). RecursiveORAM and
// ParOMap build Circuit ORAMs of exactly these shapes at every recursion
// level / shard, so this checks the production default stash size (33)
// against a mix of small, large, power-of-two-adjacent, and
// odd-multiple-of-three sizes, instead of relying only on the single large
// power-of-two tree that logs/circuit_oram_stash_batched_current.log was
// measured on.
TEST(CircuitORAM, BatchOverflowHandling) {
  using TestORAM = ODSL::CircuitORAM::ORAM<int, 2, 33, uint32_t, uint32_t>;
  for (int size : {2, 3, 5, 7, 9, 33, 40, 55, 127, 129, 543, 678, 1023, 1025,
                   2000, 71429, 126, 198, 522, 2002, 65538}) {
    TestORAM oram(size);
    std::vector<uint32_t> posMap(size);
    std::vector<int> valMap(size, 0);
    for (uint32_t i = 0; i < (uint32_t)size; ++i) {
      posMap[i] = oram.Write(i, 0);
    }
    const std::vector<bool> writeBackFlags(1, true);
    int opCount = 2e5;
    for (int r = 0; r < opCount; ++r) {
      uint32_t idx = UniformRandom(size - 1);
      uint32_t oldPos = posMap[idx];
      uint32_t newPos = oram.GetRandPos();
      int val;
      oram.BatchReadAndRemove(1, &oldPos, &idx, &val);
      ASSERT_EQ(val, valMap[idx]);
      val = ++valMap[idx];
      oram.BatchWriteBack(1, &idx, &newPos, &val, writeBackFlags);
      posMap[idx] = newPos;
    }
  }
}

TEST(CircuitORAM, StashLoad) {
  size_t memSize = 1UL << 16;
  static constexpr int Z = 2;
  static constexpr int stashSize = 50;
  static constexpr int numOrams = 32;
  auto getWindowCount = [](const char* name, size_t defaultValue) {
    const char* value = std::getenv(name);
    return value == nullptr ? defaultValue
                            : static_cast<size_t>(std::strtoull(value, nullptr, 10));
  };
  size_t warmupWindowCount =
      getWindowCount("STASH_LOAD_WARMUP_WINDOWS", 100'000ULL);
  // Calibrated from the long-window runs: about ten hours total when all four
  // variants are run sequentially on the 32-thread test setup.
  size_t windowCount =
      getWindowCount("STASH_LOAD_WINDOWS", 4'450'000'000ULL);
  size_t windowSize = 1;
  double overloadFactor = 1.0;
  size_t elementCount = memSize * overloadFactor;
  std::string variant = "current";
  if (const char* value = std::getenv("CIRCUIT_ORAM_STASH_LOAD_VARIANT")) {
    variant = value;
  }

  auto runStashLoad = [&]<typename ORAMType, bool batchedAccess = false>() {
    std::vector<std::vector<uint64_t>> oramElementDistribute(
        numOrams, std::vector<uint64_t>(stashSize + 1));

#pragma omp parallel for num_threads(numOrams) schedule(static)
    for (int oramIndex = 0; oramIndex < numOrams; ++oramIndex) {
      ORAMType oramForThread(memSize);
      std::mt19937 rng(oramIndex + 1);
      std::uniform_int_distribution<uint32_t> randomPosition(0, memSize / Z - 1);
      std::vector<uint32_t> posMap(elementCount);
      for (uint32_t i = 0; i < elementCount; ++i) {
        posMap[i] = oramForThread.Write(i, 0, randomPosition(rng));
      }
      const std::vector<bool> batchWriteBackFlags(1, true);

      auto& elementDistribute = oramElementDistribute[oramIndex];
      for (size_t i = 0; i < warmupWindowCount + windowCount; ++i) {
        int windowMaxStashLoad = 0;
        for (size_t j = 0; j < windowSize; ++j) {
          int val;
          uint32_t idx = rng() % memSize;
          uint32_t oldpos = posMap[idx];
          uint32_t pos = randomPosition(rng);
          if constexpr (batchedAccess) {
            // BatchReadAndRemove accesses oldpos, but deferred writeback does
            // not revisit it. It evicts deterministic path c once and c+1
            // twice, then advances past the three logical eviction slots.
            uint32_t uid = idx;
            oramForThread.BatchReadAndRemove(1, &oldpos, &uid, &val);
            oramForThread.BatchWriteBack(1, &uid, &pos, &val,
                                         batchWriteBackFlags);
          } else {
            pos = oramForThread.Read(oldpos, idx, val, pos);
          }
          posMap[idx] = pos;

          int stashLoad = 0;
          for (int k = 0; k < stashSize; ++k) {
            if (!oramForThread.GetStash().blocks[k].IsDummy()) {
              ++stashLoad;
            }
          }
          windowMaxStashLoad = std::max(windowMaxStashLoad, stashLoad);
        }
        if (i >= warmupWindowCount) {
          ++elementDistribute[windowMaxStashLoad];
        }
      }
    }

    std::vector<uint64_t> elementDistribute(stashSize + 1);
    for (const auto& oramDistribution : oramElementDistribute) {
      for (size_t load = 0; load < elementDistribute.size(); ++load) {
        elementDistribute[load] += oramDistribution[load];
      }
    }
    for (int i = 0; i <= stashSize; ++i) {
      printf("%d %lu\n", i, elementDistribute[i]);
    }
  };

  using CurrentORAM = ODSL::CircuitORAM::ORAM<
      int, Z, stashSize, uint32_t, uint32_t, 4096, false, 2, 2, true>;
  // No eviction on the accessed/read path, followed by two single evictions
  // on consecutive deterministic paths.
  using NoReadTwoDeterministicPathsORAM = ODSL::CircuitORAM::ORAM<
      int, Z, stashSize, uint32_t, uint32_t, 4096, false, 2, 1, false>;
  using PartialTwoPathORAM = ODSL::CircuitORAM::ORAM<
      int, Z, stashSize, uint32_t, uint32_t, 4096, false, 2, 1, true>;

  if (variant == "current") {
    runStashLoad.template operator()<CurrentORAM>();
  } else if (variant == "batched_current") {
    runStashLoad.template operator()<CurrentORAM, true>();
  } else if (variant == "original" || variant == "no_read_two_paths") {
    runStashLoad.template operator()<NoReadTwoDeterministicPathsORAM>();
  } else if (variant == "partial_two_paths") {
    runStashLoad.template operator()<PartialTwoPathORAM>();
  } else {
    throw std::runtime_error(
        "CIRCUIT_ORAM_STASH_LOAD_VARIANT must be current, batched_current, "
        "original, no_read_two_paths, or partial_two_paths");
  }
}

TEST(CircuitORAM, VariousSizes) {
  for (int memSize : {2, 3, 5, 7, 9, 33, 40, 55, 127, 129, 543, 678, 1023, 1025,
                      2000, 3000, 4567, 12345}) {
    for (size_t cacheSize = 1UL << 19; cacheSize > 1UL << 10;
         cacheSize *= 0.95) {
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
      for (int r = 0; r < 2; ++r) {
        for (uint64_t i = 0; i < memSize; i++) {
          uint64_t uid = UniformRandom(memSize - 1);
          uint64_t val = 0;
          uint64_t pos = oram.Read(posMap[uid], uid, val);
          posMap[uid] = pos;
          ASSERT_EQ(val, valMap[uid]);
        }
      }
    }
  }
}

TEST(CircuitORAM, ReadPerf) {
  int memSize = 1 << 20;

  ODSL::CircuitORAM::ORAM<TestElement, 2, 20, uint32_t, uint32_t, 4096, false>
      oram(memSize);

  uint64_t numAccesses = 3e6;
  auto start = std::chrono::system_clock::now();
  for (uint64_t i = 0; i < numAccesses; i++) {
    TestElement val;
    uint32_t readPos = rand() % memSize;
    oram.Read(readPos, 0, val);
  }
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> diff = end - start;
  printf("Time per access (%lu bytes element): %f us\n", sizeof(TestElement),
         diff.count() * 1e6 / numAccesses);
}

TEST(CircuitORAM, testInitWithReader) {
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

template <bool useCompact>
void testBatchReadWrite() {
  for (uint64_t size = 1; size < 100; ++size) {
    for (uint64_t batchSize = 1; batchSize < 200; ++batchSize) {
      ODSL::LinearORAM::ORAM<int64_t> oram(size);
      std::vector<int64_t> ref(size);
      for (uint64_t i = 0; i < size; i++) {
        oram.Write(i, i * 3);
        ref[i] = i * 3;
      }
      std::vector<uint64_t> addrs(batchSize);
      std::vector<int64_t> vals(batchSize);
      for (int r = 0; r < 10; ++r) {
        for (uint64_t& addr : addrs) {
          addr = UniformRandom(size - 1);
        }
        std::sort(addrs.begin(), addrs.end());
        if (r % 2 == 0) {
          for (ssize_t i = addrs.size() - 1; i >= 0; --i) {
            vals[i] = UniformRandom();
            ref[addrs[i]] = vals[i];
          }
          if constexpr (useCompact) {
            oram.BatchWriteBackViaCompaction(
                batchSize, &addrs[0], &vals[0],
                std::vector<bool>(batchSize, true));
          } else {
            oram.BatchWriteBackNaive(batchSize, &addrs[0], &vals[0],
                                     std::vector<bool>(batchSize, true));
          }
        } else {
          if constexpr (useCompact) {
            oram.BatchReadViaCompaction(batchSize, &addrs[0], &vals[0]);
          } else {
            oram.BatchReadNaive(batchSize, &addrs[0], &vals[0]);
          }
          for (size_t i = 0; i < addrs.size(); ++i) {
            ASSERT_EQ(vals[i], ref[addrs[i]]);
          }
        }
      }
    }
  }
}

TEST(LinearORAM, testBatchReadWriteNaive) { testBatchReadWrite<false>(); }

TEST(LinearORAM, testBatchReadWriteCompact) { testBatchReadWrite<true>(); }

template <const uint64_t element_size>
void testThreshold() {
  using Element = Bytes<element_size>;
  using LinearORAM_ = ODSL::LinearORAM::ORAM<Element, uint32_t>;
  using CircuitORAM_ =
      ODSL::CircuitORAM::ORAM<Element, 2, 20, uint32_t, uint32_t, 4096, false>;
  for (uint32_t size = 10; size < 1000; size += 10) {
    LinearORAM_ linear_oram(size);
    CircuitORAM_ circuit_oram(size);
    // time linear oram performing 1e5 updates
    auto start = std::chrono::system_clock::now();
    for (uint64_t i = 0; i < 1e5; ++i) {
      linear_oram.Update(i % size, [](Element& x) { return false; });
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start;
    printf("Linear  ORAM of size %u (%lu byte element): %f s\n", size,
           element_size, diff.count());
    // time circuit oram performing 1e5 updates
    start = std::chrono::system_clock::now();
    for (uint64_t i = 0; i < 1e5; ++i) {
      circuit_oram.Update(UniformRandom32() % size, i % size,
                          [](Element& x) { return false; });
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff2 = end - start;
    printf("Circuit ORAM of size %u (%lu byte element): %f s\n", size,
           element_size, diff2.count());
    if (diff2.count() < diff.count()) {
      break;
    }
  }
}

TEST(AdaptiveORAM, testThreshold) {
  GTEST_SKIP();
  testThreshold<16>();
  testThreshold<64>();
  testThreshold<128>();
  testThreshold<256>();
}

TEST(RecursiveORAM, testInitDefault) {
  for (uint64_t size = 1234; size < 12345; size = size * 4 / 3) {
    ODSL::RecursiveORAM<int64_t> oram(size);
    oram.InitDefault(54321);
    int64_t val;
    for (uint64_t i = 0; i < size; i += 123) {
      oram.Read(i, val);
      ASSERT_EQ(val, 54321);
    }
  }
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
    oram.Read(i, val);
    ASSERT_EQ(val, i * 3);
  }
}

TEST(RecursiveORAM, testReadAfterWriteWithInitDefault) {
  uint64_t size = 1234;
  ODSL::RecursiveORAM<int64_t> oram(size);
  oram.InitDefault(54321);
  int64_t val;
  for (uint64_t i = 0; i < size / 2; i++) {
    oram.Write(i, i * 3);
  }
  for (uint64_t i = 0; i < size / 2; i++) {
    oram.Read(i, val);
    ASSERT_EQ(val, i * 3);
  }
  for (uint64_t i = size / 2; i < size; i++) {
    oram.Read(i, val);
    ASSERT_EQ(val, 54321);
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
  oram.InitFromReader(reader);
  for (int round = 0; round < 1e5; ++round) {
    uint64_t addr = UniformRandom(size - 1);
    int64_t val;
    oram.Read(addr, val);
    ASSERT_EQ(val, ref[addr]);
  }
}

TEST(RecursiveORAM, testPerf) {
  static constexpr uint64_t blockSize = 64;
  std::vector<uint64_t> sizes = {1UL << 14, 1UL << 16, 1UL << 20};
  for (uint64_t size : sizes) {
    struct Block {
      uint8_t data[blockSize];
    };
    // set backend size to 2GB
    if (EM::Backend::g_DefaultBackend) {
      delete EM::Backend::g_DefaultBackend;
    }
    size_t BackendSize = 10e9;
    EM::Backend::g_DefaultBackend =
        new EM::Backend::MemServerBackend(BackendSize);
    auto init_start = std::chrono::system_clock::now();
    ODSL::RecursiveORAM<Block> oram(size);
    oram.InitDefault(Block());
    auto init_end = std::chrono::system_clock::now();
    std::chrono::duration<double> init_diff = init_end - init_start;
    printf("Time for initializing %lu elements of %lu bytes: %f s\n", size,
           blockSize, init_diff.count());
    uint64_t test_round = 1e5;
    // start timer
    auto start = std::chrono::system_clock::now();
    for (uint64_t i = 0; i < test_round; i++) {
      Block b;
      oram.Write(i % size, b);
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start;
    auto avgTime = diff.count() / test_round;
    printf("Time per write (%lu bytes block): %f us\n", blockSize,
           avgTime * 1e6);
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
  Vec::Reader reader(ref.begin(), ref.end());
  oram.InitFromReader(reader);
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
  for (uint64_t size = 1000; size < 12345; size = size * 3 / 2) {
    StdVector<uint64_t> ref(size);
    ODSL::RecursiveORAM<uint64_t> oram(size, MAX_CACHE_SIZE);
    for (uint64_t i = 0; i < size; i++) {
      ref[i] = UniformRandom();
    }
    StdVector<uint64_t>::Reader reader(ref.begin(), ref.end());
    oram.InitFromReader(reader);
    for (int round = 0; round < 1e2; ++round) {
      std::vector<uint64_t> addrs(rand() % 100 + 1);
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
}

TEST(RecursiveORAM, testBatchAccessDeferInitDefault) {
  for (uint64_t size = 1234; size < 12345; size = size * 3 / 2) {
    StdVector<uint64_t> ref(size);
    ODSL::RecursiveORAM<uint64_t> oram(size, MAX_CACHE_SIZE);
    for (uint64_t i = 0; i < size; i++) {
      ref[i] = 54321;
    }
    oram.InitDefault(54321);
    for (int round = 0; round < 1e2; ++round) {
      std::vector<uint64_t> addrs(rand() % 100 + 1);
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
}

TEST(RecursiveORAM, testBatchAccessDeferLarge) {
  uint32_t size = 654321;
  StdVector<uint64_t> ref(size);
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
  size_t BackendSize = 2e9;
  EM::Backend::g_DefaultBackend =
      new EM::Backend::MemServerBackend(BackendSize);
  ODSL::RecursiveORAM<TestElement> oram(size, 1UL << 24);
  for (uint32_t i = 0; i < size; i++) {
    ref[i] = UniformRandom();
  }
  EM::VirtualVector::VirtualReader<TestElement> reader(size, [&](uint64_t i) {
    TestElement e;
    e.key = ref[i];
    for (int j = 0; j < sizeof(e.payload) / sizeof(e.payload[0]); j++) {
      e.payload[j] = (uint8_t)ref[i];
    }
    return e;
  });
  oram.InitFromReader(reader);
  for (int round = 0; round < 1e3; ++round) {
    std::vector<uint32_t> addrs(UniformRandom(1, 1000));
    for (uint32_t& addr : addrs) {
      addr = UniformRandom32(size - 1);
    }
    std::sort(addrs.begin(), addrs.end());
    auto incFunc = [&](std::vector<TestElement>& xs) {
      for (size_t i = 0; i < addrs.size(); ++i) {
        if (xs[i].key != ref[addrs[i]]) {
          printf("Read wrong result for uid %u, expected %lu, got %lu\n",
                 addrs[i], ref[addrs[i]], xs[i].key);
          abort();
        }

        for (int j = 0; j < sizeof(xs[i].payload) / sizeof(xs[i].payload[0]);
             j++) {
          if (xs[i].payload[j] != (uint8_t)ref[addrs[i]]) {
            printf(
                "Read wrong payload at index %d, ref[addrs[i]] = %lu, "
                "xs[i].key = %lu, "
                "xs[i].payload[j] = %lu\n",
                j, ref[addrs[i]], xs[i].key, (uint64_t)xs[i].payload[j]);
            abort();
          }
          ++xs[i].payload[j];
        }
        ++xs[i].key;
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
      if (rand() % 10 == 0) {
        TestElement val;
        oram.Read(addrs[i], val);
        ASSERT_EQ(val.key, ref[addrs[i]]);
      }
    }
  }
}

TEST(PageORAM, testReadAfterWrite) {
  uint64_t size = 12345;
  ODSL::PageORAM<int64_t> oram(size);
  oram.InitDefault(31);
  oram.Write(0, 1235);
  int64_t val;
  oram.Read(0, val);
  ASSERT_EQ(val, 1235);
  oram.Read(1, val);
  ASSERT_EQ(val, 31);
  for (uint64_t i = 0; i < size; i++) {
    oram.Write(i, i * 3);
  }
  for (uint64_t i = 0; i < size; i++) {
    oram.Read(i, val);
    ASSERT_EQ(val, i * 3);
  }
}

TEST(PageORAM, testReadAfterInit) {
  uint64_t size = 123456;
  StdVector<int64_t> ref(size);
  ODSL::PageORAM<int64_t> oram(size);
  for (uint64_t i = 0; i < size; i++) {
    ref[i] = UniformRandom();
  }
  StdVector<int64_t>::Reader reader(ref.begin(), ref.end());
  oram.InitFromReader(reader, 1UL << 20);
  for (int round = 0; round < 1e5; ++round) {
    uint64_t addr = UniformRandom(size - 1);
    int64_t val;
    oram.Read(addr, val);
    ASSERT_EQ(val, ref[addr]);
  }
}

TEST(PageORAM, testMixed) {
  uint64_t size = 53235;
  using Vec = StdVector<TestElement>;
  Vec ref(size);
  ODSL::PageORAM<TestElement> oram(size);
  for (uint64_t i = 0; i < size; i++) {
    ref[i].key = UniformRandom();
  }
  Vec::Reader reader(ref.begin(), ref.end());
  oram.InitFromReader(reader, 1UL << 17);
  for (int round = 0; round < 1e6; ++round) {
    int op = UniformRandom(2);
    uint64_t addr = UniformRandom(size - 1);
    switch (op) {
      case 0: {
        TestElement val;
        oram.Read(addr, val);
        ASSERT_EQ(val.key, ref[addr].key);
      } break;
      case 1: {
        TestElement val;
        val.key = UniformRandom();
        oram.Write(addr, val);
        ref[addr].key = val.key;
      } break;
      case 2: {
        auto incFunc = [](TestElement& x) {
          ++x.key;
          return true;
        };
        oram.Access(addr, incFunc);
        ++ref[addr].key;
      } break;
    }
  }
  oram.PrintMemoryUsage();
  // std::cout << oram.
}
