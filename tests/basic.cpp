#include <gtest/gtest.h>

#include "common/mov_intrinsics.hpp"
#include "common/probability.hpp"
#include "oram/heap_tree.hpp"
#include "testutils.hpp"

template <const uint64_t size>
struct TestBlock {
  uint8_t data[size];
};

template <const uint64_t size>
void testObliMov() {
  TestBlock<size>* b = new TestBlock<size>();
  TestBlock<size>* b2 = new TestBlock<size>();
  printf("Test size: %lu\n", size);
  for (uint64_t i = 0; i < size; i++) {
    b->data[i] = 7 * i + 3;
    b2->data[i] = 11 * i + 5;
  }
  printf("Fake mov\n");
  obliMove(false, *b, *b2);
  for (size_t i = 0; i < size; i++) {
    ASSERT_EQ(b->data[i], (uint8_t)(7 * i + 3));
    ASSERT_EQ(b2->data[i], (uint8_t)(11 * i + 5));
  }
  printf("True mov\n");
  obliMove(true, *b, *b2);
  for (size_t i = 0; i < size; i++) {
    ASSERT_EQ(b->data[i], (uint8_t)(11 * i + 5));
    ASSERT_EQ(b2->data[i], (uint8_t)(11 * i + 5));
  }
  delete b;
  delete b2;
}

template <const uint64_t size>
void testObliSwap() {
  TestBlock<size>* b = new TestBlock<size>();
  TestBlock<size>* b2 = new TestBlock<size>();
  printf("Test size: %lu\n", size);
  for (uint64_t i = 0; i < size; i++) {
    b->data[i] = 7 * i + 3;
    b2->data[i] = 11 * i + 5;
  }
  printf("Fake swap\n");
  obliSwap(false, *b, *b2);
  for (size_t i = 0; i < size; i++) {
    ASSERT_EQ(b->data[i], (uint8_t)(7 * i + 3));
    ASSERT_EQ(b2->data[i], (uint8_t)(11 * i + 5));
  }
  printf("True swap\n");
  obliSwap(true, *b, *b2);
  for (size_t i = 0; i < size; i++) {
    ASSERT_EQ(b->data[i], (uint8_t)(11 * i + 5));
    ASSERT_EQ(b2->data[i], (uint8_t)(7 * i + 3));
  }
  delete b;
  delete b2;
}

template <const uint64_t size>
void testObliMovPerf() {
  TestBlock<size>* b = new TestBlock<size>();
  TestBlock<size>* b2 = new TestBlock<size>();
  printf("Test size: %lu\n", size);
  auto start = std::chrono::system_clock::now();
  for (uint64_t i = 0; i < size; i++) {
    b->data[i] = 7 * i + 3;
    b2->data[i] = 11 * i + 5;
  }
  for (size_t i = 0; i < 1000000; i++) {
    obliMove(false, *b, *b2);
    obliMove(true, *b, *b2);
  }
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  printf("Elapsed time: %f\n", elapsed_seconds.count());

  delete b;
  delete b2;
}

template <const uint64_t size>
void testObliSwapPerf() {
  TestBlock<size>* b = new TestBlock<size>();
  TestBlock<size>* b2 = new TestBlock<size>();
  printf("Test size: %lu\n", size);
  auto start = std::chrono::system_clock::now();
  for (uint64_t i = 0; i < size; i++) {
    b->data[i] = 7 * i + 3;
    b2->data[i] = 11 * i + 5;
  }
  for (size_t i = 0; i < 1000000; i++) {
    obliSwap(false, *b, *b2);
    obliSwap(true, *b, *b2);
  }
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  printf("Elapsed time: %f\n", elapsed_seconds.count());

  delete b;
  delete b2;
}

TEST(Basic, ObliMove) {
  testObliMov<1>();
  testObliMov<2>();
  testObliMov<3>();
  testObliMov<7>();
  testObliMov<8>();
  testObliMov<9>();
  testObliMov<12>();
  testObliMov<15>();
  testObliMov<16>();
  testObliMov<17>();
  testObliMov<30>();
  testObliMov<31>();
  testObliMov<32>();
  testObliMov<33>();
  testObliMov<63>();
  testObliMov<64>();
  testObliMov<65>();
  testObliMov<127>();
  testObliMov<128>();
  testObliMov<129>();
  testObliMov<200>();
  testObliMov<256>();
  testObliMov<257>();
  testObliMov<300>();
  testObliMov<512>();
  testObliMov<513>();
  testObliMov<600>();
  testObliMov<1023>();
  testObliMov<1024>();
  testObliMov<1025>();
  testObliMov<2000>();
}

TEST(Basic, ObliSwap) {
  testObliSwap<1>();
  testObliSwap<2>();
  testObliSwap<3>();
  testObliSwap<7>();
  testObliSwap<8>();
  testObliSwap<9>();
  testObliSwap<12>();
  testObliSwap<15>();
  testObliSwap<16>();
  testObliSwap<17>();
  testObliSwap<30>();
  testObliSwap<31>();
  testObliSwap<32>();
  testObliSwap<33>();
  testObliSwap<63>();
  testObliSwap<64>();
  testObliSwap<65>();
  testObliSwap<127>();
  testObliSwap<128>();
  testObliSwap<129>();
  testObliSwap<200>();
  testObliSwap<256>();
  testObliSwap<257>();
  testObliSwap<300>();
  testObliSwap<512>();
  testObliSwap<513>();
  testObliSwap<600>();
  testObliSwap<1023>();
  testObliSwap<1024>();
  testObliSwap<1025>();
  testObliSwap<2000>();
}

TEST(Basic, MovPerf) { testObliMovPerf<200>(); }

TEST(Basic, SwapPerf) { testObliSwapPerf<200>(); }

TEST(Basic, TestHeapTree) {
  size_t size = 11;
  int maxLevel = GetLogBaseTwo(size - 1) + 2;
  for (int level = 0; level < maxLevel; ++level) {
    for (int i = 0; i < size; ++i) {
      size_t idx = HeapTree<int>::GetIdx(i, size, level, maxLevel, maxLevel);
      printf("level: %d, i: %d, idx: %lu\n", level, i, idx);
    }
  }
}

TEST(Basic, TestHeapTree2) {
  size_t size = 13;
  int maxLevel = GetLogBaseTwo(size - 1) + 2;
  std::vector<size_t> pathIdxs(maxLevel);
  for (size_t i = 0; i < size; ++i) {
    int level = HeapTree<int, uint64_t, 32>::GetPathIdx(
        pathIdxs.begin(), pathIdxs.end(), i, size, 1);
    printf("Path %lu:\n", i);
    for (int j = 0; j < level; ++j) {
      printf("%lu ", pathIdxs[j]);
    }
    printf("\n");
  }
}

template <const size_t page_size>
void testHeapTreeCorrectness() {
  // printf(
  //     "Testing heap tree correctness for page size: %lu, packing %d
  //     levels\n", page_size, HeapTree<int, uint64_t, page_size>::packLevel);
  using HeapTree_ = HeapTree<int, uint64_t, page_size>;
  for (size_t size = 2; size < 100000; size = 1.05 * size + 1) {
    for (int cacheLevel = 1; cacheLevel < GetLogBaseTwo(size) + 2;
         ++cacheLevel) {
      int maxLevel = GetLogBaseTwo(size - 1) + 2;
      std::vector<size_t> pathIdxs1(maxLevel);
      std::vector<size_t> pathIdxs2(maxLevel);
      for (int r = 0; r < 10000; ++r) {
        size_t path1 = UniformRandom(size - 1);
        size_t path2 = UniformRandom(size - 1);
        if (path1 == path2) {
          continue;
        }
        int commonSuffixLen = std::countr_zero(path1 ^ path2);
        int level1 = HeapTree_::GetPathIdx(pathIdxs1.begin(), pathIdxs1.end(),
                                           path1, size, cacheLevel);
        for (int i = 0; i < level1; ++i) {
          ASSERT_EQ(pathIdxs1[i],
                    HeapTree_::GetIdx(path1, size, i, maxLevel, cacheLevel));
          ASSERT_LT(pathIdxs1[i], size * 2 - 1);
          if (i > 0) {
            ASSERT_LT(pathIdxs1[i - 1], pathIdxs1[i]);
          }
        }
        int level2 = HeapTree_::GetPathIdx(pathIdxs2.begin(), pathIdxs2.end(),
                                           path2, size, cacheLevel);
        for (int i = 0; i < level2; ++i) {
          ASSERT_LT(pathIdxs2[i], size * 2 - 1);
          if (i > 0) {
            ASSERT_LT(pathIdxs2[i - 1], pathIdxs2[i]);
          }
        }
        for (int i = 0; i <= commonSuffixLen; ++i) {
          ASSERT_EQ(pathIdxs1[i], pathIdxs2[i]);
        }
        for (int i = commonSuffixLen + 1; i < level1; ++i) {
          for (int j = commonSuffixLen + 1; j < level2; ++j) {
            ASSERT_NE(pathIdxs1[i], pathIdxs2[j]);
          }
        }
      }
    }
  }
}

TEST(Basic, TestHeapTreeCorrectness) {
  testHeapTreeCorrectness<3>();
  testHeapTreeCorrectness<4>();
  testHeapTreeCorrectness<6>();
  testHeapTreeCorrectness<10>();
  testHeapTreeCorrectness<16>();
  testHeapTreeCorrectness<23>();
  testHeapTreeCorrectness<32>();
  testHeapTreeCorrectness<51>();
  testHeapTreeCorrectness<72>();
  testHeapTreeCorrectness<100>();
  testHeapTreeCorrectness<500>();
  testHeapTreeCorrectness<4096>();
}

TEST(Basic, testSamplingExpectation) {
  size_t size = 100;
  std::vector<size_t> counts(size);

  SampleFromPoisson(1.0, counts.begin(), counts.end());
  // SampleFromBinomial(1e-5, 1e5, counts.begin(), counts.end());
  // sum counts
  size_t sum = 0;
  for (size_t i = 0; i < size; ++i) {
    sum += counts[i];
  }
  printf("sum: %lu\n", sum);

  // for (size_t i = 0; i < size; ++i) {
  //   printf("%lu ", counts[i]);
  // }
  // printf("\n");
}

TEST(Basic, testSamplingPerf) {
  size_t size = 1e8;
  size_t sum = 0;
  NoReplaceSampler sampler(size);
  auto start = std::chrono::system_clock::now();
  for (size_t i = 0; i < size; ++i) {
    sum += sampler.Sample();
  }
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> diff = end - start;
  ASSERT_EQ(sum, size);
  printf("Time: %f ms\n", diff.count() * 1e3);
}

TEST(Basic, testSampler) {
  size_t size = 166667;
  size_t expectedSum = 83334;
  NoReplaceSampler sampler(expectedSum, size);

  size_t sum = 0;
  for (size_t i = 0; i < size; ++i) {
    sum += sampler.Sample();
  }
  ASSERT_EQ(sum, expectedSum);
}

TEST(Basic, testBuildBottomUp) {
  for (size_t leafCount = 71429; leafCount <= 71429; ++leafCount) {
    int totalLevel = GetLogBaseTwo(leafCount - 1) + 2;
    for (int cacheLevel = 1; cacheLevel <= totalLevel; ++cacheLevel) {
      // printf("leafCount: %lu, cacheLevel: %d, packLevel: %d\n", leafCount,
      //        cacheLevel, packLevel);
      HeapTree<uint64_t, uint64_t, 72> heapTree(leafCount, cacheLevel);
      for (auto it = heapTree.beginInternal(); it != heapTree.endInternal();
           ++it) {
        *it = 0;
      }
      size_t size = leafCount * 2 - 1;
      std::function<uint64_t(uint64_t&, const uint64_t&, const uint64_t&)>
          reduceFunc =
              [](uint64_t& val, const uint64_t& i,
                 const uint64_t& j) -> uint64_t { return val = i + j; };
      std::function<uint64_t(uint64_t&, size_t)> leafFunc =
          [](uint64_t& i, size_t idx) { return i = idx + 1; };
      ASSERT_EQ(heapTree.BuildBottomUp<uint64_t>(reduceFunc, leafFunc),
                leafCount * (leafCount + 1) / 2);
      for (auto it = heapTree.beginInternal(); it != heapTree.endInternal();
           ++it) {
        ASSERT_NE(*it, 0);
      }
    }
  }
}

TEST(Basic, testBuildBottomUpPerf) {
  for (size_t leafCount = 10000; leafCount < 100000000;
       leafCount = leafCount * 5 / 4) {
    int totalLevel = GetLogBaseTwo(leafCount - 1) + 2;
    size_t size = leafCount * 2 - 1;
    HeapTree<uint64_t> heapTree(leafCount);
    std::function<uint64_t(uint64_t&, const uint64_t&, const uint64_t&)>
        reduceFunc = [](uint64_t& val, const uint64_t& i,
                        const uint64_t& j) -> uint64_t { return val = i + j; };
    std::function<uint64_t(uint64_t&, size_t)> leafFunc =
        [](uint64_t& i, size_t idx) { return i = 1UL; };
    ASSERT_EQ(heapTree.BuildBottomUp(reduceFunc, leafFunc), leafCount);
  }
}
