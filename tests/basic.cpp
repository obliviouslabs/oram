#include <gtest/gtest.h>

#include "common/mov_intrinsics.hpp"
#include "common/probability.hpp"
#include "oram/tree.hpp"
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
      size_t idx = HeapTree<int>::GetCAIdx(i, size, level, maxLevel);
      printf("level: %d, i: %d, idx: %lu\n", level, i, idx);
    }
  }
}

TEST(Basic, TestHeapTree2) {
  size_t size = 10;
  int maxLevel = GetLogBaseTwo(size - 1) + 2;
  std::vector<size_t> pathIdxs(maxLevel);
  for (size_t i = 0; i < size; ++i) {
    int level =
        HeapTree<int>::GetCAPathIdx(pathIdxs.begin(), pathIdxs.end(), i, size);
    printf("Path %lu:\n", i);
    for (int j = 0; j < level; ++j) {
      printf("%lu ", pathIdxs[j]);
    }
    printf("\n");
  }
}

TEST(Basic, TestHeapTree3) {
  for (size_t size = 2; size <= 1024; size++) {
    int maxLevel = GetLogBaseTwo(size - 1) + 2;
    for (int cacheLevel = 1; cacheLevel < maxLevel; ++cacheLevel) {
      for (int packLevel = 1; packLevel <= maxLevel; ++packLevel) {
        // printf("size: %lu, cacheLevel: %d, packLevel: %d\n", size,
        // cacheLevel,
        //        packLevel);
        std::vector<size_t> pathIdxs(maxLevel);
        std::vector<size_t> res;
        for (size_t i = 0; i < 1UL << (maxLevel - 1); ++i) {
          size_t ri = reverseBits(i, maxLevel - 1);
          if (ri >= size) {
            continue;
          }
          int level =
              HeapTree<int>::GetCAPathIdx(pathIdxs.begin(), pathIdxs.end(), ri,
                                          size, packLevel, cacheLevel);
          res.push_back(pathIdxs[level - 1]);
        }
        HeapTree<int> t(size, cacheLevel, packLevel);
        HeapTree<int>::ReverseLexLeafIndexer iter(t);
        for (size_t i = 0; i < size; ++i) {
          // printf("i: %lu\n", i);
          ASSERT_EQ(iter.getIndex(), res[i]);
          ++iter;
        }
      }
    }
  }
}

TEST(Basic, TestHeapTreeIndexerPerf) {
  size_t dummy = 0;
  // measure the time of the following code
  auto start = std::chrono::system_clock::now();
  for (size_t size = 2048; size <= 32768; size = size * 5 / 4) {
    int maxLevel = GetLogBaseTwo(size - 1) + 2;
    for (int cacheLevel = 1; cacheLevel <= maxLevel; ++cacheLevel) {
      for (int packLevel = 1; packLevel <= maxLevel; ++packLevel) {
        int totalLevel = GetLogBaseTwo(size - 1) + 2;
        for (size_t ri = 0; ri < size; ++ri) {
          size_t idx = HeapTree<int>::GetCAIdx(ri, totalLevel - 1, totalLevel,
                                               packLevel, cacheLevel);
          dummy += idx;
        }
      }
    }
  }
  auto end1 = std::chrono::system_clock::now();
  std::chrono::duration<double> diff1 = end1 - start;

  for (size_t size = 2048; size <= 8192; size = size * 5 / 4) {
    int maxLevel = GetLogBaseTwo(size - 1) + 2;
    for (int cacheLevel = 1; cacheLevel <= maxLevel; ++cacheLevel) {
      for (int packLevel = 1; packLevel <= maxLevel; ++packLevel) {
        HeapTree<int>::ReverseLexLeafIndexer iter(maxLevel, cacheLevel,
                                                  packLevel, size);
        for (size_t i = 0; i < size; ++i) {
          ++iter;
        }
      }
    }
  }
  auto end2 = std::chrono::system_clock::now();
  printf("Time naive: %f ms\n", diff1.count() * 1e3);
  std::chrono::duration<double> diff2 = end2 - end1;
  printf("Time indexer: %f ms\n", diff2.count() * 1e3);
  ASSERT_NE(dummy, 0);
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

TEST(Basic, testBuildBottomUp) {
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

TEST(Basic, testBuildBottomUpPerf) {
  for (size_t leafCount = 10000; leafCount < 100000000;
       leafCount = leafCount * 5 / 4) {
    int totalLevel = GetLogBaseTwo(leafCount - 1) + 2;
    size_t size = leafCount * 2 - 1;
    std::vector<uint64_t> heap(size);
    std::function<uint64_t(uint64_t&, const uint64_t&, const uint64_t&)>
        reduceFunc = [](uint64_t& val, const uint64_t& i,
                        const uint64_t& j) -> uint64_t { return val = i + j; };
    std::function<uint64_t(uint64_t&, size_t)> leafFunc =
        [](uint64_t& i, size_t idx) { return i = 1UL; };
    ASSERT_EQ(HeapTree<uint64_t>::BuildBottomUp(heap.begin(), heap.end(),
                                                reduceFunc, leafFunc),
              leafCount);
  }
}
