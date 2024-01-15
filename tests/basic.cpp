#include <gtest/gtest.h>

#include "common/mov_intrinsics.hpp"
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
