#include <gtest/gtest.h>

#include "common/encrypted.hpp"
#include "common/mov_intrinsics.hpp"
#include "odsl/heap_tree.hpp"
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
  for (size_t i = 0; i < size; i++) {
    ASSERT_EQ(b.data[i], b2.data[i]);
  }
  size_t i = UniformRandom(size - 1);
  // Check changing the encrypted data gets different decrypted data
  e.data[i]++;
  TestBlock b3;
  e.Decrypt(b3);
  bool isDifferent = false;
  for (size_t i = 0; i < size; i++) {
    if (b.data[i] != b3.data[i]) {
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
  for (size_t i = 0; i < size; i++) {
    if (e.data[i] != e2.data[i]) {
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
  for (size_t i = 0; i < size; i++) {
    ASSERT_EQ(b.data[i], b4.data[i]);
  }
  FreshEncrypted<TestBlock> fe2;
  fe2.Encrypt(b, iv);
  for (size_t i = 0; i < size; i++) {
    ASSERT_EQ(fe.data[i], fe2.data[i]);
  }
  for (size_t i = 0; i < MAC_SIZE; i++) {
    ASSERT_EQ(fe.tag[i], fe2.tag[i]);
  }
  // Check changing the encrypted data can be detected
  i = UniformRandom(size - 1);
  fe.data[i]++;
  TestBlock b5;
  try {
    fe.Decrypt(b5, iv);
    ASSERT_TRUE(false);
  } catch (std::exception& e) {
    ASSERT_TRUE(true);
  }
  fe.data[i]--;
  // Check changing the tag can be detected
  i = UniformRandom(MAC_SIZE - 1);
  fe.tag[i]++;
  try {
    fe.Decrypt(b5, iv);
    ASSERT_TRUE(false);
  } catch (std::exception& e) {
    ASSERT_TRUE(true);
  }
  fe.tag[i]--;
  // Check changing the iv can be detected
  i = UniformRandom(IV_SIZE - 1);
  iv[i]++;
  try {
    fe.Decrypt(b5, iv);
    ASSERT_TRUE(false);
  } catch (std::exception& e) {
    ASSERT_TRUE(true);
  }
}

template <const uint64_t size>
void testEncryptPerf() {
  for (uint64_t round = 0; round < 1000000UL; ++round) {
    struct TestBlock {
      uint8_t data[size];
    };
    TestBlock b;
    for (uint64_t i = 0; i < size; i++) {
      b.data[i] = 7 * i + round;
    }
    uint8_t iv[IV_SIZE];
    FreshEncrypted<TestBlock> e2;
    e2.Encrypt(b, iv);
  }
}

TEST(Basic, Encrypted) {
  testEncrypted<1>();
  testEncrypted<2>();
  testEncrypted<3>();
  testEncrypted<7>();
  testEncrypted<8>();
  testEncrypted<9>();
  testEncrypted<12>();
  testEncrypted<15>();
  testEncrypted<16>();
  testEncrypted<17>();
  testEncrypted<30>();
  testEncrypted<31>();
  testEncrypted<32>();
  testEncrypted<33>();
  testEncrypted<63>();
  testEncrypted<64>();
  testEncrypted<65>();
  testEncrypted<127>();
  testEncrypted<128>();
  testEncrypted<129>();
  testEncrypted<200>();
  testEncrypted<256>();
  testEncrypted<257>();
  testEncrypted<300>();
  testEncrypted<512>();
  testEncrypted<513>();
  testEncrypted<600>();
  testEncrypted<1023>();
  testEncrypted<1024>();
  testEncrypted<1025>();
  testEncrypted<2000>();
  testEncrypted<4000>();
  testEncrypted<4095>();
  testEncrypted<4096>();
  testEncrypted<4097>();
  testEncrypted<7000>();
  testEncrypted<12345>();
  testEncrypted<16384>();
  testEncrypted<16385>();
  testEncrypted<20000>();
}

TEST(Basic, EncryptPerf) { testEncryptPerf<4096>(); }

template <const size_t page_size>
void testHeapTreeCorrectness() {
  // printf(
  //     "Testing heap tree correctness for page size: %lu, packing %d
  //     levels\n", page_size, HeapTree<int, uint64_t, page_size>::packLevel);
  using HeapTree_ = HeapTree<int, uint64_t, page_size>;
  for (size_t size = 2; size < 100000; size = 1.05 * size + 1) {
    int totalLevel = GetLogBaseTwo(size - 1) + 2;
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
        int level1 = HeapTree_::GetNodeIdxArr(pathIdxs1.begin(), path1, size,
                                              totalLevel, cacheLevel);
        for (int i = 0; i < level1; ++i) {
          ASSERT_LT(pathIdxs1[i], size * 2 - 1);
          if (i > 0) {
            ASSERT_LT(pathIdxs1[i - 1], pathIdxs1[i]);
          }
        }
        int level2 = HeapTree_::GetNodeIdxArr(pathIdxs2.begin(), path2, size,
                                              totalLevel, cacheLevel);
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