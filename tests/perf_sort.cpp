#include <gtest/gtest.h>

#include <functional>
#include <unordered_map>

#include "external_memory/algorithm/ca_bucket_sort.hpp"
#include "external_memory/algorithm/kway_butterfly_sort.hpp"
#include "external_memory/algorithm/kway_distri_sort.hpp"
#include "external_memory/algorithm/randomized_shell_sort.hpp"
#include "testutils.hpp"

#define PAGE_SIZE 4096

#define test_sort(n, f, ...)                         \
  _test_sort<Vector<SortElement, PAGE_SIZE>, false>( \
      n, #f, [](auto&& arr) { f(arr); } __VA_OPT__(, ) __VA_ARGS__)
#define test_sort_cache(n, f, Vec, ...)                         \
  _test_sort<Vec<SortElement, 16320, true, true, 1024>, false>( \
      n, #f, [](auto&& arr) { f(arr); } __VA_OPT__(, ) __VA_ARGS__)
#define test_sort_internal(n, f, ...)       \
  _test_sort<StdVector<SortElement>, true>( \
      n, #f, [](auto&& arr) { f(arr); } __VA_OPT__(, ) __VA_ARGS__)
#define ENABLE_PROFILING

using namespace EM::Algorithm;
using namespace EM::NonCachedVector;
using namespace std;

void printProfile(uint64_t N, ostream& ofs, auto& diff) {
  PERFCTR_LOG_TO(ofs);
  ofs << "N: " << N << std::endl;
  ofs << "Element size (bytes): " << sizeof(SortElement) << std::endl;
  ofs << "time (s): " << setw(9) << diff.count() << std::endl;
}

template <typename Vec, const bool isInternal, typename F>
void _test_sort(size_t size, string funcname, F&& sortFunc,
                bool isPermutation = false) {
  delete EM::Backend::g_DefaultBackend;
  EM::Backend::g_DefaultBackend =
      new EM::Backend::MemServerBackend((1ULL << 10) * (size + 1024));
  srand(time(NULL));
  const uint64_t N = size;
  cout << "test " << funcname << " perf on input size " << N << endl;
  StdVector<SortElement> v(N);
  unordered_map<uint64_t, int> value_count;
  for (uint64_t i = 0; i < N; i++) {
    v[i].key = UniformRandom();
    memset(v[i].payload, (char)v[i].key, sizeof(SortElement::payload));
    ++value_count[v[i].key];
  }
  auto cmp = [](const SortElement& ele1, const SortElement& ele2) {
    return ele1.key < ele2.key;
  };
  if constexpr (!isInternal) {
    Vec vExt(N);
    if constexpr (!Vec::useStdCopy) {
      CopyOut(v.begin(), v.end(), vExt.begin());
    } else {
      std::copy(v.begin(), v.end(), vExt.begin());
    }

    PERFCTR_RESET();
    auto start = std::chrono::system_clock::now();
    sortFunc(vExt);

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start;
    printProfile(N, cout, diff);
    if constexpr (!Vec::useStdCopy) {
      CopyIn(vExt.begin(), vExt.end(), v.begin(), 1);
    } else {
      std::copy(vExt.begin(), vExt.end(), v.begin());
    }
  } else {
    PERFCTR_RESET();
    auto start = std::chrono::system_clock::now();
    sortFunc(v);

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end - start;
    printProfile(N, cout, diff);
  }

  // check it's a permutation
  for (uint64_t i = 0; i < N; ++i) {
    auto key = v[i].key;
    if (value_count[key] == 0) {
      printf("index %ld contains duplicate key %ld\n", i, key);
      Assert(false);
    }
    char keyLastByte = (char)key;
    for (size_t j = 0; j < sizeof(SortElement::payload); ++j) {
      ASSERT_EQ((char)v[i].payload[j], keyLastByte);
    }
    ASSERT_GE(--value_count[key], 0);
  }
  if (!isPermutation) {
    // check increasing order
    for (uint64_t i = 0; i < N - 1; ++i) {
      ASSERT_TRUE(v[i].key <= v[i + 1].key);
    }
  }
}

TEST(TestSort, TestKWayButterflySortPerf) {
  // RELEASE_ONLY_TEST();
  for (double N = 1; N < 10000000; N *= 5) {
    test_sort((size_t)N, KWayButterflySort);
  }
}

TEST(TestSortInternal, TestKWayButterflySortPerf) {
  // RELEASE_ONLY_TEST();
  for (double N = 1; N < 100000000; N *= 5) {
    test_sort_internal((size_t)N, KWayButterflySortInternal, false);
  }
}

TEST(TestSort, TestKWayButterflyOShufflePerf) {
  // RELEASE_ONLY_TEST();
  for (double N = 1; N < 10000000; N *= 5) {
    test_sort((size_t)N, KWayButterflyOShuffle, true);
  }
}

TEST(TestSortInternal, TestKWayButterflyOShufflePerf) {
  // RELEASE_ONLY_TEST();
  for (double N = 1; N < 100000000; N *= 5) {
    test_sort_internal((size_t)N, KWayButterflyOShuffleInternal, true);
  }
}

TEST(TestSort, TestKWayDistriSortPerf) {
  // RELEASE_ONLY_TEST();
  for (double N = 1000000; N < 10000000; N *= 5) {
    test_sort((size_t)N, KWayDistriSort);
  }
}

TEST(TestSort, TestKWayDistriSortShuffledPerf) {
  // RELEASE_ONLY_TEST();
  for (double N = 1000000; N < 10000000; N *= 5) {
    test_sort((size_t)N, KWayDistriSortShuffled);
  }
}

TEST(TestSort, TestCacheObliviouBucketSortPerf) {
  // RELEASE_ONLY_TEST();
  test_sort(200000, CABucketSort);
}

TEST(TestSort, TestExtMergeSortPerf) {
  // RELEASE_ONLY_TEST();
  test_sort(10000000, ExtMergeSort);
}

TEST(TestSort, TestCacheObliviouBucketPermutationPerf) {
  // RELEASE_ONLY_TEST();
  test_sort(200000, CABucketShuffle, true);
}

TEST(TestSort, TestBitonicObliviousSortPerf) {
  // RELEASE_ONLY_TEST();
  for (double N = 1; N < 10000000; N *= 5) {
    test_sort_cache(N, BitonicSort, EM::ExtVector::Vector);
  }
}

TEST(TestSortInternal, TestRandomizedShellSort) {
  // RELEASE_ONLY_TEST();
  for (double N = 1; N < 100000000; N *= 2) {
    test_sort_internal((size_t)N, RandomizedShellSort);
  }
}

TEST(TestSortInternal, TestBitonicObliviousSortPerf) {
  // RELEASE_ONLY_TEST();
  for (double N = 1; N < 100000000; N *= 5) {
    test_sort_internal(N, BitonicSort, false);
  }
}

TEST(TestSort, TestOrShufflePerf) {
  // RELEASE_ONLY_TEST();
  test_sort_cache(200000, OrShuffle, EM::ExtVector::Vector, true);
}

TEST(TestSort, TestBitonicShufflePerf) {
  // RELEASE_ONLY_TEST();
  test_sort_cache(20000, BitonicShuffle, EM::ExtVector::Vector, true);
}

TEST(TestSort, TestSequentialPassPerf) {
  RELEASE_ONLY_TEST();
  srand(time(NULL));
  delete EM::Backend::g_DefaultBackend;
  EM::Backend::g_DefaultBackend =
      new EM::Backend::MemServerBackend((1ULL << 34));

  uint64_t logN = 24;
  const uint64_t N = 1UL << logN;
  cout << "test bucket oblivious permutation perf " << N << endl;
  cout << logN << '\t' << N << '\t';
  SortElement defaultval;
  vector<SortElement> v(N, defaultval);
  Vector<SortElement> vExt(N);
  PERFCTR_RESET();
  auto start = std::chrono::system_clock::now();

  CopyOut(v.begin(), v.end(), vExt.begin());
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> diff = end - start;
  printProfile(N, cout, diff);
  PERFCTR_RESET();
  start = std::chrono::system_clock::now();
  uint64_t xorsum = 0;
  CopyIn(vExt.begin(), vExt.end(), v.begin());
  end = std::chrono::system_clock::now();
  diff = end - start;
  std::cout << xorsum << std::endl;
  printProfile(N, cout, diff);
}

TEST(TestSort, TestMergeSplitPerf) {
  RELEASE_ONLY_TEST();
  PERFCTR_RESET();
  auto start = std::chrono::system_clock::now();
  uint64_t bitMask = 4UL;
  for (uint64_t Z = 512; Z <= 4096; Z *= 2) {
    TaggedT<SortElement> defaultVal;
    defaultVal.setDummy();
    size_t bucketCount = 2000;
    const uint64_t N = 2 * Z * (bucketCount / 2);
    vector<TaggedT<SortElement>> v(N, defaultVal);

    for (int round = 0; round < bucketCount / 2; ++round) {
      MergeSplitTwoWay(v.begin() + 2 * round, v.begin() + 2 * round + 1, Z,
                       bitMask);
    }
  }
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> diff = end - start;
  std::cout << "time (s): " << setw(9) << diff.count() << std::endl;

  PERFCTR_LOG();
}

TEST(TestSort, TestKWayMergeSplitPerf) {
  RELEASE_ONLY_TEST();
  PERFCTR_RESET();
  for (uint64_t Z = 256; Z <= 16384; Z *= 2) {
    for (size_t way = 2; way <= 8; ++way) {
      auto start = std::chrono::system_clock::now();
      TaggedT<SortElement> defaultVal;
      defaultVal.setDummy();
      size_t bucketCount = 2000;
      const uint64_t N = way * Z * (bucketCount / way);
      vector<TaggedT<SortElement>> v(N, defaultVal);
      for (int round = 0; round < bucketCount / way; ++round) {
        vector<TaggedT<SortElement>>::iterator begins[8];
        for (size_t i = 0; i < way; ++i) {
          begins[i] = v.begin() + i * Z + Z * round * way;
        }
        MergeSplitKWay(begins, way, Z);
      }

      auto end = std::chrono::system_clock::now();
      std::chrono::duration<double> diff = end - start;
      std::cout << "way=" << way << " Z=" << Z << " time (s): " << setw(9)
                << diff.count() << std::endl;
    }
  }
}

TEST(TestSort, TestKWayInterleaveSepMarksPerf) {
  RELEASE_ONLY_TEST();
  // PERFCTR_RESET();
  PROFILER_SET(false);
  for (uint64_t Z = 256; Z <= 16384; Z *= 2) {
    for (size_t k = 3; k <= 8; ++k) {
      auto start = std::chrono::system_clock::now();

      TaggedT<SortElement> defaultVal;
      defaultVal.setDummy();

      const uint64_t N = k * Z;
      vector<TaggedT<SortElement>> v(N, defaultVal);
      for (int round = 0; round < 2000 / k; ++round) {
        const auto getMark = [k = k](const auto& a) {
          return (uint16_t)(a.tag) % k;
        };
        vector<uint8_t> marks(N);
        for (size_t i = 0; i < N; ++i) {
          marks[i] = getMark(v[i]);
        }
        Interleave(v.begin(), v.end(), marks.begin(), marks.end(), k);
      }
      auto end = std::chrono::system_clock::now();
      std::chrono::duration<double> diff = end - start;
      std::cout << "k=" << k << " Z=" << Z << " time (s): " << setw(9)
                << diff.count() << std::endl;
      // PERFCTR_LOG();
    }
  }
}

TEST(TestSort, BNSort) {
  constexpr int i = 8;
  PERFCTR_RESET();
  vector<uint64_t> v(i);
  vector<uint8_t> marks(i);
  static constexpr StaticSort<i> staticSort;
  auto iter_pair = std::make_pair(v.begin(), marks.begin());
  staticSort(iter_pair);
  PERFCTR_LOG();
}

TEST(TestSort, TestMergeSplitPerfForDifferentBlockSizes) {
  RELEASE_ONLY_TEST();

  size_t Zbegin = 200, Zend = 2000;
  uint64_t runtime[4][Zend - Zbegin];
  uint64_t bitMask = 4UL;
  for (uint64_t Z = Zbegin; Z < Zend; Z = Z * 3 / 2) {
    TaggedT<uint64_t> defaultVal;
    defaultVal.setDummy();

    const uint64_t N = 2 * Z;
    vector<TaggedT<uint64_t>> v(N, defaultVal);
    for (int method = 0; method != 4; ++method) {
      PERFCTR_RESET();
      MergeSplitInPlace(v.begin(), v.end(), bitMask, (PartitionMethod)method);
      runtime[method][Z - Zbegin] = g_PerfCounters.swapCount;
    }
  }
  ofstream ofs;
  ofs.open("mergesplitperf.out");
  ofs << "\tEvenOdd\tOrCompact\tGoodrich\tBitonic\t" << std::endl;
  for (uint64_t Z = Zbegin; Z < Zend; ++Z) {
    ofs << Z << '\t';
    for (int method = 0; method != 4; ++method) {
      ofs << runtime[method][Z - Zbegin] << "\t";
    }
    ofs << std::endl;
  }
  ofs.close();
}