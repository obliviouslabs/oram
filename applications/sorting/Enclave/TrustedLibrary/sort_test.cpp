#include "../Enclave.h"
#include "Enclave_t.h"
#include "external_memory/algorithm/ca_bucket_sort.hpp"
#include "external_memory/algorithm/kway_butterfly_sort.hpp"
#include "external_memory/algorithm/kway_distri_sort.hpp"
#ifndef ALGO
#define ALGO KWAYBUTTERFLYOSORT
#endif
#ifndef MIN_SIZE
#define MIN_SIZE (1UL << 19)
#endif
#ifndef MAX_SIZE
#define MAX_SIZE (800 * (1UL << 19))
#endif
#ifndef STEP_RATIO
#define STEP_RATIO 1.2
#endif

using namespace EM::Algorithm;
using namespace EM::NonCachedVector;
EM::Backend::MemServerBackend* EM::Backend::g_DefaultBackend = nullptr;

void ecall_sort_perf() {
  // printf("ecall_sort_perf called\n");
  dbg_printf("run in debug mode\n");
  // using TestVector = EM::ExtVector::Vector<SortElement, 4032, true, 28000>;
  constexpr bool NeedCache = ALGO == ORSHUFFLE || ALGO == BITONICSORT ||
                             ALGO == UNOPTBITONICSORT || ALGO == BITONICSHUFFLE;
  constexpr size_t IdealPageSize = 16384;
#ifdef DISK_IO
  constexpr size_t MinRandomAccessPageSize = 8192;
#else
  constexpr size_t MinRandomAccessPageSize = 2048;
#endif
  constexpr size_t CachePageSize =
      ALGO != UNOPTBITONICSORT
          ? (ALGO != KWAYDISTRIBUTIONOSORT ? IdealPageSize
                                           : MinRandomAccessPageSize)
          : ELEMENT_SIZE + 32;
  constexpr size_t CacheSize =
      ALGO != BITONICSHUFFLE && ALGO != UNOPTBITONICSORT
          ? (uint64_t)(ENCLAVE_SIZE << 19) / CachePageSize
          : 1;
  using CachedVector = EM::ExtVector::Vector<SortElement, CachePageSize - 32,
                                             true, true, CacheSize>;
  using TestVector = std::conditional<NeedCache, CachedVector,
                                      Vector<SortElement, CachePageSize>>::type;
  for (double factor = 1; factor <= MAX_SIZE / MIN_SIZE; factor *= STEP_RATIO) {
    uint64_t size = (uint64_t)(factor * MIN_SIZE);
    if (EM::Backend::g_DefaultBackend) {
      delete EM::Backend::g_DefaultBackend;
    }
    size_t BackendSize =
        ((ALGO == CABUCKETSORT || ALGO == CABUCKETSHUFFLE) ? 16 : 4) *
        sizeof(TaggedT<SortElement>) * size;
    EM::Backend::g_DefaultBackend =
        new EM::Backend::MemServerBackend(BackendSize);
    TestVector vExt(size);
    TestVector::Writer writer(vExt.begin(), vExt.end());
    for (uint64_t i = 0; i < size; ++i) {
      SortElement element = SortElement();
      element.key = UniformRandom();
      writer.write(element);
    }
    writer.flush();
    uint64_t currTime;
    ocall_measure_time(&currTime);
    if constexpr (ALGO == KWAYBUTTERFLYOSORT) {
      KWayButterflySort(vExt);
    } else if constexpr (ALGO == KWAYBUTTERFLYOSHUFFLE) {
      KWayButterflyOShuffle(vExt);
    } else if constexpr (ALGO == CABUCKETSORT) {
      CABucketSort(vExt);
    } else if constexpr (ALGO == BITONICSORT || ALGO == UNOPTBITONICSORT) {
      BitonicSort(vExt);
    } else if constexpr (ALGO == ORSHUFFLE) {
      OrShuffle(vExt);
    } else if constexpr (ALGO == CABUCKETSHUFFLE) {
      CABucketShuffle(vExt);
    } else if constexpr (ALGO == BITONICSHUFFLE) {
      BitonicShuffle(vExt);
    } else if constexpr (ALGO == KWAYDISTRIBUTIONOSORT) {
      KWayDistriSort(vExt);
    } else if constexpr (ALGO == KWAYDISTRIBUTIONOSORTSHUFFLED) {
      KWayDistriSortShuffled(vExt);
    } else if constexpr (ALGO == EXTMERGESORT) {
      ExtMergeSort(vExt);
    }

    uint64_t currTime2;
    ocall_measure_time(&currTime2);
    uint64_t timediff = currTime2 - currTime;
    printf("%ld\t%ld\t%d.%d\n", size, ELEMENT_SIZE, timediff / 1'000'000'000,
           timediff % 1'000'000'000);
    if (timediff / 1'000'000'000 > 1800) {
      break;
    }
  }
}

void ecall_pageswap_with_crypt_perf() {
  printf("ecall_pageswap_with_crypt_perf called\n");
  uint64_t size = (1UL << 24);
  uint64_t currTime, currTime2, currTime3;

  uint64_t r = 0;
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
  EM::Backend::g_DefaultBackend = new EM::Backend::MemServerBackend(1UL << 34);
  printf("prepare to create vector\n");
  Vector<SortElement> vExt(size);
  ocall_measure_time(&currTime);
  for (int round = 0; round < 8; ++round) {
    Vector<SortElement>::Writer writer(vExt.begin(), vExt.end());
    for (uint64_t i = 0; i < size; ++i) {
      SortElement element = SortElement();
      element.key = UniformRandom();
      printf("%ld\n", element.key);
      writer.write(element);
    }
  }
  ocall_measure_time(&currTime2);
  for (int round = 0; round < 8; ++round) {
    Vector<SortElement>::PrefetchReader reader(vExt.begin(), vExt.end());
    for (uint64_t i = 0; i < size; ++i) {
      r ^= reader.read().key;
    }
  }
  ocall_measure_time(&currTime3);

  printf("r = %ld\n", r);
  uint64_t timediff1 = currTime2 - currTime;
  uint64_t timediff2 = currTime3 - currTime2;
  printf("Elapsed: %d.%d s for write pages\n", timediff1 / 1'000'000'000,
         timediff1 % 1'000'000'000);
  printf("Elapsed: %d.%d s for read pages\n", timediff2 / 1'000'000'000,
         timediff2 % 1'000'000'000);
}

void ecall_pageswap_perf() {
  printf("ecall_pageswap_perf\n");
  uint64_t currTime, currTime2, currTime3, currTime4, currTime5, currTime6;

  uint64_t r = 0;
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
  EM::Backend::g_DefaultBackend = new EM::Backend::MemServerBackend(1UL << 35);
  uint8_t buf[4096] = {0};
  ocall_measure_time(&currTime);
  for (uint64_t round = 0; round < (1UL << 22); ++round) {
    ocall_Write(round * 4096UL, 4096, buf);
  }
  ocall_measure_time(&currTime2);
  for (uint64_t round = 0; round < (1UL << 22); ++round) {
    ocall_Read(round * 4096UL, 4096, buf);
    r ^= buf[0];
  }
  ocall_measure_time(&currTime3);

  printf("ocall per page r = %ld\n", r);
  uint64_t timediff1 = currTime2 - currTime;
  uint64_t timediff2 = currTime3 - currTime2;
  printf("Elapsed: %d.%d s for write pages\n", timediff1 / 1'000'000'000,
         timediff1 % 1'000'000'000);
  printf("Elapsed: %d.%d s for read pages\n", timediff2 / 1'000'000'000,
         timediff2 % 1'000'000'000);
#ifndef DISK_IO
  ocall_measure_time(&currTime);
  for (uint64_t round = 0; round < (1UL << 22); ++round) {
    EM::Backend::g_DefaultBackend->ocall_Write(round * 4096UL, 4096, buf);
  }
  ocall_measure_time(&currTime2);
  for (uint64_t round = 0; round < (1UL << 22); ++round) {
    EM::Backend::g_DefaultBackend->ocall_Read(round * 4096UL, 4096, buf);
    r ^= buf[0];
  }
  ocall_measure_time(&currTime3);

  printf("direct mem access r = %ld\n", r);
  timediff1 = currTime2 - currTime;
  timediff2 = currTime3 - currTime2;
  printf("Elapsed: %d.%d s for write pages\n", timediff1 / 1'000'000'000,
         timediff1 % 1'000'000'000);
  printf("Elapsed: %d.%d s for read pages\n", timediff2 / 1'000'000'000,
         timediff2 % 1'000'000'000);
#endif

  static const uint64_t batchSize = 16;
  uint8_t batchBuf[4096 * batchSize] = {0};
  uint64_t offsets[batchSize];
  uint64_t sizes[batchSize];
  for (size_t i = 0; i < batchSize; ++i) {
    sizes[i] = 4096;
  }
  ocall_measure_time(&currTime4);
  for (uint64_t round = 0; round < (1UL << 22); round += batchSize) {
    for (size_t i = 0; i < batchSize; ++i) {
      offsets[i] = (round + i) * 4096UL;
    }
    ocall_Write_Batch(offsets, sizes, batchBuf, batchSize, 4096 * batchSize);
  }
  ocall_measure_time(&currTime5);
  for (uint64_t round = 0; round < (1UL << 22); round += batchSize) {
    for (size_t i = 0; i < batchSize; ++i) {
      offsets[i] = (round + i) * 4096UL;
    }
    ocall_Read_Batch(offsets, sizes, batchBuf, batchSize, 4096 * batchSize);
    r ^= batchBuf[0];
  }
  ocall_measure_time(&currTime6);
  printf("pack ocall r = %ld\n", r);
  uint64_t timediff3 = currTime5 - currTime4;
  uint64_t timediff4 = currTime6 - currTime5;
  printf("Elapsed: %d.%d s for write pages in batch\n",
         timediff3 / 1'000'000'000, timediff3 % 1'000'000'000);
  printf("Elapsed: %d.%d s for read pages in batch\n",
         timediff4 / 1'000'000'000, timediff4 % 1'000'000'000);
}

void ecall_test_sanity() {
  uint64_t currTime, currTime2;
  SortElement* pa = new SortElement();
  SortElement* pb = new SortElement();
  SortElement& a = *pa;
  SortElement& b = *pb;
  a.key = 1;
  b.key = 2;
  for (size_t i = 0; i < sizeof(a.payload); ++i) {
    a.payload[i] = i % 128;
    b.payload[i] = (17 + i * i) % 128;
  }
  condSwap(false, a, b);
  if (!(a.key == 1 && b.key == 2)) {
    printf("wrong key for not swap\n");
    abort();
  }
  for (size_t i = 0; i < sizeof(a.payload); ++i) {
    if (a.payload[i] != i % 128) {
      printf("wrong a.payload for not swap\n");
      abort();
    }
    if (b.payload[i] != (17 + i * i) % 128) {
      printf("wrong b.payload for not swap\n");
      abort();
    }
  }
  uint64_t round = 1000000000;
  ocall_measure_time(&currTime);
  for (size_t i = 0; i < round; ++i) {
    condSwap(i, a, b);
  }
  ocall_measure_time(&currTime2);
  if (!(b.key == 1 && a.key == 2)) {
    printf("wrong key for swap\n");
    abort();
  }
  for (size_t i = 0; i < sizeof(a.payload); ++i) {
    if (b.payload[i] != i % 128) {
      printf("wrong b.payload for swap\n");
      abort();
    }
    if (a.payload[i] != (17 + i * i) % 128) {
      printf("wrong a.payload for swap\n");
      abort();
    }
  }
  printf("sanity test pass\n");
  uint64_t timediff = currTime2 - currTime;
  printf("%ld\t%d.%d\n", round, timediff / 1'000'000'000,
         timediff % 1'000'000'000);
}

void mergesplit_tag_perf() {
  printf("{\n{");
  for (uint64_t Z = 256; Z <= 16384; Z *= 2) {
    for (size_t way = 2; way <= 8; ++way) {
      uint64_t start, end;
      TaggedT<SortElement> defaultVal;
      defaultVal.setDummy();
      size_t bucketCount =
          (ENCLAVE_SIZE << 19) / sizeof(TaggedT<SortElement>) / Z;
      const uint64_t N = way * Z * (bucketCount / way);
      std::vector<TaggedT<SortElement>> v(N, defaultVal);
      ocall_measure_time(&start);
      for (size_t round = 0; round < 5000 / way; ++round) {
        size_t offset = round % (bucketCount / way);
        std::vector<TaggedT<SortElement>>::iterator begins[8];
        for (size_t i = 0; i < way; ++i) {
          begins[i] = v.begin() + i * Z + Z * offset * way;
        }
        MergeSplitKWay(begins, way, Z);
      }
      ocall_measure_time(&end);
      printf("%f", 1e-9 * (end - start));
      if (way != 8) {
        printf(", ");
      }
    }
    printf("}");
    if (Z != 16384) {
      printf(", \n{");
    } else {
      printf("}");
    }
  }
}

void mergesplit_pivot_perf() {
  printf("{\n{");
  for (uint64_t Z = 32768; Z <= 65536; Z *= 2) {
    for (size_t way = 2; way <= 8; ++way) {
      std::vector<Block<SortElement>> pivots(way - 1);
      uint64_t start, end;
      Block<SortElement> defaultVal;
      defaultVal.setDummy();
      size_t bucketCount =
          (ENCLAVE_SIZE << 19) / sizeof(Block<SortElement>) / Z;
      const uint64_t N = way * Z * (bucketCount / way);
      std::vector<Block<SortElement>> v(N, defaultVal);
      ocall_measure_time(&start);
      for (size_t round = 0; round < 5000 / way; ++round) {
        size_t offset = round % (bucketCount / way);
        std::vector<Block<SortElement>>::iterator begins[8];
        for (size_t i = 0; i < way; ++i) {
          begins[i] = v.begin() + i * Z + Z * offset * way;
        }

        MergeSplitKWay(begins, way, Z, pivots.begin());
      }
      ocall_measure_time(&end);
      printf("%f", 1e-9 * (end - start));
      if (way != 8) {
        printf(", ");
      }
    }
    printf("}");
    if (Z != 65536) {
      printf(", \n{");
    } else {
      printf("}");
    }
  }
}

void ecall_mergesplit_perf() {
  // mergesplit_tag_perf();
  mergesplit_pivot_perf();
}

template <class Iterator>
void multiWayUsingTwoWay(const std::vector<Iterator>& begins, size_t Z,
                         PartitionMethod method) {
  size_t k = begins.size();

  if (k == 1) {
    return;
  }
  if (k == GetNextPowerOfTwo(k)) {
    std::vector<Iterator> left(begins.begin(), begins.begin() + k / 2);
    std::vector<Iterator> right(begins.begin() + k / 2, begins.end());
    multiWayUsingTwoWay(left, Z, method);
    multiWayUsingTwoWay(right, Z, method);
    for (size_t i = 0; i < k / 2; ++i) {
      MergeSplitTwoWay(left[i], right[i], Z, k / 2, method);
    }
  } else {  // this implementation is incorrect, but the runtime should be the
            // same as a correct implementation
    MergeSplitTwoWay(begins[0], begins[0] + Z * k / 2, Z * k / 2,
                     GetNextPowerOfTwo(k) / 2, method);
    std::vector<Iterator> left(begins.begin(), begins.begin() + k / 2);
    std::vector<Iterator> right(begins.begin() + k / 2, begins.end());
    multiWayUsingTwoWay(left, Z, method);
    multiWayUsingTwoWay(right, Z, method);
  }
}

void ecall_mergesplit_compare() {
  uint64_t Z = 4096;
  for (size_t way = 2; way <= 8; ++way) {
    uint64_t start, end;
    TaggedT<SortElement> defaultVal;
    defaultVal.setDummy();
    const uint64_t N = way * Z;
    std::vector<TaggedT<SortElement>> v(N, defaultVal);
    std::vector<TaggedT<SortElement>>::iterator begins[8];
    for (size_t i = 0; i < way; ++i) {
      begins[i] = v.begin() + i * Z;
    }
    ocall_measure_time(&start);

    for (size_t round = 0; round < 10000; ++round) {
      MergeSplitKWay(begins, way, Z);
    }
    ocall_measure_time(&end);
    printf("%f", 1e-9 * (end - start));
    if (way != 8) {
      printf(", ");
    } else {
      printf("\n");
    }
  }

  for (size_t way = 2; way <= 8; ++way) {
    uint64_t start, end;
    TaggedT<SortElement> defaultVal;
    defaultVal.setDummy();
    const uint64_t N = way * Z;
    std::vector<TaggedT<SortElement>> v(N, defaultVal);
    std::vector<std::vector<TaggedT<SortElement>>::iterator> begins;
    for (size_t i = 0; i < way; ++i) {
      begins.push_back(v.begin() + i * Z);
    }
    ocall_measure_time(&start);

    for (size_t round = 0; round < 10000; ++round) {
      multiWayUsingTwoWay(begins, Z, OR_COMPACT);
    }
    ocall_measure_time(&end);
    printf("%f", 1e-9 * (end - start));
    if (way != 8) {
      printf(", ");
    } else {
      printf("\n");
    }
  }

  // bitonic
  for (size_t way = 2; way <= 8; ++way) {
    uint64_t start, end;
    TaggedT<SortElement> defaultVal;
    defaultVal.setDummy();
    const uint64_t N = way * Z;
    std::vector<TaggedT<SortElement>> v(N, defaultVal);
    std::vector<TaggedT<SortElement>> temp(N, defaultVal);

    ocall_measure_time(&start);

    for (size_t round = 0; round < 1000;
         ++round) {  // 1000 times rather than 10000 times for speed
      std::memcpy(&temp[0], &v[0], N * sizeof(TaggedT<SortElement>));
      binPacking<SortElement>(temp.begin(), temp.end(), temp.begin(), Z,
                              1UL << 63);
      std::memcpy(&v[0], &temp[0], N * sizeof(TaggedT<SortElement>));
    }

    ocall_measure_time(&end);
    printf("%f", 1e-9 * (end - start));
    if (way != 8) {
      printf(", ");
    } else {
      printf("\n");
    }
  }
}

void ecall_bitonic_perf() {
  for (uint64_t Z = 256; Z <= 65536; Z *= 2) {
    uint64_t start, end;
    TaggedT<SortElement> defaultVal;
    defaultVal.setDummy();
    size_t bucketCount =
        (ENCLAVE_SIZE << 19) / sizeof(TaggedT<SortElement>) / Z;
    const uint64_t N = Z * bucketCount;
    std::vector<TaggedT<SortElement>> v(N, defaultVal);
    ocall_measure_time(&start);
    static const auto cmp = [](const auto& t1, const auto& t2) {
      return t1.tag < t2.tag;
    };
    for (size_t round = 0; round < 5000; ++round) {
      size_t offset = round % bucketCount;
      BitonicSort(v.begin() + offset * Z, v.begin() + (offset + 1) * Z, cmp);
    }
    ocall_measure_time(&end);
    printf("%f, ", 1e-9 * (end - start));
  }
}