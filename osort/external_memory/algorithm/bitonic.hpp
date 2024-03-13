#pragma once
#include <vector>

#include "external_memory/extemvector.hpp"
#include "external_memory/stdvector.hpp"
#include "sort_def.hpp"
#include "static_sort.hpp"

namespace EM::Algorithm {

template <class Iterator, typename Compare>
void SlowBitonicSort(Iterator begin, Iterator end, Compare cmp) {
  uint64_t n = end - begin;
  Assert(GetNextPowerOfTwo(n) == n, n);
  for (uint64_t k = 2; k <= n; k *= 2) {
    for (uint64_t j = k / 2; j > 0; j /= 2) {
      for (uint64_t i = 0; i < n; ++i) {
        uint64_t l = i ^ j;
        if (l > i) {
          bool llti = cmp(*(begin + l), *(begin + i));
          bool swapFlag = (((i & k) == 0) * llti) + (((i & k) != 0) * (!llti));
          condSwap(swapFlag, *(begin + i), *(begin + l));
        }
      }
    }
  }
}

template <class Iterator, typename Compare>
void BitonicMergePow2(Iterator begin, Iterator end, Compare cmp, bool dire) {
  size_t size = end - begin;
  if (size > 1) {
    size_t halfSize = size / 2;
    Iterator leftIt = begin;
    Iterator rightIt = begin + halfSize;
    for (size_t i = 0; i < halfSize; ++i, ++leftIt, ++rightIt) {
      condSwap(dire != cmp(*leftIt, *rightIt), *leftIt, *rightIt);
    }
    BitonicMergePow2(begin, leftIt, cmp, dire);
    BitonicMergePow2(leftIt, end, cmp, dire);
  }
}

template <class Iterator, typename Compare>
void BitonicMerge(Iterator begin, Iterator end, Compare cmp, bool dire) {
  size_t size = end - begin;
  if (size > 1) {
    size_t halfSize = GetNextPowerOfTwo(size) / 2;
    Iterator leftIt = begin;
    Iterator rightIt = begin + halfSize;
    for (size_t i = 0; i < size - halfSize; ++i, ++leftIt, ++rightIt) {
      condSwap(dire != cmp(*leftIt, *rightIt), *leftIt, *rightIt);
    }
    Iterator midIt = begin + halfSize;
    BitonicMergePow2(begin, midIt, cmp, dire);
    BitonicMerge(midIt, end, cmp, dire);
  }
}

template <class Iterator, typename Compare>
void BitonicSort(Iterator begin, Iterator end, Compare cmp, bool dire) {
  size_t size = end - begin;
  if (size > 1) {
    size_t halfSize = size / 2;
    Iterator mid = begin + halfSize;
    BitonicSort(begin, mid, cmp, !dire);
    BitonicSort(mid, end, cmp, dire);
    BitonicMerge(begin, end, cmp, dire);
  }
}

template <class Iterator, typename Compare>
void BitonicSort(Iterator begin, Iterator end, Compare cmp) {
  BitonicSort(begin, end, cmp, true);
}

template <class Iterator>
void BitonicSort(Iterator begin, Iterator end) {
  auto cmp = [](const auto& a, const auto& b) { return a < b; };
  BitonicSort(begin, end, cmp, true);
}

template <typename Vec, typename Compare>
void BitonicSort(Vec& v, Compare cmp) {
  BitonicSort(v.begin(), v.end(), cmp);
}

template <typename Vec>
void BitonicSort(Vec& v) {
  auto cmp = [](const auto& a, const auto& b) { return a < b; };
  BitonicSort(v, cmp);
}

template <class Iterator>
void BitonicShuffle(Iterator begin, Iterator end) {
  using T = typename std::iterator_traits<Iterator>::value_type;
  auto cmpTag = [](const TaggedT<T>& a, const TaggedT<T>& b) {
    return a.tag < b.tag;
  };
  if constexpr (std::is_same<Iterator,
                             typename StdVector<T>::Iterator>::value) {
    std::vector<TaggedT<T>> taggedData(end - begin);
    size_t idx = 0;
    for (auto it = begin; it != end; ++it, ++idx) {
      taggedData[idx].v = *it;
      taggedData[idx].tag = UniformRandom();
    }
    BitonicSort(taggedData.begin(), taggedData.end(), cmpTag);
    auto taggedIt = taggedData.begin();
    for (auto it = begin; it != end; ++it, ++taggedIt) {
      *it = taggedIt->v;
    }
  } else {
    constexpr size_t CachePageSize =
        4096 / sizeof(TaggedT<T>) * sizeof(TaggedT<T>);
    constexpr size_t CacheSize = DEFAULT_HEAP_SIZE / CachePageSize - 1;
    EM::ExtVector::Vector<TaggedT<T>, CachePageSize, true, true, CacheSize>
        taggedData(end - begin);
    size_t idx = 0;
    for (auto it = begin; it != end; ++it, ++idx) {
      const auto& it_const = it;
      taggedData.AtForLateInit(idx) = {UniformRandom(), *it_const};
    }
    BitonicSort(taggedData.begin(), taggedData.end(), cmpTag);
    auto taggedIt = taggedData.begin();
    for (auto it = begin; it != end; ++it, ++taggedIt) {
      *it = taggedIt.derefNoWriteBack().v;
    }
  }
}

template <typename Vec>
void BitonicShuffle(Vec& v) {
  BitonicShuffle(v.begin(), v.end());
}

template <class KeyIterator, class PayloadIterator>
void BitonicMergePow2SepPayload(KeyIterator keyBegin, KeyIterator keyEnd,
                                PayloadIterator payloadBegin, bool dire) {
  size_t size = keyEnd - keyBegin;
  if (size > 1) {
    size_t halfSize = size / 2;
    KeyIterator leftKeyIt = keyBegin;
    KeyIterator rightKeyIt = keyBegin + halfSize;
    PayloadIterator leftPayloadIt = payloadBegin;
    PayloadIterator rightPayloadIt = payloadBegin + halfSize;
    for (size_t i = 0; i < halfSize;
         ++i, ++leftKeyIt, ++rightKeyIt, ++leftPayloadIt, ++rightPayloadIt) {
      bool swapFlag = dire != (*leftKeyIt < *rightKeyIt);
      condSwap(swapFlag, *leftKeyIt, *rightKeyIt);
      condSwap(swapFlag, *leftPayloadIt, *rightPayloadIt);
    }
    BitonicMergePow2SepPayload(keyBegin, leftKeyIt, payloadBegin, dire);
    BitonicMergePow2SepPayload(leftKeyIt, keyEnd, leftPayloadIt, dire);
  }
}

template <class KeyIterator, class PayloadIterator>
void BitonicMergeSepPayload(KeyIterator keyBegin, KeyIterator keyEnd,
                            PayloadIterator payloadBegin, bool dire) {
  size_t size = keyEnd - keyBegin;
  if (size > 1) {
    size_t halfSize = GetNextPowerOfTwo(size) / 2;
    KeyIterator leftKeyIt = keyBegin;
    KeyIterator rightKeyIt = keyBegin + halfSize;
    PayloadIterator leftPayloadIt = payloadBegin;
    PayloadIterator rightPayloadIt = payloadBegin + halfSize;
    for (size_t i = 0; i < size - halfSize;
         ++i, ++leftKeyIt, ++rightKeyIt, ++leftPayloadIt, ++rightPayloadIt) {
      bool swapFlag = dire != (*leftKeyIt < *rightKeyIt);
      condSwap(swapFlag, *leftKeyIt, *rightKeyIt);
      condSwap(swapFlag, *leftPayloadIt, *rightPayloadIt);
    }
    KeyIterator midKeyIt = keyBegin + halfSize;
    PayloadIterator midPayloadIt = payloadBegin + halfSize;
    BitonicMergePow2SepPayload(keyBegin, midKeyIt, payloadBegin, dire);
    BitonicMergeSepPayload(midKeyIt, keyEnd, midPayloadIt, dire);
  }
}

template <class KeyIterator, class PayloadIterator>
void BitonicSortSepPayload(KeyIterator keyBegin, KeyIterator keyEnd,
                           PayloadIterator payloadBegin, bool dire) {
  size_t size = keyEnd - keyBegin;
  if (size > 1) {
    size_t halfSize = size / 2;
    KeyIterator keyMid = keyBegin + halfSize;
    PayloadIterator payloadMid = payloadBegin + halfSize;
    BitonicSortSepPayload(keyBegin, keyMid, payloadBegin, !dire);
    BitonicSortSepPayload(keyMid, keyEnd, payloadMid, dire);
    BitonicMergeSepPayload(keyBegin, keyEnd, payloadBegin, dire);
  }
}

template <class KeyIterator, class PayloadIterator>
void BitonicSortSepPayload(KeyIterator keyBegin, KeyIterator keyEnd,
                           PayloadIterator payloadBegin) {
  BitonicSortSepPayload(keyBegin, keyEnd, payloadBegin, true);
}

template <class KeyIterator, class Swap>
void BitonicMergePow2CustomSwap(KeyIterator keyBegin, KeyIterator keyEnd,
                                const Swap& swap, uint64_t beginIdx,
                                bool dire) {
  size_t size = keyEnd - keyBegin;
  if (size > 1) {
    size_t halfSize = size / 2;
    KeyIterator leftKeyIt = keyBegin;
    KeyIterator rightKeyIt = keyBegin + halfSize;
    uint64_t leftBeginIdx = beginIdx;
    for (size_t i = 0; i < halfSize;
         ++i, ++leftKeyIt, ++rightKeyIt, ++leftBeginIdx) {
      bool swapFlag = dire != (*leftKeyIt < *rightKeyIt);
      condSwap(swapFlag, *leftKeyIt, *rightKeyIt);
      swap(swapFlag, leftBeginIdx, leftBeginIdx + halfSize);
    }
    BitonicMergePow2CustomSwap(keyBegin, leftKeyIt, swap, beginIdx, dire);
    BitonicMergePow2CustomSwap(leftKeyIt, keyEnd, swap, leftBeginIdx, dire);
  }
}

template <class KeyIterator, class Swap>
void BitonicMergeCustomSwap(KeyIterator keyBegin, KeyIterator keyEnd,
                            const Swap& swap, uint64_t beginIdx, bool dire) {
  size_t size = keyEnd - keyBegin;
  if (size > 1) {
    size_t halfSize = GetNextPowerOfTwo(size) / 2;
    KeyIterator leftKeyIt = keyBegin;
    KeyIterator rightKeyIt = keyBegin + halfSize;
    uint64_t leftBeginIdx = beginIdx;
    for (size_t i = 0; i < size - halfSize;
         ++i, ++leftKeyIt, ++rightKeyIt, ++leftBeginIdx) {
      bool swapFlag = dire != (*leftKeyIt < *rightKeyIt);
      condSwap(swapFlag, *leftKeyIt, *rightKeyIt);
      swap(swapFlag, leftBeginIdx, leftBeginIdx + halfSize);
    }
    KeyIterator midKeyIt = keyBegin + halfSize;
    uint64_t midIdx = beginIdx + halfSize;
    BitonicMergePow2CustomSwap(keyBegin, midKeyIt, swap, beginIdx, dire);
    BitonicMergeCustomSwap(midKeyIt, keyEnd, swap, midIdx, dire);
  }
}

template <class KeyIterator, class Swap>
void BitonicSortCustomSwap(KeyIterator keyBegin, KeyIterator keyEnd,
                           const Swap& swap, uint64_t beginIdx, bool dire) {
  size_t size = keyEnd - keyBegin;
  if (size > 1) {
    size_t halfSize = size / 2;
    KeyIterator keyMid = keyBegin + halfSize;
    BitonicSortCustomSwap(keyBegin, keyMid, swap, beginIdx, !dire);
    BitonicSortCustomSwap(keyMid, keyEnd, swap, beginIdx + halfSize, dire);
    BitonicMergeCustomSwap(keyBegin, keyEnd, swap, beginIdx, dire);
  }
}

template <class KeyIterator, class Swap>
void BitonicSortCustomSwap(KeyIterator keyBegin, KeyIterator keyEnd,
                           const Swap& swap) {
  BitonicSortCustomSwap(keyBegin, keyEnd, swap, 0, true);
}

template <class KeyIterator, class PayloadIterator>
void ParBitonicMergePow2SepPayload(KeyIterator keyBegin, KeyIterator keyEnd,
                                   PayloadIterator payloadBegin, int numThreads,
                                   bool dire) {
  size_t size = keyEnd - keyBegin;
  if (size > 1) {
    size_t halfSize = size / 2;
    KeyIterator midIt = keyBegin + halfSize;
    PayloadIterator midPayloadIt = payloadBegin + halfSize;
    // int chunkSize = divRoundUp(halfSize, numThreads);
    // #pragma omp for schedule(static, chunkSize)
    for (size_t i = 0; i < halfSize; ++i) {
      KeyIterator leftKeyIt = keyBegin + i;
      KeyIterator rightKeyIt = midIt + i;
      PayloadIterator leftPayloadIt = payloadBegin + i;
      PayloadIterator rightPayloadIt = midPayloadIt + i;
      bool swapFlag = dire != (*leftKeyIt < *rightKeyIt);
      condSwap(swapFlag, *leftKeyIt, *rightKeyIt);
      condSwap(swapFlag, *leftPayloadIt, *rightPayloadIt);
    }
    int leftThreads = numThreads / 2;
    int rightThreads = numThreads - leftThreads;
    if (leftThreads > 1) {
      ParBitonicMergePow2SepPayload(keyBegin, midIt, payloadBegin, leftThreads,
                                    dire);
    } else {
      BitonicMergePow2SepPayload(keyBegin, midIt, payloadBegin, dire);
    }
    if (rightThreads > 1) {
      ParBitonicMergePow2SepPayload(midIt, keyEnd, midPayloadIt, rightThreads,
                                    dire);
    } else {
      BitonicMergePow2SepPayload(midIt, keyEnd, midPayloadIt, dire);
    }
  }
}

template <class KeyIterator, class PayloadIterator>
void ParBitonicMergeSepPayload(KeyIterator keyBegin, KeyIterator keyEnd,
                               PayloadIterator payloadBegin, int numThreads,
                               bool dire) {
  Assert(numThreads > 1);
  size_t size = keyEnd - keyBegin;
  if (size > 1) {
    size_t halfSize = GetNextPowerOfTwo(size) / 2;
    KeyIterator midIt = keyBegin + halfSize;
    PayloadIterator midPayloadIt = payloadBegin + halfSize;
    size_t loopCycle = size - halfSize;
    // int chunkSize = divRoundUp(loopCycle, numThreads);
    // #pragma omp for schedule(static, chunkSize)
    for (size_t i = 0; i < loopCycle; ++i) {
      KeyIterator leftKeyIt = keyBegin + i;
      KeyIterator rightKeyIt = midIt + i;
      PayloadIterator leftPayloadIt = payloadBegin + i;
      PayloadIterator rightPayloadIt = midPayloadIt + i;
      bool swapFlag = dire != (*leftKeyIt < *rightKeyIt);
      condSwap(swapFlag, *leftKeyIt, *rightKeyIt);
      condSwap(swapFlag, *leftPayloadIt, *rightPayloadIt);
    }
    int leftThreads = numThreads * halfSize / size;
    int rightThreads = numThreads - leftThreads;
    if (rightThreads > 1) {
#pragma omp task
      {
        ParBitonicMergePow2SepPayload(keyBegin, midIt, payloadBegin,
                                      leftThreads, dire);
      }
#pragma omp task
      {
        ParBitonicMergeSepPayload(midIt, keyEnd, midPayloadIt, rightThreads,
                                  dire);
      }
#pragma omp taskwait
    } else {
      BitonicMergePow2SepPayload(keyBegin, midIt, payloadBegin, dire);
      BitonicMergeSepPayload(midIt, keyEnd, midPayloadIt, dire);
    }
  }
}

INLINE static bool parity(int x) {
  x ^= x >> 16;
  x ^= x >> 8;
  x ^= x >> 4;
  x ^= x >> 2;
  x ^= x >> 1;
  return x & 1;
}

template <class Iterator>
void GetThreadSize(Iterator begin, Iterator end, uint64_t size) {
  int numThreads = end - begin;
  if (numThreads == 0) {
    return;
  }
  if (numThreads == 1) {
    *begin = size;
    return;
  }
  size_t halfSize = size / 2;
  Iterator mid = begin + numThreads / 2;
  GetThreadSize(begin, mid, halfSize);
  GetThreadSize(mid, end, size - halfSize);
}

template <class KeyIterator, class PayloadIterator>
void ParBitonicSortSepPayloadRecursiveMerge(KeyIterator keyBegin, KeyIterator keyEnd,
                              PayloadIterator payloadBegin, int numThreads,
                              bool dire) {
  uint64_t size = keyEnd - keyBegin;
  if (numThreads > 1 && size > 1) {
    size_t halfSize = size / 2;
    KeyIterator mid = keyBegin + halfSize;
    PayloadIterator midPayload = payloadBegin + halfSize;
    int leftThreads = numThreads / 2;
    int rightThreads = numThreads - leftThreads;
    ParBitonicSortSepPayloadRecursiveMerge(keyBegin, mid, payloadBegin, leftThreads, !dire);
    ParBitonicSortSepPayloadRecursiveMerge(mid, keyEnd, midPayload, rightThreads, dire);
    BitonicMergeSepPayload(keyBegin, keyEnd, payloadBegin, dire);
  }
}

template <class KeyIterator, class PayloadIterator>
void ParBitonicSortSepPayload(KeyIterator keyBegin, KeyIterator keyEnd,
                              PayloadIterator payloadBegin, int numThreads,
                              bool dire) {
  size_t size = keyEnd - keyBegin;
  size_t maxThread = size * (sizeof(*keyBegin) + sizeof(*payloadBegin)) >> 18;
  numThreads = std::min((uint64_t)numThreads, maxThread);
  int logNumThreads = GetLogBaseTwo(numThreads);
  numThreads = 1UL << logNumThreads;
  if (numThreads > 1) {
    std::vector<uint64_t> threadSize(numThreads);
    GetThreadSize(threadSize.begin(), threadSize.end(), size);
    // prefix sum
    for (int i = 1; i < numThreads; ++i) {
      threadSize[i] += threadSize[i - 1];
    }
    #pragma omp parallel for num_threads(numThreads)
    for (int i = 0; i < numThreads; ++i) {
      size_t leftOffset = i == 0 ? 0 : threadSize[i - 1];
      size_t rightOffset = threadSize[i];
      KeyIterator threadKeyBegin = keyBegin + leftOffset;
      KeyIterator threadKeyEnd = keyBegin + rightOffset;
      PayloadIterator threadPayloadBegin = payloadBegin + leftOffset;
      BitonicSortSepPayload(threadKeyBegin, threadKeyEnd, threadPayloadBegin, dire != parity(i) ^ (logNumThreads & 1));
    }
    ParBitonicSortSepPayloadRecursiveMerge(keyBegin, keyEnd, payloadBegin, numThreads, dire);
  } else {
    BitonicSortSepPayload(keyBegin, keyEnd, payloadBegin, dire);
  }
}

template <class KeyIterator, class PayloadIterator>
void ParBitonicSortSepPayload(KeyIterator keyBegin, KeyIterator keyEnd,
                              PayloadIterator payloadBegin, int numThreads) {
  ParBitonicSortSepPayload(keyBegin, keyEnd, payloadBegin, numThreads, true);
}

template <class KeyIterator, class Swap>
void ParBitonicMergePow2CustomSwap(KeyIterator keyBegin, KeyIterator keyEnd,
                                   const Swap& swap, uint64_t beginIdx,
                                   int numThreads, bool dire) {
  size_t size = keyEnd - keyBegin;
  if (size > 1) {
    size_t halfSize = size / 2;
    KeyIterator midIt = keyBegin + halfSize;
    uint64_t midIdx = beginIdx + halfSize;
    for (size_t i = 0; i < halfSize; ++i) {
      KeyIterator leftKeyIt = keyBegin + i;
      KeyIterator rightKeyIt = midIt + i;
      uint64_t leftIdx = beginIdx + i;
      uint64_t rightIdx = midIdx + i;
      bool swapFlag = dire != (*leftKeyIt < *rightKeyIt);
      condSwap(swapFlag, *leftKeyIt, *rightKeyIt);
      swap(swapFlag, leftIdx, rightIdx);
    }
    int leftThreads = numThreads / 2;
    int rightThreads = numThreads - leftThreads;
    if (leftThreads > 1) {
      ParBitonicMergePow2CustomSwap(keyBegin, midIt, swap, beginIdx,
                                    leftThreads, dire);
    } else {
      BitonicMergePow2CustomSwap(keyBegin, midIt, swap, beginIdx, dire);
    }
    if (rightThreads > 1) {
      ParBitonicMergePow2CustomSwap(midIt, keyEnd, swap, midIdx, rightThreads,
                                    dire);
    } else {
      BitonicMergePow2CustomSwap(midIt, keyEnd, swap, midIdx, dire);
    }
  }
}

template <class KeyIterator, class Swap>
void ParBitonicMergeCustomSwap(KeyIterator keyBegin, KeyIterator keyEnd,
                               const Swap& swap, uint64_t beginIdx,
                               int numThreads, bool dire) {
  Assert(numThreads > 1);
  size_t size = keyEnd - keyBegin;
  if (size > 1) {
    size_t halfSize = GetNextPowerOfTwo(size) / 2;
    KeyIterator midIt = keyBegin + halfSize;
    uint64_t midIdx = beginIdx + halfSize;
    size_t loopCycle = size - halfSize;
    // int chunkSize = divRoundUp(loopCycle, numThreads);
    // #pragma omp for schedule(static, chunkSize)
    for (size_t i = 0; i < loopCycle; ++i) {
      KeyIterator leftKeyIt = keyBegin + i;
      KeyIterator rightKeyIt = midIt + i;
      uint64_t leftIdx = beginIdx + i;
      uint64_t rightIdx = midIdx + i;
      bool swapFlag = dire != (*leftKeyIt < *rightKeyIt);
      condSwap(swapFlag, *leftKeyIt, *rightKeyIt);
      swap(swapFlag, leftIdx, rightIdx);
    }
    int leftThreads = numThreads * halfSize / size;
    int rightThreads = numThreads - leftThreads;
    if (rightThreads > 1) {
#pragma omp task
      {
        ParBitonicMergePow2CustomSwap(keyBegin, midIt, swap, beginIdx,
                                      leftThreads, dire);
      }
#pragma omp task
      {
        ParBitonicMergeCustomSwap(midIt, keyEnd, swap, midIdx, rightThreads,
                                  dire);
      }
#pragma omp taskwait
    } else {
      BitonicMergePow2CustomSwap(keyBegin, midIt, swap, beginIdx, dire);
      BitonicMergeCustomSwap(midIt, keyEnd, swap, midIdx, dire);
    }
  }
}

template <class KeyIterator, class Swap>
void ParBitonicSortCustomSwap(KeyIterator keyBegin, KeyIterator keyEnd,
                              const Swap& swap, uint64_t beginIdx,
                              int numThreads, bool dire) {
  if (numThreads > 1) {
    size_t size = keyEnd - keyBegin;
    if (size > 1) {
      size_t halfSize = size / 2;
      KeyIterator keyMid = keyBegin + halfSize;
      uint64_t midIdx = beginIdx + halfSize;
      int leftThreads = numThreads / 2;
      int rightThreads = numThreads - leftThreads;
#pragma omp task
      {
        ParBitonicSortCustomSwap(keyBegin, keyMid, swap, beginIdx, leftThreads,
                                 !dire);
      }
#pragma omp task
      {
        ParBitonicSortCustomSwap(keyMid, keyEnd, swap, midIdx, rightThreads,
                                 dire);
      }
#pragma omp taskwait
      ParBitonicMergeCustomSwap(keyBegin, keyEnd, swap, beginIdx, numThreads,
                                dire);
    }
  } else {
    BitonicSortCustomSwap(keyBegin, keyEnd, swap, beginIdx, dire);
  }
}

template <class KeyIterator, class Swap>
void ParBitonicSortCustomSwap(KeyIterator keyBegin, KeyIterator keyEnd,
                              const Swap& swap, int numThreads) {
  ParBitonicSortCustomSwap(keyBegin, keyEnd, swap, 0, numThreads, true);
}

}  // namespace EM::Algorithm