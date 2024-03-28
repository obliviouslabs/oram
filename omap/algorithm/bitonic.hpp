#pragma once
#include <vector>

#include "element.hpp"
#include "external_memory/extemvector.hpp"
#include "external_memory/stdvector.hpp"
#include "static_sort.hpp"

namespace Algorithm {

template <class Iterator, typename Compare>
void BitonicMergePow2(Iterator begin, Iterator end, Compare cmp, bool dire) {
  size_t size = end - begin;
  if (size > 1) {
    size_t halfSize = size / 2;
    Iterator leftIt = begin;
    Iterator rightIt = begin + halfSize;
    for (size_t i = 0; i < halfSize; ++i, ++leftIt, ++rightIt) {
      obliSwap(dire != cmp(*leftIt, *rightIt), *leftIt, *rightIt);
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
      obliSwap(dire != cmp(*leftIt, *rightIt), *leftIt, *rightIt);
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
      obliSwap(swapFlag, *leftKeyIt, *rightKeyIt);
      obliSwap(swapFlag, *leftPayloadIt, *rightPayloadIt);
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
      obliSwap(swapFlag, *leftKeyIt, *rightKeyIt);
      obliSwap(swapFlag, *leftPayloadIt, *rightPayloadIt);
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
  int numThreads = (int)(end - begin);
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
void ParBitonicSortSepPayloadRecursiveMerge(KeyIterator keyBegin,
                                            KeyIterator keyEnd,
                                            PayloadIterator payloadBegin,
                                            int numThreads, bool dire) {
  uint64_t size = keyEnd - keyBegin;
  if (numThreads > 1 && size > 1) {
    size_t halfSize = size / 2;
    KeyIterator mid = keyBegin + halfSize;
    PayloadIterator midPayload = payloadBegin + halfSize;
    int leftThreads = numThreads / 2;
    int rightThreads = numThreads - leftThreads;
    ParBitonicSortSepPayloadRecursiveMerge(keyBegin, mid, payloadBegin,
                                           leftThreads, !dire);
    ParBitonicSortSepPayloadRecursiveMerge(mid, keyEnd, midPayload,
                                           rightThreads, dire);
    BitonicMergeSepPayload(keyBegin, keyEnd, payloadBegin, dire);
  }
}

template <class KeyIterator, class PayloadIterator>
void ParBitonicSortSepPayload(KeyIterator keyBegin, KeyIterator keyEnd,
                              PayloadIterator payloadBegin, int numThreads,
                              bool dire) {
  size_t size = keyEnd - keyBegin;
  size_t maxThread = size * (sizeof(*keyBegin) + sizeof(*payloadBegin)) >> 17;
  numThreads = (int)std::min((uint64_t)numThreads, maxThread);
  int logNumThreads = GetLogBaseTwo(numThreads);
  numThreads = 1 << logNumThreads;
  if (numThreads > 1) {
    std::vector<uint64_t> threadSize(numThreads);
    GetThreadSize(threadSize.begin(), threadSize.end(), size);
    // prefix sum
    for (int i = 1; i < numThreads; ++i) {
      threadSize[i] += threadSize[i - 1];
    }
#pragma omp parallel for num_threads(numThreads) schedule(static, 1)
    for (int i = 0; i < numThreads; ++i) {
      size_t leftOffset = i == 0 ? 0 : threadSize[i - 1];
      size_t rightOffset = threadSize[i];
      KeyIterator threadKeyBegin = keyBegin + leftOffset;
      KeyIterator threadKeyEnd = keyBegin + rightOffset;
      PayloadIterator threadPayloadBegin = payloadBegin + leftOffset;
      BitonicSortSepPayload(threadKeyBegin, threadKeyEnd, threadPayloadBegin,
                            (dire != parity(i)) ^ (logNumThreads & 1));
    }
    ParBitonicSortSepPayloadRecursiveMerge(keyBegin, keyEnd, payloadBegin,
                                           numThreads, dire);
  } else {
    BitonicSortSepPayload(keyBegin, keyEnd, payloadBegin, dire);
  }
}

template <class KeyIterator, class PayloadIterator>
void ParBitonicSortSepPayload(KeyIterator keyBegin, KeyIterator keyEnd,
                              PayloadIterator payloadBegin, int numThreads) {
  ParBitonicSortSepPayload(keyBegin, keyEnd, payloadBegin, numThreads, true);
}
}  // namespace Algorithm