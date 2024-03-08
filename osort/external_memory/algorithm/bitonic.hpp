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
    KeyIterator leftKeyIt = keyBegin;
    KeyIterator rightKeyIt = keyBegin + halfSize;
    PayloadIterator leftPayloadIt = payloadBegin;
    PayloadIterator rightPayloadIt = payloadBegin + halfSize;
#pragma omp parallel for num_threads(numThreads)
    for (size_t i = 0; i < halfSize; ++i) {
      bool swapFlag = dire != (*leftKeyIt < *rightKeyIt);
      condSwap(swapFlag, *leftKeyIt, *rightKeyIt);
      condSwap(swapFlag, *leftPayloadIt, *rightPayloadIt);
    }
    int leftThreads = numThreads / 2;
    int rightThreads = numThreads - leftThreads;
    if (leftThreads > 1) {
      ParBitonicMergePow2SepPayload(keyBegin, leftKeyIt, payloadBegin,
                                    leftThreads, dire);
    } else {
      BitonicMergePow2SepPayload(keyBegin, leftKeyIt, payloadBegin, dire);
    }
    if (rightThreads > 1) {
      ParBitonicMergePow2SepPayload(leftKeyIt, keyEnd, leftPayloadIt,
                                    rightThreads, dire);
    } else {
      BitonicMergePow2SepPayload(leftKeyIt, keyEnd, leftPayloadIt, dire);
    }
    ++leftKeyIt, ++rightKeyIt, ++leftPayloadIt, ++rightPayloadIt;
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
    KeyIterator leftKeyIt = keyBegin;
    KeyIterator rightKeyIt = keyBegin + halfSize;
    PayloadIterator leftPayloadIt = payloadBegin;
    PayloadIterator rightPayloadIt = payloadBegin + halfSize;
#pragma omp parallel for num_threads(numThreads)
    for (size_t i = 0; i < size - halfSize; ++i) {
      bool swapFlag = dire != (*leftKeyIt < *rightKeyIt);
      condSwap(swapFlag, *leftKeyIt, *rightKeyIt);
      condSwap(swapFlag, *leftPayloadIt, *rightPayloadIt);
    }
    KeyIterator midKeyIt = keyBegin + halfSize;
    PayloadIterator midPayloadIt = payloadBegin + halfSize;
    int leftThreads = numThreads * halfSize / size;
    int rightThreads = numThreads - leftThreads;
    if (rightThreads > 1) {
#pragma omp task
      {
        ParBitonicMergePow2SepPayload(keyBegin, midKeyIt, payloadBegin,
                                      leftThreads, dire);
      }
#pragma omp task
      {
        ParBitonicMergeSepPayload(midKeyIt, keyEnd, midPayloadIt, rightThreads,
                                  dire);
      }
#pragma omp taskwait
    } else {
      BitonicMergePow2SepPayload(keyBegin, midKeyIt, payloadBegin, dire);
      BitonicMergeSepPayload(midKeyIt, keyEnd, midPayloadIt, dire);
    }
    ++leftKeyIt, ++rightKeyIt, ++leftPayloadIt, ++rightPayloadIt;
  }
}

template <class KeyIterator, class PayloadIterator>
void ParBitonicSortSepPayload(KeyIterator keyBegin, KeyIterator keyEnd,
                              PayloadIterator payloadBegin, int numThreads,
                              bool dire) {
  if (numThreads > 1) {
    size_t size = keyEnd - keyBegin;
    if (size > 1) {
      size_t halfSize = size / 2;
      KeyIterator keyMid = keyBegin + halfSize;
      PayloadIterator payloadMid = payloadBegin + halfSize;
      int leftThreads = numThreads / 2;
      int rightThreads = numThreads - leftThreads;
#pragma omp parallel
      {
#pragma omp task
          {ParBitonicSortSepPayload(keyBegin, keyMid, payloadBegin, leftThreads,
                                    !dire);
    }
#pragma omp task
    {
      ParBitonicSortSepPayload(keyMid, keyEnd, payloadMid, rightThreads, dire);
    }
#pragma omp taskwait
  }
  ParBitonicMergeSepPayload(keyBegin, keyEnd, payloadBegin, numThreads, dire);
}
}  // namespace EM::Algorithm
else {
  BitonicSortSepPayload(keyBegin, keyEnd, payloadBegin, dire);
}
}

template <class KeyIterator, class PayloadIterator>
void ParBitonicSortSepPayload(KeyIterator keyBegin, KeyIterator keyEnd,
                              PayloadIterator payloadBegin, int numThreads) {
  ParBitonicSortSepPayload(keyBegin, keyEnd, payloadBegin, numThreads, true);
}

}  // namespace EM::Algorithm