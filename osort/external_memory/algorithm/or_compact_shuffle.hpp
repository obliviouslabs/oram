#pragma once
#include "external_memory/extemvector.hpp"

/// This file contains the implementation of the OrCompact and OrShuffle.

namespace EM::Algorithm {
template <class Iterator, class MarkIterator>
void OrOffCompactSeparateMark(Iterator begin, Iterator end,
                              const MarkIterator markBegin, size_t z) {
  size_t n = end - begin;
  Iterator mid = begin + n / 2;
  const MarkIterator markMid = markBegin + n / 2;
  size_t valMarkMid = *markMid;
  size_t m = valMarkMid - *markBegin;
  if (n == 2) {
    const MarkIterator markEnd = markMid + 1;
    bool swapFlag = (((!m) & (*markEnd - valMarkMid)) != z);
    condSwap(swapFlag, *begin, *mid);
    return;
  }
  size_t zLeft = z % (n / 2);
  size_t zRight = (z + m) % (n / 2);
  OrOffCompactSeparateMark(begin, mid, markBegin, zLeft);
  OrOffCompactSeparateMark(mid, end, markMid, zRight);
  bool s = ((zLeft + m >= n / 2) != (z >= n / 2));
  Iterator leftIt = begin, rightIt = mid;
  for (size_t i = 0; i != n / 2; ++i, ++leftIt, ++rightIt) {
    bool b = (s != (i >= zRight));
    condSwap(b, *leftIt, *rightIt);
  }
}

template <class Iterator, class MarkIterator>
void OrCompactSeparateMark(Iterator begin, Iterator end,
                           const MarkIterator markBegin) {
  size_t n = end - begin;
  if (n <= 1) {
    return;
  }
  size_t n1 = 1UL << GetLogBaseTwo(n);
  size_t n2 = n - n1;
  if (n2 == 0) {
    OrOffCompactSeparateMark(begin, end, markBegin, 0);
    return;
  }
  Iterator n2It = begin + n2;
  const MarkIterator markMid = markBegin + n2;
  size_t m = *markMid - *markBegin;  // prefix sum
  OrCompactSeparateMark(begin, n2It, markBegin);
  OrOffCompactSeparateMark(n2It, end, markMid, (n1 - n2 + m) % n1);
  Iterator leftIt = begin;
  Iterator rightIt = begin + n1;
  for (size_t i = 0; i != n2; ++i, ++leftIt, ++rightIt) {
    condSwap(i >= m, *leftIt, *rightIt);
  }
}

template <class Iterator, class MarkIterator>
void OrOffDistributeSeparateMark(Iterator begin, Iterator end,
                                 const MarkIterator markBegin, size_t z) {
  size_t n = end - begin;
  Iterator mid = begin + n / 2;
  const MarkIterator markMid = markBegin + n / 2;
  size_t valMarkMid = *markMid;
  size_t m = valMarkMid - *markBegin;
  if (n == 2) {
    const MarkIterator markEnd = markMid + 1;
    bool swapFlag = (((!m) & (*markEnd - valMarkMid)) != z);
    condSwap(swapFlag, *begin, *mid);
    return;
  }
  size_t zLeft = z % (n / 2);
  size_t zRight = (z + m) % (n / 2);
  bool s = ((zLeft + m >= n / 2) != (z >= n / 2));
  Iterator leftIt = begin, rightIt = mid;
  for (size_t i = 0; i != n / 2; ++i, ++leftIt, ++rightIt) {
    bool b = (s != (i >= zRight));
    condSwap(b, *leftIt, *rightIt);
  }
  OrOffDistributeSeparateMark(begin, mid, markBegin, zLeft);
  OrOffDistributeSeparateMark(mid, end, markMid, zRight);
}

template <class Iterator, class MarkIterator>
void OrDistributeSeparateMark(Iterator begin, Iterator end,
                              const MarkIterator markBegin) {
  size_t n = end - begin;
  if (n <= 1) {
    return;
  }
  size_t n1 = 1UL << GetLogBaseTwo(n);
  size_t n2 = n - n1;
  if (n2 == 0) {
    OrOffDistributeSeparateMark(begin, end, markBegin, 0);
    return;
  }
  Iterator n2It = begin + n2;
  const MarkIterator markMid = markBegin + n2;
  size_t m = *markMid - *markBegin;  // prefix sum
  Iterator leftIt = begin;
  Iterator rightIt = begin + n1;
  for (size_t i = 0; i != n2; ++i, ++leftIt, ++rightIt) {
    condSwap(i >= m, *leftIt, *rightIt);
  }
  OrDistributeSeparateMark(begin, n2It, markBegin);
  OrOffDistributeSeparateMark(n2It, end, markMid, (n1 - n2 + m) % n1);
}

template <typename MarkType = uint32_t, class Iterator, class Check>
void OrCompact(Iterator begin, Iterator end, const Check& isMarked) {
  size_t n = end - begin;
  Assert(n < UINT32_MAX);
  if (n <= 1) {
    return;
  }
  std::vector<MarkType> markPrefix(n + 1);
  MarkType prefixSum = 0;
  auto markIt = markPrefix.begin();
  *(markIt++) = 0;
  for (auto it = begin; it != end; ++it) {
    prefixSum += !!isMarked(*it);
    *(markIt++) = prefixSum;
  }
  OrCompactSeparateMark(begin, end, markPrefix.begin());
}

template <class Iterator, typename T>
void initIteratorVal(Iterator& it, const T& val) {
  if constexpr (std::is_same<Iterator,
                             typename std::vector<T>::iterator>::value) {
    *it = val;
  } else {
    it.derefWriteOnly() = val;
  }
}

template <class Iterator, class MarkIterator>
void OrShuffle(Iterator begin, Iterator end, MarkIterator markBegin) {
  size_t n = end - begin;
  if (n <= 1) {
    return;
  }
  if (n == 2) {
    bool flag = UniformRandomBit();
    condSwap(flag, *begin, *(begin + 1));
    return;
  }
  size_t halfSize = n / 2;
  size_t k = halfSize;

  auto it = markBegin;
  size_t currsum = 0;
  initIteratorVal(it, currsum);
  ++it;
  for (size_t N = n; N > 0; --N) {
    bool chooseFlag = UniformRandom(0, N - 1) < k;
    CMOV(chooseFlag, k, k - 1);
    CMOV(chooseFlag, currsum, currsum + 1);  // prefix sum
    initIteratorVal(it, currsum);
    ++it;
  }
  OrCompactSeparateMark(begin, end, markBegin);
  OrShuffle(begin, begin + halfSize, markBegin);
  OrShuffle(begin + halfSize, end, markBegin);
}

template <class Iterator>
void OrShuffle(Iterator begin, Iterator end) {
  size_t n = end - begin;
  if (n <= 1) {
    return;
  }
  using T = typename std::iterator_traits<Iterator>::value_type;
  if constexpr (std::is_same<Iterator,
                             typename StdVector<T>::Iterator>::value) {
    std::vector<size_t> marks(n + 1);
    OrShuffle(begin, end, marks.begin());
  } else if constexpr (std::is_same<Iterator,
                                    typename std::vector<T>::iterator>::value) {
    std::vector<size_t> marks(n + 1);
    OrShuffle(begin, end, marks.begin());
  } else {
    EM::ExtVector::Vector<size_t, 4096, true, true,
                          (uint64_t)ENCLAVE_SIZE * 16UL>
        marks(n + 1);
    OrShuffle(begin, end, marks.begin());
  }
}

template <typename Vec>
void OrShuffle(Vec& vec) {
  OrShuffle(vec.begin(), vec.end());
}

// implements https://dl.acm.org/doi/pdf/10.1145/1989493.1989555
template <class Iterator, class Check>
void GoodrichCompact(Iterator begin, Iterator end, const Check& isMarked) {
  size_t n = end - begin;
  std::vector<uint64_t> dist(n);
  uint64_t numMarked = 0;
  for (size_t i = 0; i < n; ++i) {
    dist[i] = i - numMarked;
    CMOV(isMarked(*(begin + i)), numMarked, numMarked + 1);
  }
  uint64_t level = GetLogBaseTwo(GetNextPowerOfTwo(n));
  for (uint64_t i = 0; i < level; ++i) {
    for (uint64_t j = (1UL << i); j < n; ++j) {
      bool jMarked = isMarked(*(begin + j));
      bool swapFlag = jMarked & (!!(dist[j] & ((2UL << i) - 1)));
      condSwap(swapFlag, *(begin + j), *(begin + j - (1UL << i)));
      CMOV(swapFlag, dist[j - (1UL << i)], (dist[j] >> (i + 1) << (i + 1)));
    }
  }
}
}  // namespace EM::Algorithm