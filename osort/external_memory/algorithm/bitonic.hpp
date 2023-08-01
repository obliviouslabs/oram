#pragma once
#include <vector>

#include "external_memory/extemvector.hpp"
#include "sort_def.hpp"

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
                             typename std::vector<T>::iterator>::value) {
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
    using IOVector = typename
        std::remove_reference<decltype(*(Iterator::getNullVector()))>::type;
    constexpr size_t CachePageSize =
        IOVector::item_per_page * sizeof(TaggedT<T>);
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

}  // namespace EM::Algorithm