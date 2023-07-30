#pragma once
#include "external_memory/noncachedvector.hpp"
#include "sort_building_blocks.hpp"
/// This file implements the multi-way bucket osort
/// (https://arxiv.org/pdf/2008.00332.pdf) using single core.
/// SPMS sort replaced with external-memory merge sort
// NOTE dummy flag is set at the least significant bit in this example
namespace EM::Algorithm {
using EM::NonCachedVector::Vector;
template <typename T, class Iterator>
void binPacking(Iterator begin, Iterator end, Iterator outputBegin, size_t Z,
                uint64_t bitMask) {
  size_t size = end - begin;
  size_t numBucket = size / Z;
  size_t logNumBucket = GetLogBaseTwo(numBucket);
  uint64_t baseMask = (bitMask >> logNumBucket) << 1;
  uint64_t rangeMask = baseMask * (numBucket - 1);
  uint64_t rangeMaskWithDummyFlag = rangeMask + 1UL;
  std::vector<TaggedT<T>> temp(size * 2);
  std::copy(begin, end, temp.begin());
  auto iter = temp.begin() + size;
  for (size_t i = 0; i < numBucket; ++i) {
    for (size_t j = 0; j < Z; ++j) {
      // dummyFlag at the end
      iter->tag = i * baseMask + 1UL;
      ++iter;
    }
  }
  const auto cmpTag = [rangeMaskWithDummyFlag = rangeMaskWithDummyFlag](
                          const auto& element1, const auto& element2) {
    return (element1.tag & rangeMaskWithDummyFlag) <
           (element2.tag & rangeMaskWithDummyFlag);
  };
  BitonicSort(temp, cmpTag);
  size_t count = 0;
  size_t curr = 0;
  for (auto& element : temp) {
    size_t newCurr = (element.tag & rangeMask);
    CMOV(newCurr != curr, count, 0UL);
    curr = newCurr;
    // excess real elements
    Assert(!(count >= Z && !(element.tag & 1UL)));
    // mark excessive element with -1
    CMOV(count >= Z, element.tag, -1UL);
    ++count;
  }
  BitonicSort(temp, cmpTag);  // could be replaced by OrCompact
  std::copy(temp.begin(), temp.begin() + size, outputBegin);
}

template <typename T, class Iterator>
void butterflyRouting(Iterator begin, Iterator end, Iterator outputBegin,
                      size_t Z, uint64_t bitMask, size_t gamma) {
  size_t size = end - begin;
  size_t numBucket = size / Z;
  if (numBucket <= gamma) {
    binPacking<T>(begin, end, outputBegin, Z, bitMask);
    return;
  }
  size_t num_layer = GetLogBaseTwo(numBucket);
  size_t half_num_layer = num_layer / 2;
  size_t numBucketPerBatch = 1UL << half_num_layer;
  size_t batchSize = Z * numBucketPerBatch;
  size_t batchCount = size / batchSize;
  std::vector<TaggedT<T>> nextLayer(size);
  for (size_t batchIdx = 0; batchIdx < batchCount; ++batchIdx) {
    auto batchBegin = begin + batchIdx * batchSize;
    auto batchEnd = begin + (batchIdx + 1) * batchSize;
    butterflyRouting<T>(batchBegin, batchEnd, batchBegin, Z, bitMask, gamma);
    for (size_t bucketIdx = 0; bucketIdx < numBucketPerBatch; ++bucketIdx) {
      std::copy(batchBegin + bucketIdx * Z, batchBegin + (bucketIdx + 1) * Z,
                nextLayer.begin() + bucketIdx * size / numBucketPerBatch +
                    batchIdx * Z);
    }
  }
  size_t nextLayerBucketPerBatch = 1UL << (num_layer - half_num_layer);
  size_t nextLayerBatchSize = nextLayerBucketPerBatch * Z;
  size_t nextLayerBatchCount = size / nextLayerBatchSize;
  for (size_t batchIdx = 0; batchIdx < nextLayerBatchCount; ++batchIdx) {
    butterflyRouting<T>(nextLayer.begin() + batchIdx * nextLayerBatchSize,
                        nextLayer.begin() + (batchIdx + 1) * nextLayerBatchSize,
                        outputBegin + batchIdx * nextLayerBatchSize, Z,
                        bitMask >> half_num_layer, gamma);
  }
}

template <typename IOIterator,
          typename T = typename std::iterator_traits<IOIterator>::value_type,
          typename IOVector = std::remove_reference<
              decltype(*(IOIterator::getNullVector()))>::type>
Vector<TaggedT<T>, 4096> tagAndPad(IOIterator begin, IOIterator end,
                                   uint64_t Z) {
  size_t inputSize = end - begin;
  uint64_t intermidiateSize = std::max(GetNextPowerOfTwo(2 * (end - begin)), Z);
  Vector<TaggedT<T>, 4096> tv(intermidiateSize);
  uint64_t numBucket = intermidiateSize / Z;
  uint64_t numRealPerBucket = divRoundUp(inputSize, numBucket);
  typename IOVector::PrefetchReader inputReader(begin, end);
  typename Vector<TaggedT<T>, 4096>::Writer taggedTWriter(tv.begin(), tv.end());

  for (uint64_t bucketIdx = 0; bucketIdx < numBucket; ++bucketIdx) {
    for (uint64_t offset = 0; offset < Z; ++offset) {
      TaggedT<T> tt;
      tt.tag = UniformRandom() &
               -2UL;  // dummy flag is set at the least significant bit
      if (offset < numRealPerBucket && !inputReader.eof()) {
        tt.v = inputReader.read();
      } else {
        tt.tag |= 1UL;
      }
      taggedTWriter.write(tt);
    }
  }
  taggedTWriter.flush();
  return tv;
}

template <typename T>
void externalButterflyRouting(
    typename Vector<TaggedT<T>, 4096>::Iterator begin,
    typename Vector<TaggedT<T>, 4096>::Iterator end,
    typename Vector<TaggedT<T>, 4096>::Iterator outputBegin, size_t Z,
    uint64_t bitMask, uint64_t heapSize = DEFAULT_HEAP_SIZE) {
  size_t size = end - begin;
  const size_t gamma = GetLogBaseTwo(size);
  if ((size * 2 + gamma * Z) * sizeof(TaggedT<T>) <= heapSize) {
    std::vector<TaggedT<T>> mem(size);
    CopyIn(begin, end, mem.begin());
    butterflyRouting<T>(mem.begin(), mem.end(), mem.begin(), Z, bitMask, gamma);
    CopyOut(mem.begin(), mem.end(), outputBegin);
    return;
  }
  size_t num_layer = GetLogBaseTwo(size / Z);
  size_t half_num_layer = num_layer / 2;
  size_t numBucketPerBatch = 1UL << half_num_layer;
  size_t batchSize = Z * numBucketPerBatch;
  size_t batchCount = size / batchSize;
  Vector<TaggedT<T>, 4096> nextLayer(size);
  std::vector<TaggedT<T>> temp(Z);
  for (size_t batchIdx = 0; batchIdx < batchCount; ++batchIdx) {
    auto batchBegin = begin + batchIdx * batchSize;
    auto batchEnd = begin + (batchIdx + 1) * batchSize;
    externalButterflyRouting<T>(batchBegin, batchEnd, batchBegin, Z, bitMask,
                                heapSize);
    for (size_t bucketIdx = 0; bucketIdx < numBucketPerBatch; ++bucketIdx) {
      CopyIn(batchBegin + bucketIdx * Z, batchBegin + (bucketIdx + 1) * Z,
             temp.begin());
      CopyOut(temp.begin(), temp.end(),
              nextLayer.begin() + bucketIdx * size / numBucketPerBatch +
                  batchIdx * Z);
    }
  }
  size_t nextLayerBucketPerBatch = 1UL << (num_layer - half_num_layer);
  size_t nextLayerBatchSize = nextLayerBucketPerBatch * Z;
  size_t nextLayerBatchCount = size / nextLayerBatchSize;
  for (size_t batchIdx = 0; batchIdx < nextLayerBatchCount; ++batchIdx) {
    externalButterflyRouting<T>(
        nextLayer.begin() + batchIdx * nextLayerBatchSize,
        nextLayer.begin() + (batchIdx + 1) * nextLayerBatchSize,
        outputBegin + batchIdx * nextLayerBatchSize, Z,
        bitMask >> half_num_layer, heapSize);
  }
}

template <class IOIterator>
void CABucketShuffle(IOIterator begin, IOIterator end,
                     uint64_t heapSize = DEFAULT_HEAP_SIZE) {
  using T = typename std::iterator_traits<IOIterator>::value_type;
  using IOVector =
      std::remove_reference<decltype(*(IOIterator::getNullVector()))>::type;
  const uint64_t Z = 512;
  Vector<TaggedT<T>, 4096> tv = tagAndPad(begin, end, Z);
  externalButterflyRouting<T>(tv.begin(), tv.end(), tv.begin(), Z, 1UL << 63,
                              heapSize);
  typename IOVector::Writer outputWriter(begin, end, 1);
  std::vector<TaggedT<T>> temp(Z);
  uint64_t prev = 0;
  for (size_t bucketIdx = 0; bucketIdx < tv.size() / Z; ++bucketIdx) {
    CopyIn(tv.begin() + bucketIdx * Z, tv.begin() + (bucketIdx + 1) * Z,
           temp.begin());
    BitonicSort(
        temp.begin(), temp.end(),
        [](const TaggedT<T>& a, const TaggedT<T>& b) { return a.tag < b.tag; });
    for (const TaggedT<T>& element : temp) {
      if (!(element.tag & 1UL)) {
        outputWriter.write(element.v);
        Assert(element.tag >= prev && (prev = element.tag || true));
      }
    }
  }
  outputWriter.flush();
}

template <class IOIterator>
void CABucketSort(IOIterator begin, IOIterator end,
                  uint64_t heapSize = DEFAULT_HEAP_SIZE) {
  CABucketShuffle(begin, end, heapSize);
  ExtMergeSort<false>(begin, end, heapSize, 1);
  // baseline algorithm, does not update counter for authentication here
}

// due to other overhead, cannot utilize the full capacity of the heap
template <typename Vec>
void CABucketSort(Vec& vec, uint64_t heapSize = DEFAULT_HEAP_SIZE * 14 / 15) {
  CABucketSort(vec.begin(), vec.end(), heapSize);
}

template <typename Vec>
void CABucketShuffle(Vec& vec,
                     uint64_t heapSize = DEFAULT_HEAP_SIZE * 14 / 15) {
  CABucketShuffle(vec.begin(), vec.end(), heapSize);
}
}  // namespace EM::Algorithm