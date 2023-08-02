#pragma once
#include "external_memory/dynamicvector.hpp"
#include "external_memory/noncachedvector.hpp"
#include "param_select.hpp"
#include "sort_building_blocks.hpp"

/// This file implements the flex-way butterfly osort and oshuffle algorithms.

namespace EM::Algorithm {

/// @brief KWayButterflySort is an oblivious external memory sorting algorithm.
/// The implementation incurs approximately 2.23NlogN exchanges and
/// 2N/B log_{M/B} N/B page transfers, where N is the size of the input array, B
/// is the page size, and M is the size of the memory.
/// @tparam Iterator only supports NonCachedVector::Iterator
/// @param begin the begin iterator of the input array
/// @param end the end iterator of the input array
/// @param inAuth the authentication counter of the input array
/// @param heapSize the size of available memory in bytes
template <class Iterator>
void KWayButterflySort(Iterator begin, Iterator end, uint32_t inAuth = 0,
                       uint64_t heapSize = DEFAULT_HEAP_SIZE);
template <typename Vec>
void KWayButterflySort(Vec& vec, uint32_t inAuth = 0,
                       uint64_t heapSize = DEFAULT_HEAP_SIZE);

/// @brief KWayButterflyShuffling is an oblivious external memory shuffling
/// algorithm (i.e., random permutation).
/// The implementation incurs approximately 2NlogN exchanges and
/// N/B log_{M/B} N/B page transfers, where N is the size of the input array, B
/// is the page size, and M is the size of the memory.
/// @tparam Iterator only supports NonCachedVector::Iterator
/// @param begin the begin iterator of the input array
/// @param end the end iterator of the input array
/// @param inAuth the authentication counter of the input array
/// @param heapSize the size of available memory in bytes
template <class Iterator>
void KWayButterflyOShuffle(Iterator begin, Iterator end, uint32_t inAuth = 0,
                           uint64_t heapSize = DEFAULT_HEAP_SIZE);
template <typename Vec>
void KWayButterflyOShuffle(Vec& vec, uint32_t inAuth = 0,
                           uint64_t heapSize = DEFAULT_HEAP_SIZE);

/// @brief A manager class for flex-way butterfly o-sort.
/// @tparam T the type of elements to sort
/// @tparam WrappedT the wrapped type of elements to sort
template <typename IOIterator, SortMethod task>
class ButterflySorter {
 private:
  using T = typename std::iterator_traits<IOIterator>::value_type;
  using WrappedT = TaggedT<T>;
  using IOVector = typename
      std::remove_reference<decltype(*(IOIterator::getNullVector()))>::type;
  uint64_t Z;                 // bucket size
  uint64_t numTotalBucket;    // total number of buckets
  uint64_t numRealPerBucket;  // number of real elements per bucket

  uint64_t numBucketFit;   // number of buckets that can fit in the heap
  uint64_t numElementFit;  // number of elements that can fit in the heap

  KWayButterflyParams KWayParams =
      {};  // parameters for flex-way butterfly o-sort

  Vector<T>
      mergeSortFirstLayer;  // the first layer of external-memory merge sort
  std::vector<
      std::pair<typename Vector<T>::Iterator, typename Vector<T>::Iterator>>
      mergeSortRanges;  // pairs of iterators that specifies each sorted range
                        // in the first layer of external-memory merge sort

  // writer for the first layer of external-memory merge sort
  typename Vector<T>::Writer mergeSortFirstLayerWriter;

  typename IOVector::PrefetchReader inputReader;  // input reader
  typename IOVector::Writer outputWriter;         // output writer
  WrappedT* batch;                                // batch for sorting

 public:
  /// @brief Construct a new Butterfly Sorter object
  /// @param inputBeginIt the begin iterator of the input array
  /// @param inputEndIt the end iterator of the input array
  /// @param inAuth the counter of the input array for authentication
  /// @param _heapSize the heap size in bytes
  ButterflySorter(IOIterator inputBeginIt, IOIterator inputEndIt,
                  uint32_t inAuth = 0, uint64_t _heapSize = DEFAULT_HEAP_SIZE)
      : inputReader(inputBeginIt, inputEndIt, inAuth),
        mergeSortFirstLayer(
            task == KWAYBUTTERFLYOSORT ? inputEndIt - inputBeginIt : 0),
        numElementFit(_heapSize / sizeof(WrappedT)) {
    size_t size = inputEndIt - inputBeginIt;
    batch = new WrappedT[numElementFit];
    KWayParams = bestKWayButterflyParams(size, numElementFit, sizeof(T));
    Z = KWayParams.Z;
    Assert(numElementFit > 8 * Z);
    numElementFit -=
        divRoundUp(Z * 8 * (sizeof(WrappedT) + 2), sizeof(WrappedT));
    // reserve space for 8 buckets of elements and marks
    numBucketFit = numElementFit / Z;
    numTotalBucket = KWayParams.totalBucket;
    numRealPerBucket = 1 + (size - 1) / numTotalBucket;

    if constexpr (task == KWAYBUTTERFLYOSORT) {
      mergeSortFirstLayerWriter.init(mergeSortFirstLayer.begin(),
                                     mergeSortFirstLayer.end());
    } else {
      static_assert(task == KWAYBUTTERFLYOSHUFFLE);
      outputWriter.init(inputBeginIt, inputEndIt, inAuth + 1);
    }
  }

  ~ButterflySorter() {
    if (batch != NULL) {
      delete[] batch;
    }
  }

  const auto& getMergeSortBatchRanges() { return mergeSortRanges; }

  /// @brief Base case of flex-way butterfly o-sort when input fits in memory
  /// @tparam Iterator should support random access
  /// @param begin begin iterator of the input array
  /// @param end end iterator of the input array
  /// @param ioLayer current layer of butterfly network by page swap passes
  /// @param innerLayer current layer of butterfly network within the batch
  template <class Iterator>
  void KWayButterflySortBasic(Iterator begin, Iterator end, size_t ioLayer,
                              size_t innerLayer) {
    uint64_t numElement = end - begin;
    uint64_t numBucket = numElement / Z;
    uint64_t way = KWayParams.ways[ioLayer][innerLayer];
    Assert(numElement % Z == 0);
    Assert(numBucket % way == 0);
    uint64_t waySize = numElement / way;
    uint64_t wayBucket = numBucket / way;
    if (innerLayer > 0) {
      for (uint64_t i = 0; i < way; ++i) {
        KWayButterflySortBasic(begin + i * waySize, begin + (i + 1) * waySize,
                               ioLayer, innerLayer - 1);
      }
    } else {
      Assert(numBucket == way);
      if (ioLayer == 0) {
        Assert(waySize == Z);
        // tag and pad input
        for (uint64_t i = 0; i < way; ++i) {
          auto it = begin + i * waySize;
          for (uint64_t offset = 0; offset < Z; ++offset, ++it) {
            if (offset < numRealPerBucket && !inputReader.eof()) {
              it->setData(inputReader.read());
            } else {
              it->setDummy();
            }
          }
        }
      }
    }
    Iterator KWayIts[8];
    uint8_t* marks = new uint8_t[8 * Z];
    for (uint64_t j = 0; j < wayBucket; ++j) {
      for (uint64_t i = 0; i < way; ++i) {
        KWayIts[i] = begin + (i * wayBucket + j) * Z;
      }
      MergeSplitKWay(KWayIts, way, Z, batch + numElementFit, marks);
    }
    delete[] marks;
  }

  template <class Iterator>
  void KWayButterflySort(Iterator begin, Iterator end) {
    KWayButterflySort(begin, end, KWayParams.ways.size() - 1);
  }

  template <class Iterator>
  void KWayButterflySort(Iterator begin, Iterator end, size_t ioLayer) {
    bool isLastLayer = ioLayer == KWayParams.ways.size() - 1;
    size_t numInternalWay = getVecProduct(KWayParams.ways[ioLayer]);
    size_t fetchInterval = 1;
    for (size_t layer = 0; layer < ioLayer; ++layer) {
      fetchInterval *= getVecProduct(KWayParams.ways[layer]);
    }
    size_t size = end - begin;
    Assert(size % numInternalWay == 0);
    size_t subSize = size / numInternalWay;
    if (ioLayer > 0) {
      for (size_t i = 0; i < numInternalWay; ++i) {
        KWayButterflySort(begin + i * subSize, begin + (i + 1) * subSize,
                          ioLayer - 1);
      }
    }
    size_t batchSize = numInternalWay * Z;
    size_t batchPerEnclave = 1;
    if (task == KWAYBUTTERFLYOSORT && isLastLayer) {
      batchPerEnclave = numElementFit / batchSize;
      // maximize the chunksize at the last layer
    }
    size_t batchCount = divRoundUp(size, batchSize);

    for (uint64_t batchIdx = 0; batchIdx < batchCount; ++batchIdx) {
      uint64_t intBatchIdx = batchIdx % batchPerEnclave;
      auto batchBegin = batch + intBatchIdx * batchSize;
      if (ioLayer) {  // fetch from intermediate ext vector
        for (uint64_t bucketIdx = 0; bucketIdx < numInternalWay; ++bucketIdx) {
          auto extBeginIt = begin + (batchIdx + fetchInterval * bucketIdx) * Z;
          auto intBeginIt = batchBegin + bucketIdx * Z;
          CopyIn(extBeginIt, extBeginIt + Z, intBeginIt, ioLayer - 1);
        }
      }
      KWayButterflySortBasic(batchBegin, batchBegin + batchSize, ioLayer,
                             KWayParams.ways[ioLayer].size() - 1);
      if (isLastLayer) {
        // last layer, combine with bitonic sort and output
        const auto cmpTag = [](const auto& a, const auto& b) {
          return a.tag < b.tag;
        };
        for (size_t i = 0; i < numInternalWay; ++i) {
          auto it = batchBegin + i * Z;
          Assert(it + Z <= batch + numElementFit);
          BitonicSort(it, it + Z, cmpTag);
          // for shuffling, output directly
          if constexpr (task == KWAYBUTTERFLYOSHUFFLE) {
            for (auto fromIt = it; fromIt != it + Z; ++fromIt) {
              if (!fromIt->isDummy()) {
                outputWriter.write(fromIt->getData());
              }
            }
          }
        }
        if (task == KWAYBUTTERFLYOSORT && (intBatchIdx == batchPerEnclave - 1 ||
                                           batchIdx == batchCount - 1)) {
          // sort the batch and write to first layer of merge sort
          const auto cmpVal = [](const auto& a, const auto& b) {
            return a.v < b.v;
          };
          size_t batchLocalCount = batchIdx == batchCount - 1
                                       ? batchIdx % batchPerEnclave + 1
                                       : batchPerEnclave;
          Assert(batchSize * batchLocalCount <= numElementFit);
          auto realEnd =
              partitionDummy(batch, batch + batchSize * batchLocalCount);
          // partition dummies to the end
          Assert(realEnd <= batch + numElementFit);
          std::sort(batch, realEnd, cmpVal);
          auto mergeSortReaderBeginIt = mergeSortFirstLayerWriter.it;
          for (auto it = batch; it != realEnd; ++it) {
            mergeSortFirstLayerWriter.write(it->getData());
          }

          mergeSortRanges.emplace_back(mergeSortReaderBeginIt,
                                       mergeSortFirstLayerWriter.it);
        }
      } else {  // not last layer, write to intermediate ext vector
        Assert(batchPerEnclave == 1);
        for (uint64_t bucketIdx = 0; bucketIdx < numInternalWay; ++bucketIdx) {
          auto extBeginIt = begin + (batchIdx + fetchInterval * bucketIdx) * Z;
          auto intBeginIt = batch + bucketIdx * Z;
          CopyOut(intBeginIt, intBeginIt + Z, extBeginIt, ioLayer);
        }
      }
    }

    if (isLastLayer) {
      if constexpr (task == KWAYBUTTERFLYOSORT) {
        mergeSortFirstLayerWriter.flush();
      } else {
        outputWriter.flush();
      }
    }
  }

  void sort() {
    if (batch == NULL) {
      printf("sorter called twice\n");
      abort();
    }
    EM::DynamicPageVector::Vector<TaggedT<T>> v(getOutputSize(),
                                                getBucketSize());
    KWayButterflySort(v.begin(), v.end());
    delete[] batch;
    batch = NULL;
  }

  size_t getOutputSize() { return numTotalBucket * Z; }

  size_t getBucketSize() { return Z; }
};

template <class Iterator>
void KWayButterflyOShuffle(Iterator begin, Iterator end, uint32_t inAuth,
                           uint64_t heapSize) {
  using T = typename std::iterator_traits<Iterator>::value_type;
  size_t N = end - begin;
  if (N <= 512) {
    std::vector<T> Mem(N);
    CopyIn(begin, end, Mem.begin(), inAuth);
    OrShuffle(Mem);
    CopyOut(Mem.begin(), Mem.end(), begin, inAuth + 1);
    return;
  }
  ButterflySorter<Iterator, KWAYBUTTERFLYOSHUFFLE> sorter(begin, end, inAuth,
                                                          heapSize);
  sorter.sort();
}

template <class Iterator>
void KWayButterflySort(Iterator begin, Iterator end, uint32_t inAuth,
                       uint64_t heapSize) {
  using T = typename std::iterator_traits<Iterator>::value_type;
  const uint64_t N = end - begin;
  if (N <= 512) {
    std::vector<T> Mem(N);
    CopyIn(begin, end, Mem.begin(), inAuth);
    BitonicSort(Mem);
    CopyOut(Mem.begin(), Mem.end(), begin, inAuth + 1);
    return;
  }
  ButterflySorter<Iterator, KWAYBUTTERFLYOSORT> sorter(begin, end, inAuth,
                                                       heapSize);
  sorter.sort();
  const auto& mergeRanges = sorter.getMergeSortBatchRanges();
  ExtMergeSort(begin, end, mergeRanges, inAuth + 1,
               heapSize / (sizeof(T) * Vector<T>::item_per_page * 2));
}

template <typename Vec>
void KWayButterflySort(Vec& vec, uint32_t inAuth, uint64_t heapSize) {
  KWayButterflySort(vec.begin(), vec.end(), inAuth, heapSize);
}

template <typename Vec>
void KWayButterflyOShuffle(Vec& vec, uint32_t inAuth, uint64_t heapSize) {
  KWayButterflyOShuffle(vec.begin(), vec.end(), inAuth, heapSize);
}

}  // namespace EM::Algorithm