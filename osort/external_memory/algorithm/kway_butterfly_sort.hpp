#pragma once
#include "external_memory/dynamicvector.hpp"
#include "external_memory/noncachedvector.hpp"
#include "external_memory/stdvector.hpp"
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
/// @tparam Reader the reader type of the input array
/// @tparam Writer the writer type of the output array
/// @tparam WrappedT the wrapped type of elements to sort
template <class Reader, class Writer, SortMethod task, const bool singleBatch>
class ButterflySorter {
 private:
  using T = typename Reader::value_type;
  using WrappedT = TaggedT<T>;
  uint64_t Z;                 // bucket size
  uint64_t numTotalBucket;    // total number of buckets
  uint64_t numRealPerBucket;  // number of real elements per bucket
  uint64_t numElementFit;     // number of elements that fit in memory
  KWayButterflyParams KWayParams =
      {};  // parameters for flex-way butterfly o-sort

  static constexpr bool needsExtMerge =
      !singleBatch &&
      task == KWAYBUTTERFLYOSORT;  // whether to use external merge sort
  uint64_t maxBatchBucket = 0;     // maximum number of buckets in a batch

  std::vector<std::pair<typename Writer::iterator_type,
                        typename Writer::iterator_type>>
      mergeSortRanges;  // pairs of iterators that specifies each sorted range
                        // in the first layer of external-memory merge sort

  Reader& inputReader;           // input reader
  Writer& outputWriter;          // output writer
  size_t batchSize;              // batch size in number of elements
  WrappedT* batch;               // batch for sorting
  WrappedT* mergeSplitTemp;      // temporary space for merge-splitting
  uint8_t* mergeSplitMarksTemp;  // temporary space for merge-splitting marks

 public:
  /// @brief Construct a new Butterfly Sorter object
  /// @param reader the reader of the input array
  /// @param writer the writer of the output array
  /// @param _KWayParams the parameters for flex-way butterfly o-sort
  /// @param _heapSize the heap size in bytes
  ButterflySorter(Reader& reader, Writer& writer,
                  const KWayButterflyParams& _KWayParams, uint64_t _heapSize)
      : inputReader(reader),
        outputWriter(writer),
        KWayParams(_KWayParams),
        Z(_KWayParams.Z),
        numElementFit((_heapSize - 16 * Z) / sizeof(WrappedT) - 8 * Z) {
    size_t size = reader.size();
    Z = KWayParams.Z;
    for (auto& ways : KWayParams.ways) {
      maxBatchBucket = std::max(maxBatchBucket, getVecProduct(ways));
    }
    batchSize = singleBatch ? maxBatchBucket * Z : numElementFit;
    batch = new WrappedT[batchSize];
    mergeSplitTemp = new WrappedT[8 * Z];
    mergeSplitMarksTemp = new uint8_t[8 * Z];

    // reserve space for 8 buckets of elements and marks
    numTotalBucket = KWayParams.totalBucket;
    numRealPerBucket = 1 + (size - 1) / numTotalBucket;
  }

  ~ButterflySorter() {
    if (batch) {
      delete[] batch;
    }
    if (mergeSplitTemp) {
      delete[] mergeSplitTemp;
    }
    if (mergeSplitMarksTemp) {
      delete[] mergeSplitMarksTemp;
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

    for (uint64_t j = 0; j < wayBucket; ++j) {
      for (uint64_t i = 0; i < way; ++i) {
        KWayIts[i] = begin + (i * wayBucket + j) * Z;
      }
      MergeSplitKWay(KWayIts, way, Z, mergeSplitTemp, mergeSplitMarksTemp);
    }
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
    size_t localBatchSize = numInternalWay * Z;
    size_t batchPerEnclave = 1;
    if (task == KWAYBUTTERFLYOSORT && isLastLayer) {
      batchPerEnclave = numElementFit / localBatchSize;
      // maximize the chunksize at the last layer
    }
    size_t batchCount = divRoundUp(size, localBatchSize);

    for (uint64_t batchIdx = 0; batchIdx < batchCount; ++batchIdx) {
      uint64_t intBatchIdx = batchIdx % batchPerEnclave;
      auto batchBegin = batch + intBatchIdx * localBatchSize;
      if (ioLayer) {  // fetch from intermediate ext vector
        for (uint64_t bucketIdx = 0; bucketIdx < numInternalWay; ++bucketIdx) {
          auto extBeginIt = begin + (batchIdx + fetchInterval * bucketIdx) * Z;
          auto intBeginIt = batchBegin + bucketIdx * Z;
          CopyIn(extBeginIt, extBeginIt + Z, intBeginIt, ioLayer - 1);
        }
      }
      KWayButterflySortBasic(batchBegin, batchBegin + localBatchSize, ioLayer,
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
        if constexpr (task == KWAYBUTTERFLYOSORT) {
          if (intBatchIdx == batchPerEnclave - 1 ||
              batchIdx == batchCount - 1) {
            // sort the batch and write to first layer of merge sort
            const auto cmpVal = [](const auto& a, const auto& b) {
              return a.v < b.v;
            };
            size_t BNCountPerBatch = batchIdx == batchCount - 1
                                         ? batchIdx % batchPerEnclave + 1
                                         : batchPerEnclave;
            Assert(localBatchSize * BNCountPerBatch <= batchSize);
            auto realEnd =
                partitionDummy(batch, batch + localBatchSize * BNCountPerBatch);
            // partition dummies to the end
            Assert(realEnd <= batch + batchSize);
            std::sort(batch, realEnd, cmpVal);
            if constexpr (needsExtMerge) {
              auto mergeSortReaderBeginIt = outputWriter.it;
              for (auto it = batch; it != realEnd; ++it) {
                outputWriter.write(it->getData());
              }

              mergeSortRanges.emplace_back(mergeSortReaderBeginIt,
                                           outputWriter.it);
            } else {
              for (auto it = batch; it != realEnd; ++it) {
                outputWriter.write(it->getData());
              }
            }
          }
        }
      } else {  // not last layer, write to intermediate ext vector
        Assert(!singleBatch);
        Assert(batchPerEnclave == 1);
        for (uint64_t bucketIdx = 0; bucketIdx < numInternalWay; ++bucketIdx) {
          auto extBeginIt = begin + (batchIdx + fetchInterval * bucketIdx) * Z;
          auto intBeginIt = batch + bucketIdx * Z;
          CopyOut(intBeginIt, intBeginIt + Z, extBeginIt, ioLayer);
        }
      }
    }

    if (isLastLayer) {
      outputWriter.flush();
    }
  }

  void sort() {
    if (batch == NULL) {
      abort();
    }
    if (singleBatch) {
      typename EM::DynamicPageVector::Vector<TaggedT<T>>::Iterator dummyBegin =
          0;
      typename EM::DynamicPageVector::Vector<TaggedT<T>>::Iterator dummyEnd =
          getOutputSize();
      KWayButterflySort(dummyBegin, dummyEnd);
    } else {
      EM::DynamicPageVector::Vector<TaggedT<T>> v(getOutputSize(),
                                                  getBucketSize());
      KWayButterflySort(v.begin(), v.end());
    }

    delete[] batch;
    batch = NULL;
    delete[] mergeSplitTemp;
    mergeSplitTemp = NULL;
    delete[] mergeSplitMarksTemp;
    mergeSplitMarksTemp = NULL;
  }

  size_t getOutputSize() { return numTotalBucket * Z; }

  size_t getBucketSize() { return Z; }
};

template <class Reader, class Writer>
void KWayButterflyOShuffle(Reader& reader, Writer& writer,
                           uint64_t heapSize = DEFAULT_HEAP_SIZE) {
  using T = typename Reader::value_type;
  using WrappedT = TaggedT<T>;
  size_t N = reader.size();
  if (N <= 512) {
    StdVector<T> Mem(N);
    for (auto& m : Mem) {
      m = reader.read();
    }
    OrShuffle(Mem);
    for (auto& m : Mem) {
      writer.write(m);
    }
    writer.flush();
    return;
  }
  heapSize = std::min(heapSize, std::max(N, 4096UL) * sizeof(WrappedT) * 2);
  uint64_t numElementFit = heapSize / sizeof(WrappedT);
  WrappedT* batch = new WrappedT[numElementFit];  // avoid fragmentation
  const KWayButterflyParams& KWayParams =
      bestKWayButterflyParams(N, numElementFit, sizeof(T));
  delete[] batch;
  bool singleBatch = KWayParams.ways.size() == 1;
  if (singleBatch) {
    ButterflySorter<Reader, Writer, KWAYBUTTERFLYOSHUFFLE, true> sorter(
        reader, writer, KWayParams, heapSize);
    sorter.sort();
  } else {
    ButterflySorter<Reader, Writer, KWAYBUTTERFLYOSHUFFLE, false> sorter(
        reader, writer, KWayParams, heapSize);
    sorter.sort();
  }
}

template <class Iterator>
void KWayButterflyOShuffle(Iterator begin, Iterator end, uint32_t inAuth,
                           uint64_t heapSize) {
  using IOVector = typename Iterator::vector_type;
  using Reader = typename IOVector::PrefetchReader;
  using Writer = typename IOVector::Writer;
  if constexpr (Iterator::random_access) {
    if (end - begin <= 512) {
      OrShuffle(begin, end);
      return;
    }
  }
  Reader reader(begin, end, inAuth);
  Writer writer(begin, end, inAuth + 1);
  KWayButterflyOShuffle(reader, writer, heapSize);
}

template <class Reader, class Writer>
void KWayButterflySort(Reader& reader, Writer& writer,
                       uint64_t heapSize = DEFAULT_HEAP_SIZE) {
  using T = typename Reader::value_type;
  using WrappedT = TaggedT<T>;
  const uint64_t N = reader.size();
  if (N <= 512) {
    StdVector<T> Mem(N);
    for (auto& m : Mem) {
      m = reader.read();
    }
    BitonicSort(Mem);
    for (auto& m : Mem) {
      writer.write(m);
    }
    writer.flush();
    return;
  }
  heapSize = std::min(heapSize, std::max(N, 4096UL) * sizeof(WrappedT) * 2);
  uint64_t numElementFit = heapSize / sizeof(WrappedT);
  WrappedT* batch = new WrappedT[numElementFit];  // avoid fragmentation
  const KWayButterflyParams& KWayParams =
      bestKWayButterflyParams(N, numElementFit, sizeof(T));
  delete[] batch;
  bool singleBatch = KWayParams.ways.size() == 1;

  if (singleBatch) {
    ButterflySorter<Reader, Writer, KWAYBUTTERFLYOSORT, true> sorter(
        reader, writer, KWayParams, heapSize);
    sorter.sort();
  } else {
    using ShuffleVector = typename EM::NonCachedVector::Vector<T>;
    ShuffleVector mergeSortFirstLayer(N);
    using ShuffleWriter = typename ShuffleVector::Writer;
    ShuffleWriter mergeSortFirstLayerWriter(mergeSortFirstLayer.begin(),
                                            mergeSortFirstLayer.end());
    ButterflySorter<Reader, ShuffleWriter, KWAYBUTTERFLYOSORT, false> sorter(
        reader, mergeSortFirstLayerWriter, KWayParams, heapSize);
    sorter.sort();
    const auto& mergeRanges = sorter.getMergeSortBatchRanges();
    ExtMergeSort(
        writer, mergeRanges,
        heapSize /
            (sizeof(T) * EM::NonCachedVector::Vector<T>::item_per_page * 2));
  }
}

template <class Iterator>
void KWayButterflySort(Iterator begin, Iterator end, uint32_t inAuth,
                       uint64_t heapSize) {
  using IOVector = typename Iterator::vector_type;
  using Reader = typename IOVector::PrefetchReader;
  using Writer = typename IOVector::Writer;
  if constexpr (Iterator::random_access) {
    if (end - begin <= 512) {
      BitonicSort(begin, end);
      return;
    }
  }
  Reader reader(begin, end, inAuth);
  Writer writer(begin, end, inAuth + 1);
  KWayButterflySort(reader, writer, heapSize);
}

template <typename Vec>
void KWayButterflySort(Vec& vec, uint32_t inAuth, uint64_t heapSize) {
  KWayButterflySort(vec.begin(), vec.end(), inAuth, heapSize);
}

template <typename Vec>
void KWayButterflyOShuffle(Vec& vec, uint32_t inAuth, uint64_t heapSize) {
  KWayButterflyOShuffle(vec.begin(), vec.end(), inAuth, heapSize);
}

template <typename Iterator>
void KWayButterflyOShuffleInternal(Iterator begin, Iterator end) {
  KWayButterflyOShuffle(begin, end, 0, 512UL << 30);
}

template <typename Vec>
void KWayButterflyOShuffleInternal(Vec& vec) {
  KWayButterflyOShuffleInternal(vec.begin(), vec.end());
}

template <typename Iterator>
void KWayButterflySortInternal(Iterator begin, Iterator end) {
  KWayButterflySort(begin, end, 0, 512UL << 30);
}

template <typename Vec>
void KWayButterflySortInternal(Vec& vec) {
  KWayButterflySortInternal(vec.begin(), vec.end());
}

}  // namespace EM::Algorithm