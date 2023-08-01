#pragma once

#include "external_memory/algorithm/kway_butterfly_sort.hpp"
#include "external_memory/dynamicvector.hpp"
#include "external_memory/noncachedvector.hpp"

/// This file implements the flex-way distribution o-sort algorithm.

namespace EM::Algorithm {
// using NonCachedVector by default
using EM::NonCachedVector::Vector;
/// @brief KWayDistriSort is an oblivious external memory sorting algorithm for
/// data that has been randomly permuted. The implementation incurs
/// approximately 2NlogN exchanges and N/B log_{M/B} N/B page transfers, where N
/// is the size of the input array, B is the page size, and M is the size of the
/// memory.
/// @tparam Iterator only supports NonCachedVector::Iterator
/// @tparam preshuffled whether the input array is preshuffled
/// @param begin the begin iterator of the input array
/// @param end the end iterator of the input array
/// @param inAuth the authentication counter of the input array
/// @param heapSize the size of available memory in bytes
template <class Iterator, bool preshuffled = false>
void KWayDistriSort(Iterator begin, Iterator end, uint32_t inAuth,
                    uint64_t heapSize);

template <class Iterator>
void KWayDistriSortShuffled(Iterator begin, Iterator end, uint32_t inAuth,
                            uint64_t heapSize);

template <typename Vec>
void KWayDistriSort(Vec& vec, uint32_t inAuth = 0,
                    uint64_t heapSize = DEFAULT_HEAP_SIZE);

template <typename Vec>
void KWayDistriSortShuffled(Vec& vec, uint32_t inAuth = 0,
                            uint64_t heapSize = DEFAULT_HEAP_SIZE);

template <class Vec>
struct PagePermReader : public Vec::PrefetchReader {
  Vector<uint64_t> perm;
  Vector<uint64_t>::PrefetchReader permReader;
  uint64_t nextPageIdx = 0;
  uint64_t translatedPageIdx = 0;
  using Iterator = typename Vec::Iterator;
  PagePermReader(Iterator _begin, Iterator _end, uint32_t _counter = 0,
                 size_t _heapSize = DEFAULT_HEAP_SIZE)
      : perm(divRoundUp(_end - _begin, Vec::item_per_page)) {
    Vec::PrefetchReader::init(_begin, _end, _counter);
    if (sizeof(uint64_t) * perm.size() < _heapSize) {
      std::vector<uint64_t> permVec(perm.size());
      for (size_t i = 0; i < perm.size(); ++i) {
        permVec[i] = i;
      }
      // the last page is not permuted
      fisherYatesShuffle(permVec.begin(), permVec.end() - 1);
      CopyOut(permVec.begin(), permVec.end(), perm.begin());
    } else {
      Vector<uint64_t>::Writer permWriter(perm.begin(), perm.end());
      size_t permSize = perm.size() - 1;  // last page not permuted
      Vector<TaggedT<uint64_t>> intermediateVec(permSize);
      size_t batchSize = _heapSize / sizeof(TaggedT<uint64_t>);
      std::vector<std::pair<uint64_t, uint64_t>> mergeRanges;
      {
        std::vector<TaggedT<uint64_t>> batch(batchSize);
        for (size_t i = 0; i < permSize; i += batchSize) {
          size_t batchEnd = std::min(i + batchSize, permSize);
          for (size_t j = i; j < batchEnd; ++j) {
            batch[j - i].v = j;
            batch[j - i].tag = UniformRandom();
          }
          const auto cmpTag = [](const TaggedT<uint64_t>& a,
                                 const TaggedT<uint64_t>& b) {
            return a.tag < b.tag;
          };
          std::sort(batch.begin(), batch.begin() + (batchEnd - i), cmpTag);
          CopyOut(batch.begin(), batch.begin() + (batchEnd - i),
                  intermediateVec.begin() + i);
          mergeRanges.emplace_back(i, batchEnd);
        }
      }
      using MergeReader =
          typename Vector<TaggedT<uint64_t>>::LazyPrefetchReader;
      std::vector<MergeReader> mergeReaders;
      mergeReaders.reserve(mergeRanges.size());
      for (const auto& range : mergeRanges) {
        mergeReaders.emplace_back(intermediateVec.begin() + range.first,
                                  intermediateVec.begin() + range.second);
      }
      std::vector<std::pair<MergeReader*, uint64_t>> heap;
      heap.reserve(mergeRanges.size() + 1);
      for (auto& reader : mergeReaders) {
        reader.init();
        heap.push_back({&reader, reader.get().tag});
      }
      auto cmpmerge = [](const auto& a, const auto& b) {
        return b.second < a.second;
      };
      std::make_heap(heap.begin(), heap.end(), cmpmerge);
      while (!heap.empty()) {
        MergeReader* top = heap[0].first;
        permWriter.write(top->read().getData());
        if (!top->eof()) {
          heap.emplace_back(top, top->get().tag);
          // add a top at the end, which will be swapped to the top by
          // pop_heap
        }
        std::pop_heap(heap.begin(), heap.end(), cmpmerge);
        heap.resize(heap.size() - 1);
      }
      permWriter.write(permSize);
      permWriter.flush();
    }
    permReader.init(perm.begin(), perm.end());
  }

  INLINE uint64_t translatePageIdx(uint64_t pageIdx) {
    if (pageIdx == nextPageIdx) {
      ++nextPageIdx;
    } else {
      Assert(false);  // this is undesired behavior
      permReader.curr = (uint64_t*)UINT64_MAX;
      permReader.it = permReader.it + (pageIdx - nextPageIdx);
      nextPageIdx = pageIdx + 1;
    }
    translatedPageIdx = permReader.read();
    return translatedPageIdx;
  }

  // return the real index of the latest queried item
  uint64_t getLastIdx() {
    return translatedPageIdx * Vec::item_per_page +
           (Vec::PrefetchReader::it - 1).get_page_offset();
  }
};

/// @brief Get the p-quantiles of the input array.
/// @tparam Iterator should support random access
/// @tparam T the type of the input array elements
/// @param begin the begin iterator of the input array
/// @param end the end iterator of the input array
/// @param p the number of quantiles
/// @param inAuth the authentication counter of the input array
/// @return the p-quantiles of the input array
template <typename Iterator,
          typename T = typename std::iterator_traits<Iterator>::value_type>
void getPQuantile(Iterator begin, Iterator end, std::vector<T>& quantiles,
                  uint32_t inAuth) {
  size_t p = quantiles.size() + 1;
  size_t n = end - begin;
  for (size_t k = 1; k < p; ++k) {
    auto it = (begin + (n * k) / p);
    quantiles[k - 1] = it.derefAuth(inAuth);
  }
}

template <typename Iterator,
          typename T = typename std::iterator_traits<Iterator>::value_type>
void getPQuantile(Iterator begin, Iterator end, std::vector<T>& quantiles) {
  size_t p = quantiles.size() + 1;
  size_t n = end - begin;
  for (size_t k = 1; k < p; ++k) {
    auto it = (begin + (n * k) / p);
    quantiles[k - 1] = *it;
  }
}

/// @brief Sample the input array and return the sampled array in ascending
/// order, wrapped in Block<T>.
/// @tparam T the type of the input array elements
/// @tparam reservoir whether to use reservoir sampling
/// @param begin the begin iterator of the input array
/// @param end the end iterator of the input array
/// @param M the size of the memory in number of elements
/// @param alpha the sampling rate
/// @param slackSampling the slack sampling rate for each batch
/// @return Vector<Block<T>> the sampled array in ascending order
template <typename IOIterator, const bool reservoir = false,
          typename T = typename std::iterator_traits<IOIterator>::value_type>
std::vector<Block<T>> sampleForPivots(IOIterator begin, IOIterator end,
                                      size_t M, double alpha,
                                      double slackSampling, size_t p) {
  Assert(alpha > 0);
  Assert(slackSampling * alpha < 1);
  Assert(slackSampling > 1);
  Assert(begin < end);
  using IOVector = typename
      std::remove_reference<decltype(*(IOIterator::getNullVector()))>::type;
  size_t N = end - begin;
  size_t expectedSampleSize = alpha * (double)N;
  size_t sampleSize = 0;
  std::vector<Block<T>> pivots(p - 1);
  using SampleVector =
      Vector<Block<T>, std::max((1UL << 14) - 32, 4 * sizeof(T)), true, true,
             true>;
  M -= 8 + 32 * SampleVector::item_per_page +
       2 * (sampleSize / SampleVector::item_per_page + 2) / sizeof(Block<T>) /
           8;  // subtract the space for the sample array
  Block<T>* Mem = new Block<T>[M];
  if constexpr (reservoir) {
    sampleSize = expectedSampleSize;
  }
  size_t numBatch = divRoundUp(N, M);

  size_t sampleArrSize = (size_t)(slackSampling * M * alpha) * numBatch;
  SampleVector Gamma(sampleArrSize, DUMMY<Block<T>>());
  auto outIt = Gamma.begin();

  size_t Np = N;
  size_t np = expectedSampleSize;
  size_t count = 0;
  typename IOVector::PrefetchReader inputReader(begin, end);
  auto isMarked = [](const Block<T>& element) { return !element.isDummy(); };
  for (size_t i = 0; i < numBatch; ++i) {
    size_t batchSize = std::min(M, N - i * M);
    for (size_t j = 0; j < batchSize; ++j, ++count) {
      if constexpr (reservoir) {
        uint64_t z = UniformRandom(Np - 1);
        bool chosen = z < np;
        CMOV(!chosen, Mem[j], DUMMY<Block<T>>());
        auto realData = Block<T>();
        realData.setData(inputReader.read(), count);
        CMOV(chosen, Mem[j], realData);
        np -= chosen;
        --Np;
      } else {
        uint64_t z = UniformRandom(1, N);
        bool chosen = z <= expectedSampleSize;
        auto realData = Block<T>();
        realData.setData(inputReader.read(), count);
        CMOV(chosen, Mem[j], realData);
        CMOV(!chosen, Mem[j], DUMMY<Block<T>>());
        CMOV(chosen, sampleSize, sampleSize + 1);
      }
    }
    OrCompact(Mem, Mem + batchSize, isMarked);
    size_t len = std::min(batchSize, (size_t)(slackSampling * alpha * M));
    CopyOut(Mem, Mem + len, outIt);
    outIt += len;
  }
  // next, sort the samples
  size_t prefixLen = outIt - Gamma.begin();
  if (prefixLen < M) {
    // can be solved in memory
    CopyIn(Gamma.begin(), outIt, Mem);
    BitonicSort(Mem, Mem + prefixLen);
    getPQuantile(Mem, Mem + sampleSize, pivots);
    delete[] Mem;
    return pivots;
  }
  delete[] Mem;

  KWayButterflySort(Gamma.begin(), outIt, 0, M * sizeof(Block<T>));
  getPQuantile(Gamma.begin(), Gamma.begin() + sampleSize, pivots, 1);
  return pivots;
}

/// @brief A manager class for flex-way distribution o-sort.
/// @tparam T the type of elements to sort
/// @tparam WrappedT the wrapped type of elements to sort
template <typename IOIterator, bool preshuffled>
class KWayDistriSorter {
 private:
  using T = typename std::iterator_traits<IOIterator>::value_type;
  using WrappedT = Block<T>;
  using IOVector = typename
      std::remove_reference<decltype(*(IOIterator::getNullVector()))>::type;
  uint64_t Z;                 // bucket size
  uint64_t numTotalBucket;    // total number of buckets
  uint64_t numPartition;      // number of partitions
  uint64_t partitionSize;     // size of each partition
  uint64_t numRealPerBucket;  // number of real elements per bucket

  uint64_t numBucketFit;   // number of buckets that can fit in the heap
  uint64_t numElementFit;  // number of elements that can fit in the heap
  uint64_t inputCount;  // a counter to track how many elements have been read

  DistriParams distriParams = {};  // parameters for distribution osort

  std::vector<WrappedT> allPivots;  // all pivots for distribution osort
  WrappedT* temp;                   // temporary buffer for mergesplit
  uint8_t* marks;                   // marks for mergesplit
  WrappedT pivots[7];               // pivots for mergesplit
  std::conditional_t<preshuffled, typename IOVector::PrefetchReader,
                     PagePermReader<IOVector>>
      inputReader;                         // input reader
  typename IOVector::Writer outputWriter;  // output writer

 public:
  /// @brief Construct a new KWayDistriSorter object
  /// @param inputBeginIt begin iterator of the input array
  /// @param inputEndIt end iterator of the input array
  /// @param inAuth the counter of the input array for authentication
  /// @param _heapSize the heap size in bytes
  KWayDistriSorter(IOIterator inputBeginIt, IOIterator inputEndIt,
                   uint32_t inAuth = 0, uint64_t _heapSize = DEFAULT_HEAP_SIZE)
      : inputReader(inputBeginIt, inputEndIt, inAuth, _heapSize),
        numElementFit(_heapSize / sizeof(WrappedT)) {
    size_t size = inputEndIt - inputBeginIt;
    distriParams =
        bestDistriParams(size, numElementFit,
                         preshuffled ? 1 : IOVector::item_per_page, sizeof(T));
    Z = distriParams.Z;
    numTotalBucket = distriParams.totalBucket;
    numPartition = distriParams.totalPartition;
    partitionSize = numTotalBucket / numPartition;
    numRealPerBucket = 1 + (size - 1) / numTotalBucket;

    Assert(numTotalBucket % numPartition == 0);
    Assert(numTotalBucket / numPartition * Z <= numElementFit);

    allPivots = sampleForPivots(
        inputBeginIt, inputEndIt,
        numElementFit * sizeof(WrappedT) / (sizeof(WrappedT) + 4),
        distriParams.samplingRatio, distriParams.slackSampling, numPartition);
    temp = new WrappedT[8 * Z];  // temporary array for mergesplit
    marks = new uint8_t[8 * Z];  // temporary array for marks in mergesplit
    // orcompact has 4 byte memory overhead per element
    numElementFit -=
        divRoundUp(Z * 8 * (sizeof(WrappedT) + 1), sizeof(WrappedT)) +
        numPartition;  // subtract the 8 buckets for mergesplit and pivots
    numRealPerBucket = 1 + (size - 1) / numTotalBucket;
    // increment authentication counter when writing back
    outputWriter.init(inputBeginIt, inputEndIt, inAuth + 1);
  }

  ~KWayDistriSorter() {
    delete[] temp;
    delete[] marks;
  }

  /// @brief Base case of flex-way distribution osort when input fits in memory
  /// @tparam Iterator should support random access
  /// @param begin begin iterator of the input array
  /// @param end end iterator of the input array
  /// @param ioLayer current layer of butterfly network by page swap passes
  /// @param innerLayer current layer of butterfly network within the batch
  /// @param pivotStride the stride of pivots for this batch
  /// @param pivotOffset the offset position of the first pivot for this batch
  /// @param stride the stride when fetching buckets in this batch
  template <class Iterator>
  void KWayDistriSortBasic(Iterator begin, Iterator end, size_t ioLayer,
                           size_t innerLayer, size_t pivotStride,
                           size_t pivotOffset, size_t stride = 1) {
    auto& innerWays = distriParams.ways[ioLayer];
    uint64_t way = innerWays[innerLayer];
    if (innerLayer > 0) {
      for (uint64_t i = 0; i < way; ++i) {
        KWayDistriSortBasic(begin + i * stride * Z, end, ioLayer,
                            innerLayer - 1, pivotStride, pivotOffset,
                            stride * way);
      }
    } else {
      if (ioLayer == 0) {
        // tag and pad input
        for (uint64_t i = 0; i < way; ++i) {
          auto it = begin + i * stride * Z;
          for (uint64_t offset = 0; offset < Z; ++offset, ++it) {
            if (offset < numRealPerBucket && !inputReader.eof()) {
              auto& data = inputReader.read();
              uint64_t inputIdx;
              if constexpr (preshuffled) {
                inputIdx = inputCount++;
              } else {
                inputIdx = inputReader.getLastIdx();
              }

              it->setData(data, inputIdx);
              // this count should be the same for the same element in the
              // samples
            } else {
              it->setDummy();
            }
          }
        }
      }
    }
    Iterator KWayIts[8];
    size_t strideCount = 0;
    for (Iterator strideBegin = begin; strideBegin < end;
         strideBegin += stride * way * Z, ++strideCount) {
      for (uint64_t i = 0; i < way; ++i) {
        KWayIts[i] = strideBegin + i * stride * Z;
      }
      for (uint64_t i = 1; i < way; ++i) {
        uint64_t pivotIdx =
            ((strideCount * way + i) * stride + pivotOffset) * pivotStride - 1;
        pivots[i - 1] = allPivots[pivotIdx];
      }
      MergeSplitKWay(KWayIts, way, Z, temp, marks, pivots);
    }
  }

  /// @brief Flex-way distribution osort
  /// @tparam Iterator should support CopyIn and CopyOut
  /// @param begin begin iterator of the input array
  /// @param end end iterator of the input array
  template <class Iterator>
  void KWayDistriSort(Iterator begin, Iterator end) {
    KWayDistriSort(begin, end, distriParams.ways.size());
  }

  template <class Iterator>
  void KWayDistriSort(Iterator begin, Iterator end, size_t ioLayer,
                      size_t stride = 1) {
    bool isLastLayer = ioLayer == distriParams.ways.size();
    size_t numInternalWay = isLastLayer
                                ? numTotalBucket / numPartition
                                : getVecProduct(distriParams.ways[ioLayer]);
    if (ioLayer > 0) {
      for (size_t i = 0; i < numInternalWay; ++i) {
        KWayDistriSort(begin + i * stride * Z, end, ioLayer - 1,
                       stride * numInternalWay);
      }
    }

    size_t strideCount = 0;
    WrappedT* batch = new WrappedT[numInternalWay * Z];
    auto batchBegin = batch;
    auto batchEnd = batch + numInternalWay * Z;
    for (Iterator strideBegin = begin; strideBegin < end;
         strideBegin += stride * numInternalWay * Z, ++strideCount) {
      if (ioLayer) {  // fetch from intermediate ext vector
        for (uint64_t bucketIdx = 0; bucketIdx < numInternalWay; ++bucketIdx) {
          auto extBeginIt = strideBegin + bucketIdx * stride * Z;
          auto intBeginIt = batchBegin + bucketIdx * Z;
          CopyIn(extBeginIt, extBeginIt + Z, intBeginIt, ioLayer - 1);
        }
      }

      if (isLastLayer) {
        // last layer, combine with bitonic sort and output
        static const auto notDummyMark = [](const auto& element) {
          return !element.isDummy();
        };
        static const auto notDummyComp = [](const auto& element,
                                            const auto& unused) {
          return !element.isDummy();
        };
        OrCompact(batchBegin, batchEnd, notDummyMark);
        auto realEnd = std::lower_bound(batchBegin, batchEnd, 0, notDummyComp);
        if constexpr (IO_ROUND == 0) {  // mock
          realEnd = batchBegin + uint64_t((batchEnd - batchBegin) *
                                          (numRealPerBucket - 1) / Z);
        }
        BitonicSort(batchBegin, realEnd);
        for (auto it = batchBegin; it != realEnd; ++it) {
          outputWriter.write(it->getData());
        }

      } else {
        KWayDistriSortBasic(batchBegin, batchEnd, ioLayer,
                            distriParams.ways[ioLayer].size() - 1,
                            stride / partitionSize,
                            strideCount * numInternalWay);
        for (uint64_t bucketIdx = 0; bucketIdx < numInternalWay; ++bucketIdx) {
          auto extBeginIt = strideBegin + bucketIdx * stride * Z;
          auto intBeginIt = batchBegin + bucketIdx * Z;
          CopyOut(intBeginIt, intBeginIt + Z, extBeginIt, ioLayer);
        }
      }
    }
    delete[] batch;

    if (isLastLayer) {
      outputWriter.flush();
    }
  }

  void sort() {
    EM::DynamicPageVector::Vector<Block<T>> v(getOutputSize(), getBucketSize());
    KWayDistriSort(v.begin(), v.end());
  }

  size_t getOutputSize() { return numTotalBucket * Z; }

  size_t getBucketSize() { return Z; }
};

template <class Iterator, bool preshuffled>
void KWayDistriSort(Iterator begin, Iterator end, uint32_t inAuth,
                    uint64_t heapSize) {
  using T = typename std::iterator_traits<Iterator>::value_type;
  size_t N = end - begin;
  if (N <= heapSize / sizeof(T)) {
    std::vector<T> Mem(N);
    CopyIn(begin, end, Mem.begin(), inAuth);
    BitonicSort(Mem);
    CopyOut(Mem.begin(), Mem.end(), begin, inAuth + 1);
    return;
  }
  KWayDistriSorter<Iterator, preshuffled> sorter(begin, end, inAuth, heapSize);
  sorter.sort();
}

template <class Iterator>
void KWayDistriSortShuffled(Iterator begin, Iterator end, uint32_t inAuth,
                            uint64_t heapSize) {
  KWayDistriSort<Iterator, true>(begin, end, inAuth, heapSize);
}

template <typename Vec>
void KWayDistriSort(Vec& vec, uint32_t inAuth, uint64_t heapSize) {
  KWayDistriSort(vec.begin(), vec.end(), inAuth, heapSize);
}

template <typename Vec>
void KWayDistriSortShuffled(Vec& vec, uint32_t inAuth, uint64_t heapSize) {
  KWayDistriSortShuffled(vec.begin(), vec.end(), inAuth, heapSize);
}

}  // namespace EM::Algorithm
