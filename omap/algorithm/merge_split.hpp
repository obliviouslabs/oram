#pragma once
#include <vector>

#include "edge_rec.hpp"
#include "element.hpp"
#include "external_memory/stdvector.hpp"
#include "static_sort.hpp"

/**
 * @brief This file implements two-way and multi-way mergesplit algorithms,
 * described in https://eprint.iacr.org/2023/1258.pdf.
 * Specifically, the algorithm may perform 8-way mergesplit with (0.5 +
 * o(1))NlogN oblivious swaps.
 */

namespace Algorithm {

/// @brief A special case of interleave. Two way means we need not generate the
/// Euler graph, but can swap elements in one pass. It also supports non-power
/// of two input size.
/// @tparam Iterator should support random access
/// @tparam Check type of isMarked
/// @param begin begin iterator
/// @param end end iterator
/// @param isMarked a lambda expression that specifies whether an element has
/// key 0.
template <class Iterator, class Check>
void InterleaveTwoWay(Iterator begin, Iterator end, const Check& isMarked) {
  size_t size = end - begin;
  Assert(size % 2 == 0);
  if (size == 2) {
    bool swapFlag = !isMarked(*begin);
    obliSwap(swapFlag, *begin, *(begin + 1));
    return;
  }
  size_t leftHalfSize = (size + 2) / 4 * 2;
  bool targetBit;
  Iterator mid = begin + leftHalfSize;
  Iterator left = begin, right = mid;
  targetBit = !isMarked(*left);
  switch (size % 4) {
    case 0:
      ++left;
      ++right;
      break;
    case 2:
      left += 2;
      break;
  }
  // every time left tag bit != right tag bit, we swap so that the left tag bit
  // becomes target bit, then we negate the target bit. Eventually the number of
  // 0 and 1 will be the same on both half.
  for (; right != end; ++left, ++right) {
    bool leftTag = isMarked(*left);
    bool rightTag = isMarked(*right);
    bool swapFlag = leftTag != targetBit;
    targetBit = rightTag != swapFlag;
    obliSwap(swapFlag, *left, *right);
  }
  InterleaveTwoWay(begin, mid, isMarked);
  InterleaveTwoWay(mid, end, isMarked);
}

/// @brief A special case of interleave. Similar to the overloaded function
/// above, but allows the input to come from two separate "buckets," which
/// minimizes data movement.
/// @tparam Iterator should support random access
/// @tparam Check type of isMarked
/// @param beginLeft begin iterator of the left bucket
/// @param beginRight begin iterator of the right bucket
/// @param Z bucket size
template <class Iterator, class Check>
void InterleaveTwoWay(Iterator beginLeft, Iterator beginRight, size_t Z,
                      const Check& isMarked) {
  Assert(Z % 2 == 0);
  bool targetBit;
  Iterator endLeft = beginLeft + Z;
  Iterator endRight = beginRight + Z;
  Iterator left = beginLeft;
  Iterator right = beginRight;

  targetBit = isMarked(*right);
  ++left;
  ++right;

  // every time left tag bit != right tag bit, we swap so that the left tag bit
  // becomes target bit, then we negate the target bit. Eventually the number of
  // 0 and 1 will be the same on both half.
  for (; right != endRight; ++left, ++right) {
    bool leftTag = isMarked(*left);
    bool rightTag = isMarked(*right);
    bool swapFlag = leftTag != targetBit;
    targetBit = rightTag != swapFlag;
    obliSwap(swapFlag, *left, *right);
  }

  InterleaveTwoWay(beginLeft, endLeft, isMarked);
  InterleaveTwoWay(beginRight, endRight, isMarked);
}

/// @brief Permute sorts elements in [begin, end) according to their marks.
/// Although the algorithm is O(nlogn) in theory, it is slower than a general
/// sorting network when n is small (e.g., <= 8).
/// @tparam Iterator should support random access
/// @tparam MarkIterator should support random access
/// @param begin begin iterator
/// @param end end iterator
/// @param marksBegin begin iterator of the marks
/// @param marksEnd end iterator of the marks
template <class Iterator, class MarkIterator>
void Permute(Iterator begin, Iterator end, MarkIterator marksBegin,
             MarkIterator marksEnd, int level = 0) {
  uint64_t size = marksEnd - marksBegin;
  Assert(size <= 16);
  Assert(marksEnd - marksBegin == end - begin);
  if (size == 1) {
    return;
  } else if (size == 2) {
    bool swapFlag = *marksBegin & 0xFU;
    obliSwap(swapFlag, *begin, *(begin + 1));
    obliSwap(swapFlag, *marksBegin, *(marksBegin + 1));
    return;
  }
  size_t halfSize = (size + 1) / 2;
  EdgeRec<uint64_t> rec(halfSize);
  Iterator mid = begin + halfSize;
  const MarkIterator marksMid = marksBegin + halfSize;
  MarkIterator marksLeft = marksBegin, marksRight = marksMid;
  for (; marksLeft != marksMid - 1; ++marksLeft, ++marksRight) {
    uint8_t leftLargeMask = -(uint8_t)((*marksLeft & 0xF) >= halfSize);
    *marksLeft -= leftLargeMask & halfSize;

    uint8_t rightLargeMask = -(uint8_t)((*marksRight & 0xF) >= halfSize);
    *marksRight -= rightLargeMask & halfSize;
    rec.flipEdge(*marksLeft & 0xF, *marksRight & 0XF);
    *marksLeft |= (0X10U << level) & leftLargeMask;
    *marksRight |= (0X10U << level) & rightLargeMask;
  }

  uint8_t leftLargeMask = -(uint8_t)((*marksLeft & 0xF) >= halfSize);
  *marksLeft -= leftLargeMask & halfSize;
  *marksLeft |= (0X10U << level) & leftLargeMask;

  uint8_t lastLeft = *marksLeft & 0xFU;
  uint8_t lastRight;
  if (size & 1) {
    lastRight = size - halfSize;
  } else {
    uint8_t rightLargeMask = -(uint8_t)((*marksRight & 0xF) >= halfSize);
    *marksRight -= rightLargeMask & halfSize;
    lastRight = *marksRight & 0xF;
    *marksRight |= (0X10U << level) & rightLargeMask;
  }
  rec.flipEdge(lastLeft, lastRight);
  EdgeRec<uint64_t> path = rec.EulerPath(halfSize);
  // if the last edge is in the wrong direction, reverse the path orientation
  // rather than do the swap
  uint64_t reversePathMask = path.retrieveAndFlipEdge(lastLeft, lastRight);
  path.flipEdge(-reversePathMask);
  Iterator left = begin, right = mid;
  for (marksLeft = marksBegin, marksRight = marksMid; marksLeft != marksMid - 1;
       ++marksLeft, ++marksRight, ++left, ++right) {
    uint8_t leftMark = *marksLeft & 0xFU;
    uint8_t rightMark = *marksRight & 0xFU;
    bool swapFlag = path.retrieveAndFlipEdge(leftMark, rightMark);
    obliSwap(swapFlag, *marksLeft, *marksRight);
    obliSwap(swapFlag, *left, *right);
  }
  Permute(begin, mid, marksBegin, marksMid, level + 1);
  Permute(mid, end, marksMid, marksEnd, level + 1);
  left = begin, right = mid;
  for (marksLeft = marksBegin, marksRight = marksMid; marksRight != marksEnd;
       ++marksLeft, ++marksRight, ++left, ++right) {
    bool swapFlag = !(*marksRight & (0X10U << level));

    obliSwap(swapFlag, *left, *right);
    obliSwap(swapFlag, *marksLeft, *marksRight);
  }
}

/// @brief InterleaveBaseCase sorts elements according to their marks
/// @tparam Iterator should support random access
/// @tparam MarkIterator should support random access
/// @param begin begin iterator
/// @param end end iterator
/// @param marksBegin begin iterator of the marks
/// @param marksEnd end iterator of the marks
template <class Iterator, class MarkIterator>
INLINE void InterleaveBaseCase(Iterator begin, Iterator end,
                               MarkIterator marksBegin, MarkIterator marksEnd,
                               const uint16_t k) {
  static constexpr StaticSort<3> staticSort3;
  static constexpr StaticSort<4> staticSort4;
  static constexpr StaticSort<5> staticSort5;
  static constexpr StaticSort<6> staticSort6;
  static constexpr StaticSort<7> staticSort7;
  static constexpr StaticSort<8> staticSort8;
  auto iter_pair = std::make_pair(begin, marksBegin);
  // k doesn't need to be oblivious
  switch (k) {
    case 3:
      staticSort3(iter_pair);
      break;
    case 4:
      staticSort4(iter_pair);
      break;
    case 5:
      staticSort5(iter_pair);
      break;
    case 6:
      staticSort6(iter_pair);
      break;
    case 7:
      staticSort7(iter_pair);
      break;
    case 8:
      staticSort8(iter_pair);
      break;
    default:
      throw std::runtime_error("Unsupported way of MergeSplit");
  }
}

/// @brief Interleave rearranges elements in [begin, end) so that their
/// corresponding marks appear like 0 1 2 ... k-1 0 1 2 ... k-1 0 ...
/// @tparam Iterator should support random access
/// @tparam MarkIterator should support random access
/// @param begin begin iterator
/// @param end end iterator
/// @param marksBegin begin iterator of the marks
/// @param marksEnd end iterator of the marks
/// @param k number of ways
template <class Iterator, class MarkIterator>
void Interleave(Iterator begin, Iterator end, MarkIterator marksBegin,
                MarkIterator marksEnd, const uint16_t k) {
  uint64_t size = marksEnd - marksBegin;
  Assert(size % k == 0);
  Assert(size / k == GetNextPowerOfTwo(size / k));
  Assert(marksEnd - marksBegin == end - begin);
  if (size == k) {
    // 23% of runtime
    InterleaveBaseCase(begin, end, marksBegin, marksEnd, k);
    return;
  }
  EdgeRec<uint64_t> rec(k);
  Iterator mid = begin + size / 2;
  MarkIterator marksMid = marksBegin + size / 2;

  // 3% of runtime
  for (MarkIterator marksLeft = marksBegin, marksRight = marksMid;
       marksRight != marksEnd; ++marksLeft, ++marksRight) {
    rec.flipEdge(*marksLeft, *marksRight);
  }

  // 13% of runtime
  EdgeRec<uint64_t> path = rec.EulerPath(size / 2);

  // if the first edge is in the wrong direction, reverse the path orientation
  // rather than do the swap
  uint64_t reversePathMask = path.retrieveAndFlipEdge(*marksBegin, *marksMid);
  path.flipEdge(-reversePathMask);

  Iterator left = begin + 1, right = mid + 1;

  for (MarkIterator marksLeft = marksBegin + 1, marksRight = marksMid + 1;
       marksRight != marksEnd; ++marksLeft, ++marksRight, ++left, ++right) {
    bool swapFlag = path.retrieveAndFlipEdge(*marksLeft, *marksRight);

    // 5%
    obliSwap(swapFlag, *marksLeft, *marksRight);

    // 56%
    obliSwap(swapFlag, *left, *right);
  }
  Interleave(begin, mid, marksBegin, marksMid, k);
  Interleave(mid, end, marksMid, marksEnd, k);
}

/// @brief Partition elements into two parts according to their marks using
/// Interleave
/// @tparam Iterator should support random access
/// @tparam Check type of the check function
/// @param begin begin iterator
/// @param end end iterator
/// @param isMarked lambda function that returns true if the element is marked
template <class Iterator, class Check>
void InterleaveTwoWayPartition(Iterator begin, Iterator end,
                               const Check& isMarked) {
  uint64_t size = (end - begin);
  Assert(size % 2 == 0);
  uint64_t Z = size / 2;
  InterleaveTwoWay(begin, end, isMarked);
  for (auto it1 = begin + 1, it2 = begin + Z + Z % 2; it2 < end;
       it1 += 2, it2 += 2) {
    swap(*it1, *it2);
  }
}

/// @brief Partition elements into two parts according to their marks, elements
/// are originally split into two parts
/// @tparam Iterator should support random access
/// @tparam Check type of the check function
/// @param beginLeft begin iterator of the left part
/// @param beginRight begin iterator of the right part
/// @param Z size of either part
/// @param isMarked lambda function that returns true if the element is marked
template <class Iterator, class Check>
void InterleaveTwoWayPartition(Iterator beginLeft, Iterator beginRight,
                               size_t Z, const Check& isMarked) {
  InterleaveTwoWay(beginLeft, beginRight, Z, isMarked);
  auto endRight = beginRight + Z;
  for (auto it1 = beginLeft + 1, it2 = beginRight + Z % 2; it2 < endRight;
       it1 += 2, it2 += 2) {
    std::swap(*it1, *it2);
  }
}

/**
 * Partition array such that all elements with tag & bitMask == 0 are at the
 * left half, and all elements with tag & bitMask == 1 are at the right half
 * (i.e., 0000...1111...). Note that the left and right half of the subarrays
 * doesn't need to be consecutive. This allows us to merge split in place.
 * @param[in] bitMask the mask of a single bit to partition
 * @param[in] beginLeft the beginning iterator (inclusive) of the left half to
 * partition
 * @param[in] beginRight the beginning iterator (inclusive) of the right half to
 * partition
 * @tparam Z the bucketsize
 * @pre size of the (sub)array must be a power of two
 *
 * */

/// @brief Partition array such that all elements with tag & bitMask == 0 are at
/// the left half, and all elements with tag & bitMask == 1 are at the right
/// half (i.e., 0000...1111...).
/// @tparam Iterator should support random access
/// @tparam Indicator is either a uint64_t bit mask or the type of pivot
/// @param begin begin iterator
/// @param end end iterator
/// @param indicator bit mask or pivot
/// @param method partition method
template <class Iterator, typename Indicator>
void MergeSplitInPlace(Iterator begin, Iterator end, Indicator indicator) {
  // Balance the number of 0 and 1 at the indicator bit
  uint64_t markCount = 0;
  uint64_t n = end - begin;
  uint64_t Z = n / 2;
  for (auto it = begin; it != end; ++it) {
    markCount += it->setAndGetMarked(indicator);
    // for struct with separate flag, this
    // should also update the flag
  }
  const bool dir = markCount > Z;
  uint64_t diff = Z - markCount;
  obliMove(dir, diff, -diff);

  for (auto it = begin; it != end; ++it) {
    bool isMarked = it->isMarked(indicator);
    bool isDummy = it->IsDummy();
    bool changeMark = isDummy & (!!diff) & (isMarked == dir);
    it->condChangeMark(changeMark, indicator);
    diff -= (uint64_t)changeMark;
  }

  if (diff) {
    throw std::runtime_error(
        "MergeSplit fails. Consider padding more filler elements.");
  }

  const auto isMarked = [indicator = indicator](const auto& element) {
    return element.isMarked(indicator);
  };
  // Partition 0 and 1 to even and odd bits

  InterleaveTwoWayPartition(begin, end, isMarked);
}

/// @brief A special case of multi-way mergesplit. Note that the left and right
/// half of the subarrays can be separate. This allows us to mergesplit in
/// place.
/// @tparam Iterator should support random access
/// @tparam Indicator is either a uint64_t bit mask or the type of pivot
/// @param beginLeft begin iterator of the left part
/// @param beginRight begin iterator of the right part
/// @param Z size of either part
/// @param indicator bit mask or pivot
/// @param method partition method
template <class Iterator, typename Indicator>
void MergeSplitTwoWay(Iterator beginLeft, Iterator beginRight, size_t Z,
                      Indicator indicator) {
  // Balance the number of 0 and 1 at the indicator bit
  uint64_t markCount = 0;
  Iterator endLeft = beginLeft + Z, endRight = beginRight + Z;
  Iterator end = endLeft;
  for (auto it = beginLeft;; ++it) {
    if (it == end) {
      if (it == endRight) {
        break;
      }
      it = beginRight;
      end = endRight;
    }
    markCount +=
        it->setAndGetMarked(indicator);  // for struct with separate flag, this
                                         // should also update the flag
  }
  const bool dir = markCount > Z;
  uint64_t diff = Z - markCount;
  obliMove(dir, diff, -diff);

  end = endLeft;
  for (auto it = beginLeft;; ++it) {
    if (it == end) {
      if (it == endRight) {
        break;
      }
      it = beginRight;
      end = endRight;
    }
    bool isMarked = it->isMarked(indicator);
    bool isDummy = it->IsDummy();
    bool changeMark = isDummy & (!!diff) & (isMarked == dir);
    it->condChangeMark(changeMark, indicator);
    diff -= (uint64_t)changeMark;
  }

  if (diff) {
    throw std::runtime_error(
        "MergeSplit fails. Consider padding more filler elements.");
  }

  const auto isMarked = [indicator = indicator](const auto& element) {
    return element.isMarked(indicator);
  };
  // Partition 0 and 1 to even and odd bits

  InterleaveTwoWayPartition(beginLeft, beginRight, Z, isMarked);
}

/// @brief Multi-way mergesplit algorithm.
/// @tparam Iterator should support random access
/// @tparam PivotIterator should support random access
/// @param begins specify the beginning of all the k buckets.
/// @param k is the number of ways.
/// @param Z is the bucket size.
/// @param temp temporary array for mergesplit
/// @param marks temporary array for marks
/// @param pivotBegin is the beginning iterator of the vector of k-1 pivots
template <class Iterator, class PivotIterator = void*,
          typename T = typename std::iterator_traits<Iterator>::value_type>
void MergeSplitKWay(const Iterator* begins, const size_t k, const size_t Z,
                    T* temp, uint8_t* marks, PivotIterator pivotBegin = NULL) {
  Assert(k >= 2);
  Assert(k <= 8);
  if (k == 2) {
    // special case
    // avoids Euler tour search and also minimizes data movement
    if constexpr (std::is_same<PivotIterator, void*>::value) {
      MergeSplitTwoWay(begins[0], begins[1], Z, 1UL);
      // update the tags
      for (auto it = begins[0]; it != begins[0] + Z; ++it) {
        it->tag = (int64_t)it->tag >> 1;
        // arithmetic shift to keep the dummy flag
      }
      for (auto it = begins[1]; it != begins[1] + Z; ++it) {
        it->tag = (int64_t)it->tag >> 1;
      }
    } else {
      MergeSplitTwoWay(begins[0], begins[1], Z, *pivotBegin);
    }

    return;
  }
  Assert(temp);
  Assert(marks != nullptr);
  m256i remainCounts = mm256_set1_epi32(Z);  // remaining counts of each key
  auto* tempIt = temp;
  auto* marksIt = marks;
  for (size_t i = 0; i < k; ++i) {
    auto begin = begins[i];
    for (auto it = begin; it != begin + Z; ++it, ++marksIt, ++tempIt) {
      std::memcpy(tempIt, &(*it), sizeof(T));
      uint8_t isDummy = tempIt->IsDummy();
      uint8_t mark;
      if constexpr (std::is_same<PivotIterator, void*>::value) {
        mark = tempIt->getMarkAndUpdate(k);
      } else {
        mark = 0;
        for (auto pivotIt = pivotBegin; pivotIt != pivotBegin + (k - 1);
             ++pivotIt) {
          mark += *pivotIt < *tempIt;
        }
      }
      mark |= -isDummy;  // set to -1 if it's dummy

      Assert(isDummy <= 1);
      *marksIt = mark;
      remainCounts = mm256_decrement_epi32_var_indx(remainCounts, mark);
    }
  }

  if (mm256_contain_le_zero(remainCounts)) {
    throw std::runtime_error(
        "MergeSplit fails. Consider padding more filler elements.");
  }

  int currMark = 0;
  uint32_t currRemain = mm256_extract_epi32_var_indx(remainCounts, 0);
  // assume that there's at least one dummy for each mark
  for (marksIt = marks; marksIt != marks + Z * k; ++marksIt) {
    uint8_t isDummy = (*marksIt == (uint8_t)-1);
    obliMove(isDummy, *marksIt, (uint8_t)currMark);
    currRemain -= isDummy;
    currMark += !currRemain;
    uint32_t remainCount = mm256_extract_epi32_var_indx(remainCounts, currMark);
    obliMove(!currRemain, currRemain, remainCount);
    Assert(currRemain > 0);
  }
  Interleave(temp, temp + k * Z, marks, marks + k * Z, k);
  for (size_t i = 0; i < Z; ++i) {
    for (size_t j = 0; j < k; ++j) {
      memcpy(&*(begins[j] + i), temp + (i * k + j), sizeof(T));
    }
  }
}

/// @brief Multi-way mergesplit algorithm.
/// @tparam Iterator should support random access
/// @tparam PivotIterator should support random access
/// @param begins specify the beginning of all the k buckets.
/// @param k is the number of ways.
/// @param Z is the bucket size.
/// @param pivotBegin is the beginning iterator of the vector of k-1 pivots
template <class Iterator, class PivotIterator = void*>
void MergeSplitKWay(const Iterator* begins, const size_t k, const size_t Z,
                    PivotIterator pivotBegin = NULL) {
  Assert(k >= 2);
  Assert(k <= 8);
  using T = typename std::iterator_traits<Iterator>::value_type;
  if (k != 2) {
    T* temp = new T[k * Z];               // temporary array for mergesplit
    uint8_t* marks = new uint8_t[k * Z];  // temporary array for marks
    MergeSplitKWay(begins, k, Z, temp, marks, pivotBegin);
    delete[] marks;
    delete[] temp;
  } else {
    MergeSplitKWay(begins, k, Z, (T*)NULL, NULL, pivotBegin);
  }
}
}  // namespace Algorithm