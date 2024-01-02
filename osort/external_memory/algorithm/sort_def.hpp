#pragma once
#include "common/defs.hpp"
#include "common/dummy.hpp"
#include "common/encrypted.hpp"
#include "common/mov_intrinsics.hpp"
#include "defs_for_test.hpp"

/// This file defines some macros and data structures for sorting and shuffling
/// algorithms.

// size of each element in bytes
#ifndef ELEMENT_SIZE
#define ELEMENT_SIZE 128
#endif

namespace EM::Algorithm {
enum SortMethod {
  CABUCKETSORT,                   // cache agnostic bucket sort
  BITONICSORT,                    // bitonic sort
  ORSHUFFLE,                      // ORShuffle
  CABUCKETSHUFFLE,                // cache agnostic bucket shuffle
  BITONICSHUFFLE,                 // bitonic shuffle
  KWAYDISTRIBUTIONOSORT,          // flex-way distribution o-sort
  KWAYDISTRIBUTIONOSORTSHUFFLED,  // flex-way distribution o-sort when input
                                  // is preshuffled
  KWAYBUTTERFLYOSORT,             // flex-way butterfly o-sort
  KWAYBUTTERFLYOSHUFFLE,          // flex-way butterfly o-shuffle
  UNOPTBITONICSORT,               // unoptimized bitonic sort
  EXTMERGESORT,                   // external memory merge sort
  OTHER                           // other algorithms
};

enum PartitionMethod {
  INTERLEAVE_PARTITION,  // partition by interleaving
  OR_COMPACT,            // partition by OR compact
  GOODRICH_COMPACT,      // partition using goodrich's method (not the external
                         // memory version)
  BITONIC                // partition using bitonic sort
};

/// @brief conditionally swap two values
/// @tparam perf if true, increment swap counter
/// @param cond condition for swapping
/// @param v1 first value
/// @param v2 second value
template <bool perf = true>
INLINE void condSwap(const auto& cond, auto& v1, auto& v2) {
  if constexpr (perf) {
    PERFCTR_INCREMENT(swapCount);
  }

  obliSwap(cond, v1, v2);
}

/// @brief swap two values
template <bool perf = true>
INLINE void swap(auto& v1, auto& v2) {
  if constexpr (perf) {
    PERFCTR_INCREMENT(swapCount);
  }
  std::swap(v1, v2);
}

/// @brief check if a vector is sorted
template <template <typename> class Vector2, typename T, typename Compare>
bool IsSorted(Vector2<T>& v, Compare cmp) {
  bool ret = true;
  for (uint64_t i = 1; i < v.size(); i++) {
    ret = ret * (cmp(v[i - 1], v[i]) + (!cmp(v[i], v[i - 1])));
  }
  return ret;
}

/// @brief Wrapper for flex-way distribution o-sort
/// @tparam T type of elements
template <typename T>
  requires(IS_POD<T>())
struct Block {
#if defined(__AVX512VL__) || defined(__AVX2__)
  static constexpr size_t paddingSize = sizeof(T) % 32 == 16 ? 8 : 0;
#else
  static constexpr size_t paddingSize = 0;
#endif
  T data;
  uint32_t tag;               // tie breaker
  bool dummyFlag;             // a flag to mark if the element is dummy
  bool lessFlag;              // a flag to mark if the comparison result is less
  char padding[paddingSize];  // additional padding to make the swap faster
  auto operator==(const Block& other) const { return data == other.data; }
  bool operator<(const Block& other) const {
    return (data < other.data) | ((data == other.data) & (tag < other.tag));
  }
  inline void setData(const T& _data) {
    data = _data;
    tag = UniformRandom32();
    dummyFlag = false;
  }

  inline void setData(const T& _data, uint32_t i) {
    data = _data;
    tag = i;
    dummyFlag = false;
  }

  inline const T& getData() const { return data; }

  static consteval inline Block DUMMY() {
    return Block{T::DUMMY(), 0, true, false};
  }
  inline bool isDummy() const { return dummyFlag; }
  inline bool setAndGetMarked(const Block& pivot) {
    return lessFlag = !(pivot < *this);
  }
  inline bool isMarked(const Block& unused) const { return lessFlag; }

  inline bool isLess() const { return lessFlag; }
  inline void setLessFlag(bool flag) { this->lessFlag = flag; }
  inline void condChangeMark(bool cond, const Block& unused) {
    this->lessFlag ^= cond;
  }
  inline void setDummy() { setDummyFlag(true); }
  inline void setDummyFlag(bool flag) { this->dummyFlag = flag; }
  inline void setDummyFlagCond(bool cond, bool flag) {
    CMOV(cond, this->dummyFlag, flag);
  }
};  // struct Block

/// @brief Wrapper for flex-way butterfly o-sort
/// @tparam T type of elements
template <typename T>
struct TaggedT {
#if defined(__AVX512VL__) || defined(__AVX2__)
  static constexpr size_t paddingSize = sizeof(T) % 32 == 16 ? 8 : 0;
#else
  static constexpr size_t paddingSize = 0;
#endif
  uint64_t tag;  // random label except that the most significant bit is the
                 // flag to mark if the element is dummy
  T v;
  char padding[paddingSize];

  inline void setData(const T& _data) {
    v = _data;
    tag = UniformRandom() & 0x7fff'ffff'ffff'ffffUL;
  }

  inline const T& getData() const { return v; }

  inline bool isDummy() const { return tag >> 63; }

  inline void setDummy() { tag |= 0x8000'0000'0000'0000UL; }

  inline void setTag(uint64_t _tag) { tag = _tag & 0x7fff'ffff'ffff'ffffUL; }

  inline bool setAndGetMarked(uint64_t bitMask) const {  // same as isMarked
    return isMarked(bitMask);
  }

  inline bool isMarked(uint64_t bitMask) const { return !(tag & bitMask); }

  inline void condChangeMark(bool cond, uint64_t bitMask) {
    CMOV(cond, tag, tag ^ bitMask);
  }

  inline uint8_t getMarkAndUpdate(uint64_t k) {
    uint64_t realTag = tag & 0x7fff'ffff'ffff'ffffUL;
    tag &= 0x8000'0000'0000'0000UL;
    tag |= realTag / k;
    uint8_t mark = realTag % k;
    return mark;
  }
};
}  // namespace EM::Algorithm

/// @brief Example of a sort element
struct SortElement {
  uint64_t key;  // key for comparison
  char payload[ELEMENT_SIZE -
               sizeof(key)];  // a payload that is typically longer
  static consteval inline SortElement DUMMY() {
    return SortElement{static_cast<uint64_t>(-1)};
  }
  bool operator==(const SortElement& other) const { return key == other.key; }
  bool operator<(const SortElement& other) const { return key < other.key; }
#ifndef ENCLAVE_MODE
  friend std::ostream& operator<<(std::ostream& o, const SortElement& x) {
    o << x.key;
    return o;
  }
#endif
};

INLINE void CMOV(const uint64_t& condition, SortElement& A,
                 const SortElement& B) {
  CMOV(condition, A.key, B.key);
  for (uint64_t i = 0; i < sizeof(SortElement::payload); i += 8) {
    CMOV(condition, *(uint64_t*)&(A.payload[i]), *(uint64_t*)&(B.payload[i]));
  }
}

template <typename T>
INLINE void CMOV(const uint64_t& condition, EM::Algorithm::Block<T>& A,
                 const EM::Algorithm::Block<T>& B) {
  CMOV(condition, A.data, B.data);
  CMOV(condition, A.tag, B.tag);
  CMOV(condition, A.dummyFlag, B.dummyFlag);
  CMOV(condition, A.lessFlag, B.lessFlag);
}

template <typename T>
INLINE void CMOV(const uint64_t& condition, EM::Algorithm::TaggedT<T>& A,
                 const EM::Algorithm::TaggedT<T>& B) {
  CMOV(condition, A.tag, B.tag);
  CMOV(condition, A.v, B.v);
}

OVERLOAD_TSET_CXCHG(EM::Algorithm::Block<T>, typename T)
OVERLOAD_TSET_CXCHG(EM::Algorithm::TaggedT<T>, typename T);