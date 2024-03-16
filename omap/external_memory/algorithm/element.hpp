#pragma once
#include "common/defs.hpp"
#include "common/dummy.hpp"
#include "common/encrypted.hpp"
#include "common/mov_intrinsics.hpp"

/// This file defines some macros and data structures for sorting and shuffling
/// algorithms.

// size of each element in bytes
#ifndef ELEMENT_SIZE
#define ELEMENT_SIZE 128
#endif

namespace EM::Algorithm {
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