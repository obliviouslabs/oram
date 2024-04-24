#pragma once
#include "common/cmp_intrinsics.hpp"
#include "common/defs.hpp"
#include "common/dummy.hpp"
#include "common/mov_intrinsics.hpp"
#include "common/utils.hpp"

/// This file defines some macros and data structures for sorting and shuffling
/// algorithms.

// size of each element in bytes
#ifndef ELEMENT_SIZE
#define ELEMENT_SIZE 128
#endif

namespace Algorithm {
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

  inline bool IsDummy() const { return tag >> 63; }

  inline void setDummy() { tag |= 0x8000'0000'0000'0000UL; }

  inline void setTag(uint64_t _tag) { tag = _tag & 0x7fff'ffff'ffff'ffffUL; }

  inline bool setAndGetMarked(uint64_t bitMask) const {  // same as isMarked
    return isMarked(bitMask);
  }

  inline bool isMarked(uint64_t bitMask) const { return !(tag & bitMask); }

  inline void condChangeMark(bool cond, uint64_t bitMask) {
    obliMove(cond, tag, tag ^ bitMask);
  }

  inline uint8_t getMarkAndUpdate(uint64_t k) {
    uint64_t realTag = tag & 0x7fff'ffff'ffff'ffffUL;
    tag &= 0x8000'0000'0000'0000UL;
    tag |= realTag / k;
    uint8_t mark = realTag % k;
    return mark;
  }
};
}  // namespace Algorithm

/// @brief Example of a sort element
struct TestElement {
  uint64_t key;  // key for comparison
  uint8_t payload[ELEMENT_SIZE -
                  sizeof(key)];  // a payload that is typically longer
  static consteval inline TestElement DUMMY() {
    return TestElement{static_cast<uint64_t>(-1), {}};
  }
  bool operator==(const TestElement& other) const { return key == other.key; }
  bool operator<(const TestElement& other) const { return key < other.key; }
#ifndef ENCLAVE_MODE
  friend std::ostream& operator<<(std::ostream& o, const TestElement& x) {
    o << x.key;
    return o;
  }
#endif
};

template <const size_t size>
struct Bytes {
  uint8_t data[size];

  Bytes() = default;
  
  explicit Bytes(const void* _data) { memcpy(data, _data, size); }

  bool operator==(const Bytes<size>& other) const {
    return obliCheckEqual<size>(data, other.data);
  }

  bool operator!=(const Bytes<size>& other) const { return !(*this == other); }

  bool operator<(const Bytes<size>& other) const {
    return obliCheckLess<size>(data, other.data);
  }

  void SetRand() { read_rand(data, size); }

  const uint8_t* GetData() const { return data; }

  static consteval inline Bytes DUMMY() {
    Bytes ret;
    for (size_t i = 0; i < size; ++i) {
      ret.data[i] = 0xff;
    }
    return ret;
  }

// out stream
#ifndef ENCLAVE_MODE
  friend std::ostream& operator<<(std::ostream& o, const Bytes<size>& x) {
    for (size_t i = 0; i < size; ++i) {
      o << std::hex << std::setw(2) << std::setfill('0') << (int)x.data[i]
        << std::dec;
    }
    return o;
  }
#endif
};

namespace std {
template <const size_t size>
struct hash<Bytes<size>> {
  std::size_t operator()(const Bytes<size>& bytes) const {
    return std::hash<std::string_view>()(
        std::string_view((const char*)bytes.GetData(), size));
  }
};
}  // namespace std