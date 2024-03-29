#pragma once
#include <immintrin.h>
#include <inttypes.h>

#include <cstring>
#include <typeinfo>
#include <utility>

#include "cpp_extended.hpp"

template <const uint64_t sz>
INLINE bool CE_internal(const uint8_t* a, const uint8_t* b) {
  static_assert(sz <= 8);
  if constexpr (sz == 8) {
    return *reinterpret_cast<const uint64_t*>(a) ==
           *reinterpret_cast<const uint64_t*>(b);
  } else if constexpr (sz == 4) {
    return *reinterpret_cast<const uint32_t*>(a) ==
           *reinterpret_cast<const uint32_t*>(b);
  } else if constexpr (sz == 2) {
    return *reinterpret_cast<const uint16_t*>(a) ==
           *reinterpret_cast<const uint16_t*>(b);
  } else if constexpr (sz == 1) {
    return *reinterpret_cast<const uint8_t*>(a) ==
           *reinterpret_cast<const uint8_t*>(b);
  } else {
    uint64_t aCopy = 0;
    uint64_t bCopy = 0;
    memcpy(&aCopy, a, sz);
    memcpy(&bCopy, b, sz);
    return aCopy == bCopy;
  }
}

template <const uint64_t sz>
INLINE bool CL_internal(const uint8_t* a, const uint8_t* b) {
  static_assert(sz <= 8);
  if constexpr (sz == 8) {
    return *reinterpret_cast<const uint64_t*>(a) <
           *reinterpret_cast<const uint64_t*>(b);
  } else if constexpr (sz == 4) {
    return *reinterpret_cast<const uint32_t*>(a) <
           *reinterpret_cast<const uint32_t*>(b);
  } else if constexpr (sz == 2) {
    return *reinterpret_cast<const uint16_t*>(a) <
           *reinterpret_cast<const uint16_t*>(b);
  } else if constexpr (sz == 1) {
    return *reinterpret_cast<const uint8_t*>(a) <
           *reinterpret_cast<const uint8_t*>(b);
  } else {
    uint64_t aCopy = 0;
    uint64_t bCopy = 0;
    memcpy(&aCopy, a, sz);
    memcpy(&bCopy, b, sz);
    return aCopy < bCopy;
  }
}

template <const uint64_t sz>
INLINE bool obliCheckEqual(const uint8_t* a, const uint8_t* b) {
  bool res = true;

  for (uint64_t i = 0; i + 8 <= sz; i += 8) {
    res &= CE_internal<8>(a + i, b + i);
  }
  static constexpr uint64_t rem = sz % 8;
  if constexpr (rem) {
    static constexpr uint64_t last = sz - rem;
    res &= CE_internal<rem>(a + last, b + last);
  }
  return res;
}

// Check if a is less than b, if the system is big-endian, the result is
// lexicographical order, otherwise, it only guarantees to follows a
// deterministic ordering.
template <const uint64_t sz>
INLINE bool obliCheckLess(const uint8_t* a, const uint8_t* b) {
  bool eq = true;
  bool less = false;

  for (uint64_t i = 0; i + 8 < sz; i += 8) {
    less |= eq & CL_internal<8>(a + i, b + i);
    eq &= CE_internal<8>(a + i, b + i);
  }
  static constexpr uint64_t rem = (sz - 1) % 8 + 1;
  static constexpr uint64_t last = sz - rem;
  return less | (eq & CL_internal<rem>(a + last, b + last));
}