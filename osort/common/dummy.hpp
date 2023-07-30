#pragma once
#include "common/cpp_extended.hpp"

template <typename T>
INLINE consteval T DUMMY() {
  return T::DUMMY();
}

template <>
INLINE consteval uint8_t DUMMY<uint8_t>() {
  return static_cast<uint8_t>(-1);
}

template <>
INLINE consteval uint16_t DUMMY<uint16_t>() {
  return static_cast<uint16_t>(-1);
}

template <>
INLINE consteval uint32_t DUMMY<uint32_t>() {
  return static_cast<uint32_t>(-1);
}

template <>
INLINE consteval uint64_t DUMMY<uint64_t>() {
  return static_cast<uint64_t>(-1);
}