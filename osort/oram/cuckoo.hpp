#pragma once
#include "recursive_oram.hpp"
namespace ODSL {

template <typename PositionType = uint64_t>
struct CuckooHashMapIndexer {
  static constexpr int saltLength = 16;
  uint8_t salts[2][saltLength];
  PositionType _size;
  CuckooHashMapIndexer(PositionType size) : _size(size) {
    for (int i = 0; i < 2; ++i) {
      read_rand(salts[i], saltLength);
    }
  }
};
}  // namespace ODSL