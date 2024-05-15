#pragma once
#include "block.hpp"

/// @brief This file contains the definition of bucket used in oram.

namespace ODSL {

/// @brief A bucket is a wrapper around a set of blocks.
template <typename T, const int Z, typename PositionType = uint64_t,
          typename UidType = uint64_t>
struct Bucket {
  using BlockType = Block<T, PositionType, UidType>;
  static_assert(std::is_trivially_copyable_v<T>,
                "T must be trivially copyable");
  static_assert(std::is_trivially_copyable_v<BlockType>,
                "Block must be trivially copyable");
  BlockType blocks[Z];
  Bucket() = default;
  Bucket(const Bucket& other) = default;
  Bucket(Bucket&& other) = default;
  Bucket& operator=(const Bucket& other) = default;
  Bucket& operator=(Bucket&& other) = default;
  static inline Bucket DUMMY() {
    Bucket b;
    for (int i = 0; i < Z; i++) {
      b.blocks[i] = BlockType();
    }
    return b;
  }

#ifndef ENCLAVE_MODE
  // cout
  friend std::ostream& operator<<(std::ostream& os, const Bucket& b) {
    os << "[ ";
    for (int i = 0; i < Z; i++) {
      os << b.blocks[i] << " ";
    }
    os << "]";
    return os;
  }
#endif
};

template <typename BucketType>
struct FreshORAMNode {
  BucketType bucket;
  uint64_t leftNonce = 0;
  uint64_t rightNonce = 0;
  static inline FreshORAMNode DUMMY() {
    FreshORAMNode node;
    node.bucket = BucketType::DUMMY();
    return node;
  }
};

template <typename BucketType>
struct NonFreshORAMNode {
  BucketType bucket;
  static inline NonFreshORAMNode DUMMY() {
    NonFreshORAMNode node;
    node.bucket = BucketType::DUMMY();
    return node;
  }
};
};  // namespace ODSL