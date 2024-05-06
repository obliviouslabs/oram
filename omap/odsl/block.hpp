#pragma once
#include "common/cpp_extended.hpp"
#include "common/dummy.hpp"
#include "common/mov_intrinsics.hpp"
/// @brief This file contains the definition of the blocks used in oram.

namespace ODSL {

/// @brief A block is a wrapper around a data element, with a position and a
/// unique id.
/// @tparam T The type of the data
/// @tparam PositionType The type of the position, default to uint64_t
/// @tparam UidType The type of the unique id, default to uint64_t
template <typename T, typename PositionType = uint64_t,
          typename UidType = uint64_t>
struct Block {
  T data;
  PositionType position;
  UidType uid = DUMMY<UidType>();
  Block() = default;
  Block(const T& data, PositionType position, UidType uid)
      : data(data), position(position), uid(uid) {}

  Block(PositionType position, UidType uid)
      : uid(uid), position(position), data(DUMMY<T>()) {}
  Block(const Block& other) = default;
  Block(Block&& other) = default;
  Block& operator=(const Block& other) = default;
  Block& operator=(Block&& other) = default;
  bool operator==(const Block& other) const { return uid == other.uid; }
  bool operator!=(const Block& other) const { return !(*this == other); }
  static_assert(std::is_trivially_copyable_v<T>,
                "T must be trivially copyable");

  INLINE bool ReadAndRemove(const UidType& uid, T& out) {
    bool match = this->uid == uid;
    obliMove(match, this->uid, DUMMY<UidType>());
    obliMove(match, out, data);
    return match;
  }

  INLINE bool Insert(bool cond, const Block& data) {
    return obliMove(cond & IsDummy(), *this, data);
  }

  INLINE bool IsDummy() const { return uid == DUMMY<UidType>(); }

#ifndef ENCLAVE_MODE
  friend std::ostream& operator<<(std::ostream& os, const Block& block) {
    os << "Block{data: " << block.data << ", position: " << block.position
       << ", uid: ";
    if (block.IsDummy()) {
      os << "-1";
    } else {
      os << block.uid;
    }
    os << "}";
    return os;
  }
#endif
};

/// @brief A block with a unique id, but no position. Useful for linear oram.
/// @tparam T The type of the data
/// @tparam UidType The type of the unique id, default to uint64_t
template <typename T, typename UidType = uint64_t>
struct UidBlock {
  T data;
  UidType uid = DUMMY<UidType>();
  UidBlock() = default;
  UidBlock(const T& data, UidType uid) : data(data), uid(uid) {}

  explicit UidBlock(UidType uid) : uid(uid), data(DUMMY<T>()) {}
  UidBlock(const UidBlock& other) = default;
  UidBlock(UidBlock&& other) = default;
  UidBlock& operator=(const UidBlock& other) = default;
  UidBlock& operator=(UidBlock&& other) = default;
  bool operator==(const UidBlock& other) const { return uid == other.uid; }
  bool operator!=(const UidBlock& other) const { return !(*this == other); }

  // define less
  bool operator<(const UidBlock& other) const { return uid < other.uid; }

  INLINE bool IsDummy() const { return uid == DUMMY<UidType>(); }
};

};  // namespace ODSL