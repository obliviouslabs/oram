#pragma once
#include "common/cpp_extended.hpp"
#include "common/dummy.hpp"
#include "common/mov_intrinsics.hpp"
namespace ORAM {
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

  INLINE bool checkUidMatch(const UidType& uid) const {
    return this->uid == uid;
  }

  INLINE bool invalidateAndCopyDataIfUidMatch(const UidType& uid, T& out) {
    bool match = this->uid == uid;
    obliMove(match, this->uid, DUMMY<UidType>());
    obliMove(match, out, data);
    return match;
  }

  INLINE bool conditionalFillIfDummy(bool cond, const Block& data) {
    return obliMove(cond & isDummy(), *this, data);
  }

  INLINE bool isDummy() const { return uid == DUMMY<UidType>(); }
};

template <typename T, typename UidType = uint64_t>
struct UidBlock {
  T data;
  UidType uid = DUMMY<UidType>();
  UidBlock() = default;
  UidBlock(const T& data, UidType uid) : data(data), uid(uid) {}

  UidBlock(UidType uid) : uid(uid), data(DUMMY<T>()) {}
  UidBlock(const UidBlock& other) = default;
  UidBlock(UidBlock&& other) = default;
  UidBlock& operator=(const UidBlock& other) = default;
  UidBlock& operator=(UidBlock&& other) = default;
  bool operator==(const UidBlock& other) const { return uid == other.uid; }
  bool operator!=(const UidBlock& other) const { return !(*this == other); }

  // define less
  bool operator<(const UidBlock& other) const { return uid < other.uid; }

  INLINE bool checkUidMatch(const UidType& uid) const {
    return this->uid == uid;
  }

  INLINE bool isDummy() const { return uid == DUMMY<UidType>(); }
};

};  // namespace ORAM::PathORAM