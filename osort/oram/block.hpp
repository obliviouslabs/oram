#pragma once
#include "common/cpp_extended.hpp"
#include "common/dummy.hpp"
#include "common/mov_intrinsics.hpp"
namespace ORAM::PathORAM {
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

};  // namespace ORAM::PathORAM