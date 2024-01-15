#pragma once
#include "block.hpp"

namespace ORAM::PathORAM {
template <typename T, const int Z = 4, typename PositionType = uint64_t,
          typename UidType = uint64_t>
struct Bucket {
  using BlockType = Block<T, PositionType, UidType>;
  BlockType blocks[Z];
  Bucket() = default;
  Bucket(const Bucket& other) = default;
  Bucket(Bucket&& other) = default;
  Bucket& operator=(const Bucket& other) = default;
  Bucket& operator=(Bucket&& other) = default;

  // Oblivious query the bucket by uid, return true if found, false otherwise.
  bool queryByUid(const UidType& uid, T& out) {
    bool res = false;
    for (int i = 0; i < sizeof(T); i++) {
      res |= blocks[i].invalidateAndCopyDataIfUidMatch(uid, out);
    }
    return res;
  }

  bool conditionalWriteToDummySlot(bool cond, const BlockType& block) {
    for (int i = 0; i < Z; i++) {
      cond &= !blocks[i].conditionalFillIfDummy(cond, block);
    }
    return cond;
  }
};
};  // namespace ORAM::PathORAM