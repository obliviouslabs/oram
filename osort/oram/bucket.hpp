#pragma once
#include "block.hpp"

namespace ODSL {
template <typename T, const int Z = 4, typename PositionType = uint64_t,
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
};  // namespace ODSL