#pragma once

#include <cinttypes>
#include <vector>

#include "common/lock_util.hpp"

namespace EM {
// LargeBlockAllocator
// This is a slow allocator used only in oblivious structures constructors.
// This is not optimized for speed, it is meant to not be used very frequently
// to connect oblivious structure memory managers to the external server
// singleton.
//
struct LargeBlockAllocator {
  typedef uint64_t Size_t;
  struct AllocatorSlot {
    Size_t base;
    Size_t size;
  };

  Size_t size;
  Lock lock;

  // $Invariant: IsSorted(freeList)
  //
  std::vector<AllocatorSlot> freeList;

  explicit LargeBlockAllocator(Size_t _size)
      : size(_size), freeList({AllocatorSlot{.base = 0, .size = _size}}) {
    freeList.reserve(1000);
  }

  const AllocatorSlot Allocate(Size_t _size) {
    Critical section(lock);
    TRACE_FUNCTION(_size);

    uint64_t bestIdx = -1;
    uint64_t bestSize = -1;
    for (uint64_t i = 0; i < freeList.size(); i++) {
      if (freeList[i].size >= _size) {
        if (freeList[i].size < bestSize) {
          bestSize = freeList[i].size;
          bestIdx = i;
          if (bestSize == _size) break;
        }
      }
    }

    // Assert(bestIdx != (uint64_t)-1);
    if (bestIdx == (uint64_t)-1) {
      throw std::runtime_error("LargeBlockAllocator: Out of memory");
    }
    AllocatorSlot ret = freeList[bestIdx];
    ret.size = _size;
    freeList[bestIdx].size -= _size;
    freeList[bestIdx].base += _size;
    if (freeList[bestIdx].size == 0) {
      freeList.erase(freeList.begin() + bestIdx);
    }
    return ret;
  }

  void Free(const AllocatorSlot& slot) {
    Critical section(lock);
    TRACE_FUNCTION(slot.base, slot.size);

    if (slot.size == 0) return;

    AllocatorSlot tmp = slot;
    // In place insert:
    //
    for (uint64_t i = 0; i < freeList.size(); i++) {
      if (freeList[i].base > tmp.base) {
        std::swap(freeList[i], tmp);
      }
    }
    freeList.push_back(tmp);

    // Now we should coallesce the blocks:
    //
    uint64_t j = 0;
    for (uint64_t i = 1; i < freeList.size(); i++) {
      if (freeList[j].base + freeList[j].size == freeList[i].base) {
        freeList[j].size += freeList[i].size;
      } else {
        j++;
        if (i != j) {
          freeList[j] = freeList[i];
        }
      }
    }
    freeList.resize(j + 1);
  }
};
}  // namespace EM