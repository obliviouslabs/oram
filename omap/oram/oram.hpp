#pragma once

#include "circuit_oram.hpp"
#include "common/lock_util.hpp"
#include "linear_oram.hpp"

/**
 * @brief An ORAM manager that switches between linear and tree ORAM based on
 * size.
 *
 */
namespace ODSL {
template <typename T, typename PositionType = uint64_t,
          typename UidType = uint64_t>
struct ORAM {
  using LinearORAM_ = LinearORAM::LinearORAM<T, UidType>;
  using ORAM_ = CircuitORAM::ORAM<T, 2, 20, PositionType, UidType, 4096, 2>;
  LinearORAM_* linearOram = NULL;
  ORAM_* treeOram = NULL;
  UidType nextUid = 0;  // uid 0 is reserved for dummy
  bool isLinear = false;
  static constexpr PositionType linear_oram_threshold = 100;

  ORAM() {}

  ORAM(PositionType size) { SetSize(size); }

  ORAM(PositionType size, size_t cacheBytes) { SetSize(size, cacheBytes); }

  ~ORAM() {
    if (linearOram) {
      delete linearOram;
      linearOram = NULL;
    }
    if (treeOram) {
      delete treeOram;
      treeOram = NULL;
    }
  }

  /**
   * @brief Initialize the oram from an existing ram array.
   * Reader must be sorted and has uid [0, reader.size())
   * Writer will write out the position map for tree oram as (uid, pos) pairs.
   *
   * @tparam Reader
   * @tparam Writer
   * @param reader
   * @param writer
   */
  template <typename Reader, typename Writer>
  void InitFromReader(Reader& reader, Writer& writer) {
    nextUid = reader.size();
    if (isLinear) {
      linearOram->InitFromReader(reader);
    } else {
      treeOram->InitFromReader(reader, writer);
    }
  }

  // TODO rename to capacity
  PositionType size() const {
    if (isLinear) {
      Assert(linearOram);
      return linearOram->size();
    } else {
      Assert(treeOram);
      return treeOram->size();
    }
  }

  void SetSize(PositionType size, size_t cacheBytes = 1UL << 62,
               int numThreads = 0) {
    if (linearOram || treeOram) {
      throw std::runtime_error("SetSize can only be called on empty oram");
    }
    isLinear = size <= linear_oram_threshold;
    if (isLinear) {
      if (LinearORAM_::GetMemoryUsage(size) > cacheBytes) {
        throw std::runtime_error(
            "LinearORAM_::GetMemoryUsage(size) > cacheBytes");
      }
      linearOram = new LinearORAM_(size);
    } else {
      treeOram = new ORAM_(size, cacheBytes);
    }
  }

  static size_t GetMemoryUsage(size_t size, size_t cacheBytes = 1UL << 62) {
    if (size <= linear_oram_threshold) {
      return LinearORAM_::GetMemoryUsage(size);
    } else {
      return ORAM_::GetMemoryUsage(size, cacheBytes);
    }
  }

  size_t GetMemoryUsage() const {
    if (isLinear) {
      Assert(linearOram);
      return linearOram->GetMemoryUsage();
    } else {
      Assert(treeOram);
      return treeOram->GetMemoryUsage();
    }
  }

  /**
   * @brief Read from the oram.
   *
   * @param pos position (path idx) of the element
   * @param uid uid of the element
   * @param out output
   * @return PositionType new position of the element
   */
  PositionType Read(PositionType pos, const UidType& uid, T& out) {
    if (isLinear) {
      linearOram->Read(uid, out);
      return 0;
    } else {
      return treeOram->Read(pos, uid, out);
    }
  }

  /**
   * @brief Write an element to the oram.
   *
   * @param uid uid of the element to write
   * @param in element to write
   * @return PositionType position of the element
   */
  PositionType Write(const UidType& uid, const T& in) {
    if (isLinear) {
      linearOram->Write(uid, in);
      return 0;
    } else {
      return treeOram->Write(uid, in);
    }
  }

  /**
   * @brief Update the element at pos with uid.
   *
   * @param pos position (path idx) of the element
   * @param uid uid of the element
   * @param updateFunc may modify the element in place, return true to keep the
   * element, false to delete it.
   * @return PositionType new position of the element
   */
  template <class Func>
  PositionType Update(PositionType pos, const UidType& uid,
                      const Func& updateFunc) {
    if (isLinear) {
      linearOram->Update(uid, updateFunc);
      return 0;
    } else {
      return treeOram->Update(pos, uid, updateFunc);
    }
  }

  // allow read the updated value
  template <class Func>
  PositionType Update(PositionType pos, const UidType& uid,
                      const Func& updateFunc, T& out) {
    if (isLinear) {
      linearOram->Update(uid, updateFunc, out);
      return 0;
    } else {
      return treeOram->Update(pos, uid, updateFunc, out);
    }
  }

  // allow update the uid
  template <class Func>
  PositionType Update(PositionType pos, const UidType& uid,
                      const Func& updateFunc, T& out,
                      const UidType& updatedUid) {
    if (isLinear) {
      linearOram->Update(uid, updateFunc, out, updatedUid);
      return 0;
    } else {
      return treeOram->Update(pos, uid, updateFunc, out, updatedUid);
    }
  }

  PositionType Read(PositionType pos, const UidType& uid, T& out,
                    PositionType newPos) {
    if (isLinear) {
      linearOram->Read(uid, out);
      return 0;
    } else {
      return treeOram->Read(pos, uid, out, newPos);
    }
  }

  PositionType Write(const UidType& uid, const T& in, PositionType newPos) {
    if (isLinear) {
      linearOram->Write(uid, in);
      return 0;
    } else {
      return treeOram->Write(uid, in, newPos);
    }
  }

  template <class Func>
  PositionType Update(PositionType pos, const UidType& uid, PositionType newPos,
                      const Func& updateFunc) {
    if (isLinear) {
      linearOram->Update(uid, updateFunc);
      return 0;
    } else {
      return treeOram->Update(pos, uid, newPos, updateFunc);
    }
  }

  template <class Func>
  PositionType Update(PositionType pos, const UidType& uid, PositionType newPos,
                      const Func& updateFunc, T& out) {
    if (isLinear) {
      linearOram->Update(uid, updateFunc, out);
      return 0;
    } else {
      return treeOram->Update(pos, uid, newPos, updateFunc, out);
    }
  }

  template <class Func>
  PositionType Update(PositionType pos, const UidType& uid, PositionType newPos,
                      const Func& updateFunc, T& out,
                      const UidType& updatedUid) {
    if (isLinear) {
      linearOram->Update(uid, updateFunc, out, updatedUid);
      return 0;
    } else {
      return treeOram->Update(pos, uid, newPos, updateFunc, out, updatedUid);
    }
  }

  void BatchReadAndRemove(uint64_t batchSize, PositionType* pos,
                          const UidType* uid, T* out) {
    if (isLinear) {
      linearOram->BatchReadAndRemove(batchSize, uid, out);
    } else {
      treeOram->BatchReadAndRemove(batchSize, pos, uid, out);
    }
  }

  void BatchWriteBack(uint64_t batchSize, UidType* uid,
                      const PositionType* newPos, const T* in,
                      const std::vector<bool>& writeBackFlags) {
    if (isLinear) {
      linearOram->BatchWriteBack(batchSize, uid, in, writeBackFlags);
    } else {
      treeOram->BatchWriteBack(batchSize, uid, newPos, in, writeBackFlags);
    }
  }

  template <class Func>
  void BatchUpdate(uint64_t batchSize, PositionType* pos, const UidType* uid,
                   const PositionType* newPos, const Func& updateFunc, T* out) {
    if (isLinear) {
      linearOram->BatchUpdate(batchSize, uid, updateFunc, out);
    } else {
      treeOram->BatchUpdate(batchSize, pos, uid, newPos, updateFunc, out);
    }
  }

  template <class Func>
  void BatchUpdate(uint64_t batchSize, PositionType* pos, const UidType* uid,
                   const PositionType* newPos, const Func& updateFunc) {
    if (isLinear) {
      linearOram->BatchUpdate(batchSize, uid, updateFunc);
    } else {
      treeOram->BatchUpdate(batchSize, pos, uid, newPos, updateFunc);
    }
  }

  template <class Func>
  std::vector<PositionType> BatchUpdate(uint64_t batchSize, PositionType* pos,
                                        const UidType* uid,
                                        const Func& updateFunc) {
    if (isLinear) {
      linearOram->BatchUpdate(batchSize, uid, updateFunc);
      return std::vector<PositionType>(batchSize, 0);
    } else {
      return treeOram->BatchUpdate(batchSize, pos, uid, updateFunc);
    }
  }

  // returns the next unique id, if real is false, returns the dummy id
  UidType GetNextUid(bool real = true) {
    UidType res = DUMMY<UidType>();
    obliMove(real, res, nextUid);
    nextUid += (UidType)real;
    return res;
  }
};

}  // namespace ODSL