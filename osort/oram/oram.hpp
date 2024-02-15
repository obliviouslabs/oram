#pragma once

#include "circuit_oram.hpp"
#include "common/lock_util.hpp"
#include "linear_oram.hpp"
#include "path_oram.hpp"

/**
 * @brief An ORAM manager that switches between linear and tree ORAM based on
 * size.
 *
 */
namespace ODSL {
template <typename T, typename PositionType = uint64_t,
          typename UidType = uint64_t>
struct ORAM {
  using LinearORAM_ = LinearORAM::LinearORAM<T, PositionType, UidType>;
  using ORAM_ = CircuitORAM::ORAM<T, 2, 50, PositionType, UidType>;
  // using ORAM_ = PathORAM::ORAM<T, 5, 64, PositionType, UidType>;
  LinearORAM_* linearOram = NULL;
  ORAM_* treeOram = NULL;
  UidType nextUid = 0;  // uid 0 is reserved for dummy
  Lock lock;
  bool isLinear = false;
  static constexpr PositionType linear_oram_threshold = 100;

  ORAM() {}

  ORAM(PositionType size) { SetSize(size); }

  ORAM(PositionType size, size_t cacheBytes) { SetSize(size, cacheBytes); }

  // TODO add copy constructor

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
    Critical section(lock);
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

  void SetSize(PositionType size, size_t cacheBytes = 1UL << 62) {
    Critical section(lock);
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
    Critical section(lock);
    if (isLinear) {
      return linearOram->Read(pos, uid, out);
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
    Critical section(lock);
    if (isLinear) {
      return linearOram->Write(uid, in);
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
  PositionType Update(PositionType pos, const UidType& uid,
                      std::function<bool(T&)> updateFunc) {
    Critical section(lock);
    if (isLinear) {
      return linearOram->Update(pos, uid, updateFunc);
    } else {
      return treeOram->Update(pos, uid, updateFunc);
    }
  }

  // allow read the updated value
  PositionType Update(PositionType pos, const UidType& uid,
                      std::function<bool(T&)> updateFunc, T& out) {
    Critical section(lock);
    if (isLinear) {
      return linearOram->Update(pos, uid, updateFunc, out);
    } else {
      return treeOram->Update(pos, uid, updateFunc, out);
    }
  }

  // allow update the uid
  PositionType Update(PositionType pos, const UidType& uid,
                      std::function<bool(T&)> updateFunc, T& out,
                      const UidType& updatedUid) {
    Critical section(lock);
    if (isLinear) {
      return linearOram->Update(pos, uid, updateFunc, out, updatedUid);
    } else {
      return treeOram->Update(pos, uid, updateFunc, out, updatedUid);
    }
  }

  // update a batch of elements atomically, return the new positions
  std::vector<PositionType> BatchUpdate(
      const std::vector<PositionType>& pos, const std::vector<UidType>& uid,
      std::function<std::vector<bool>(std::vector<T>&)> updateFunc) {
    Critical section(lock);
    if (isLinear) {
      return linearOram->BatchUpdate(uid, updateFunc);
    } else {
      return treeOram->BatchUpdate(pos, uid, updateFunc);
    }
  }

  void BatchUpdate(
      const std::vector<PositionType>& pos, const std::vector<UidType>& uid,
      const std::vector<PositionType>& newPos,
      std::function<std::vector<bool>(std::vector<T>&)> updateFunc) {
    Critical section(lock);
    if (isLinear) {
      linearOram->BatchUpdate(uid, updateFunc);
    } else {
      treeOram->BatchUpdate(pos, uid, newPos, updateFunc);
    }
  }

  void BatchUpdate(const std::vector<PositionType>& pos,
                   const std::vector<UidType>& uid,
                   const std::vector<PositionType>& newPos,
                   std::function<std::vector<bool>(std::vector<T>&)> updateFunc,
                   std::vector<T>& out) {
    Critical section(lock);
    if (isLinear) {
      linearOram->BatchUpdate(uid, updateFunc, out);
    } else {
      treeOram->BatchUpdate(pos, uid, newPos, updateFunc, out);
    }
  }

  PositionType Read(PositionType pos, const UidType& uid, T& out,
                    PositionType newPos) {
    Critical section(lock);
    if (isLinear) {
      return linearOram->Read(pos, uid, out);
    } else {
      return treeOram->Read(pos, uid, out, newPos);
    }
  }

  PositionType Write(const UidType& uid, const T& in, PositionType newPos) {
    Critical section(lock);
    if (isLinear) {
      return linearOram->Write(uid, in);
    } else {
      return treeOram->Write(uid, in, newPos);
    }
  }

  PositionType Update(PositionType pos, const UidType& uid, PositionType newPos,
                      std::function<bool(T&)> updateFunc) {
    Critical section(lock);
    if (isLinear) {
      return linearOram->Update(pos, uid, updateFunc);
    } else {
      return treeOram->Update(pos, uid, newPos, updateFunc);
    }
  }

  PositionType Update(PositionType pos, const UidType& uid, PositionType newPos,
                      std::function<bool(T&)> updateFunc, T& out) {
    Critical section(lock);
    if (isLinear) {
      return linearOram->Update(pos, uid, updateFunc, out);
    } else {
      return treeOram->Update(pos, uid, newPos, updateFunc, out);
    }
  }

  PositionType Update(PositionType pos, const UidType& uid, PositionType newPos,
                      std::function<bool(T&)> updateFunc, T& out,
                      const UidType& updatedUid) {
    Critical section(lock);
    if (isLinear) {
      return linearOram->Update(pos, uid, updateFunc, out, updatedUid);
    } else {
      return treeOram->Update(pos, uid, newPos, updateFunc, out, updatedUid);
    }
  }

  // returns the next unique id, if real is false, returns the dummy id
  UidType GetNextUid(bool real = true) {
    Critical section(lock);
    UidType res = DUMMY<UidType>();
    obliMove(real, res, nextUid);
    nextUid += (UidType)real;
    return res;
  }
};

}  // namespace ODSL