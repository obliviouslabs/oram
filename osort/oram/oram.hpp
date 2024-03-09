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
          typename UidType = uint64_t, const bool parallel = false>
struct ORAM {
  using LinearORAM_ = LinearORAM::LinearORAM<T, PositionType, UidType>;
  using ORAM_ =
      CircuitORAM::ORAM<T, 2, 20, PositionType, UidType, 4096, 2, parallel>;
  // using ORAM_ = PathORAM::ORAM<T, 5, 64, PositionType, UidType>;
  LinearORAM_* linearOram = NULL;
  ORAM_* treeOram = NULL;
  UidType nextUid = 0;  // uid 0 is reserved for dummy
  Lock _lock;
  bool isLinear = false;
  static constexpr PositionType linear_oram_threshold = 100;

  ORAM() {}

  ORAM(PositionType size) { SetSize(size); }

  ORAM(PositionType size, size_t cacheBytes) { SetSize(size, cacheBytes); }

  ORAM(PositionType size, size_t cacheBytes, int numThreads) {
    static_assert(parallel, "Non-parallel ORAM too many arguments");
    SetSize(size, cacheBytes, numThreads);
  }

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

  void lock() { _lock.lock(); }
  void unlock() { _lock.unlock(); }

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
    if constexpr (parallel) {
      if (numThreads == 0) {
        numThreads = omp_get_max_threads();
      }
    }
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
      if constexpr (parallel) {
        treeOram = new ORAM_(size, cacheBytes, numThreads);
      } else {
        treeOram = new ORAM_(size, cacheBytes);
      }
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
  template <class Func>
  PositionType Update(PositionType pos, const UidType& uid,
                      const Func& updateFunc) {
    if (isLinear) {
      return linearOram->Update(pos, uid, updateFunc);
    } else {
      return treeOram->Update(pos, uid, updateFunc);
    }
  }

  // allow read the updated value
  template <class Func>
  PositionType Update(PositionType pos, const UidType& uid,
                      const Func& updateFunc, T& out) {
    if (isLinear) {
      return linearOram->Update(pos, uid, updateFunc, out);
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
      return linearOram->Update(pos, uid, updateFunc, out, updatedUid);
    } else {
      return treeOram->Update(pos, uid, updateFunc, out, updatedUid);
    }
  }

  // update a batch of elements atomically, return the new positions
  std::vector<PositionType> BatchUpdate(
      const std::vector<PositionType>& pos, const std::vector<UidType>& uid,
      std::function<std::vector<bool>(std::vector<T>&)> updateFunc) {
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
    if (isLinear) {
      linearOram->BatchUpdate(uid, updateFunc, out);
    } else {
      treeOram->BatchUpdate(pos, uid, newPos, updateFunc, out);
    }
  }

  PositionType Read(PositionType pos, const UidType& uid, T& out,
                    PositionType newPos) {
    if (isLinear) {
      return linearOram->Read(pos, uid, out);
    } else {
      return treeOram->Read(pos, uid, out, newPos);
    }
  }

  PositionType Write(const UidType& uid, const T& in, PositionType newPos) {
    if (isLinear) {
      return linearOram->Write(uid, in);
    } else {
      return treeOram->Write(uid, in, newPos);
    }
  }

  template <class Func>
  PositionType Update(PositionType pos, const UidType& uid, PositionType newPos,
                      const Func& updateFunc) {
    if (isLinear) {
      return linearOram->Update(pos, uid, updateFunc);
    } else {
      return treeOram->Update(pos, uid, newPos, updateFunc);
    }
  }

  template <class Func>
  PositionType Update(PositionType pos, const UidType& uid, PositionType newPos,
                      const Func& updateFunc, T& out) {
    if (isLinear) {
      return linearOram->Update(pos, uid, updateFunc, out);
    } else {
      return treeOram->Update(pos, uid, newPos, updateFunc, out);
    }
  }

  template <class Func>
  PositionType Update(PositionType pos, const UidType& uid, PositionType newPos,
                      const Func& updateFunc, T& out,
                      const UidType& updatedUid) {
    if (isLinear) {
      return linearOram->Update(pos, uid, updateFunc, out, updatedUid);
    } else {
      return treeOram->Update(pos, uid, newPos, updateFunc, out, updatedUid);
    }
  }

  template <class Func>
  void ParBatchUpdate(std::vector<PositionType>& pos,
                      const std::vector<UidType>& uid,
                      const std::vector<PositionType>& newPos,
                      const Func& updateFunc, std::vector<T>& out,
                      int numThreads = 0) {
    if (isLinear) {
      linearOram->ParBatchUpdate(uid, updateFunc, out, numThreads);
    } else {
      treeOram->ParBatchUpdate(pos, uid, newPos, updateFunc, out, numThreads);
    }
  }

  template <class Func>
  void ParBatchUpdate(std::vector<PositionType>& pos, std::vector<UidType>& uid,
                      const std::vector<PositionType>& newPos,
                      const Func& updateFunc, int numThreads = 0) {
    if (isLinear) {
      linearOram->ParBatchUpdate(uid, updateFunc, numThreads);
    } else {
      treeOram->ParBatchUpdate(pos, uid, newPos, updateFunc, numThreads);
    }
  }

  template <class Func>
  std::vector<PositionType> ParBatchUpdate(std::vector<PositionType>& pos,
                                           std::vector<UidType>& uid,
                                           const Func& updateFunc,
                                           int numThreads = 0) {
    if (isLinear) {
      linearOram->ParBatchUpdate(pos, uid, updateFunc, numThreads);
      return std::vector<PositionType>(uid.size(), 0);
    } else {
      return treeOram->ParBatchUpdate(uid, updateFunc, numThreads);
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