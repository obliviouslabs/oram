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
  std::vector<PositionType> BatchUpdate(uint64_t batchSize,
      const PositionType* pos, const UidType* uid,
      std::function<std::vector<bool>(uint64_t, T*)> updateFunc) {
    if (isLinear) {
      return linearOram->BatchUpdate(batchSize, uid, updateFunc);
    } else {
      return treeOram->BatchUpdate(batchSize, pos, uid, updateFunc);
    }
  }

  void BatchUpdate(
    uint64_t batchSize,
      const PositionType* pos, const UidType* uid,
      const PositionType* newPos,
      std::function<std::vector<bool>(uint64_t, T*)> updateFunc) {
    if (isLinear) {
      linearOram->BatchUpdate(batchSize, uid, updateFunc);
    } else {
      treeOram->BatchUpdate(batchSize, pos, uid, newPos, updateFunc);
    }
  }

  void BatchUpdate(uint64_t batchSize, const PositionType* pos,
                   const UidType* uid,
                   const PositionType* newPos,
                   std::function<std::vector<bool>(uint64_t, T*)> updateFunc,
                   T* out) {
    if (isLinear) {
      linearOram->BatchUpdate(batchSize, uid, updateFunc, out);
    } else {
      treeOram->BatchUpdate(batchSize, pos, uid, newPos, updateFunc, out);
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

  void BatchReadAndRemove(uint64_t batchSize, PositionType* pos,
                          const UidType* uid,
                          T* out) {
    if (isLinear) {
      linearOram->BatchReadAndRemove(batchSize, uid, out);
    } else {
      treeOram->BatchReadAndRemove(batchSize, pos, uid, out);
    }
  }

  void BatchWriteBack(uint64_t batchSize, UidType* uid,
                      const PositionType* newPos,
                      const T* in,
                      const std::vector<bool>& writeBackFlags) {
    if (isLinear) {
      linearOram->BatchWriteBack(batchSize, uid, in, writeBackFlags);
    } else {
      treeOram->BatchWriteBack(batchSize, uid, newPos, in, writeBackFlags);
    }
  }

  template <class Func>
  void BatchUpdateWithDup(uint64_t batchSize, PositionType* pos,
                          const UidType* uid,
                          const PositionType* newPos,
                          const Func& updateFunc, T* out) {
    if (isLinear) {
      linearOram->BatchUpdateWithDup(batchSize, uid, updateFunc, out);
    } else {
      treeOram->BatchUpdateWithDup(batchSize, pos, uid, newPos, updateFunc, out);
    }
  }

  template <class Func>
  void BatchUpdateWithDup(uint64_t batchSize, PositionType* pos,
                          const UidType* uid,
                          const PositionType* newPos,
                          const Func& updateFunc) {
    if (isLinear) {
      linearOram->BatchUpdateWithDup(batchSize, uid, updateFunc);
    } else {
      treeOram->BatchUpdateWithDup(batchSize, pos, uid, newPos, updateFunc);
    }
  }

  template <class Func>
  std::vector<PositionType> BatchUpdateWithDup(uint64_t batchSize, PositionType* pos,
                                               const UidType* uid,
                                               const Func& updateFunc) {
    if (isLinear) {
      linearOram->BatchUpdateWithDup(batchSize, uid, updateFunc);
      return std::vector<PositionType>(batchSize, 0);
    } else {
      return treeOram->BatchUpdateWithDup(batchSize, pos, uid, updateFunc);
    }
  }

  template <class Func>
  void ParBatchUpdate(uint64_t batchSize, PositionType* pos,
                      const UidType* uid,
                      const PositionType* newPos,
                      const Func& updateFunc, T* out,
                      int numThreads = 0) {
    if (isLinear) {
      linearOram->ParBatchUpdate(batchSize, uid, updateFunc, out, numThreads);
    } else {
      treeOram->ParBatchUpdate(batchSize, pos, uid, newPos, updateFunc, out, numThreads);
    }
  }

  template <class Func>
  void ParBatchUpdate(uint64_t batchSize, PositionType* pos, UidType* uid,
                      const PositionType* newPos,
                      const Func& updateFunc, int numThreads = 0) {
    if (isLinear) {
      linearOram->ParBatchUpdate(batchSize, uid, updateFunc, numThreads);
    } else {
      treeOram->ParBatchUpdate(batchSize, pos, uid, newPos, updateFunc, numThreads);
    }
  }

  template <class Func>
  std::vector<PositionType> ParBatchUpdate(uint64_t batchSize, PositionType* pos,
                                           UidType* uid,
                                           const Func& updateFunc,
                                           int numThreads = 0) {
    if (isLinear) {
      linearOram->ParBatchUpdate(batchSize, pos, uid, updateFunc, numThreads);
      return std::vector<PositionType>(uid.size(), 0);
    } else {
      return treeOram->ParBatchUpdate(batchSize, uid, updateFunc, numThreads);
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