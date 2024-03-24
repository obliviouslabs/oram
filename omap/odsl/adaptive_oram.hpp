#pragma once

#include "circuit_oram.hpp"
#include "common/lock_util.hpp"
#include "linear_oram.hpp"

/**
 * @brief An ORAM manager that switches between linear and tree ORAM based on
 * size.
 *
 */
namespace ODSL::AdaptiveORAM {
template <typename T, typename PositionType = uint64_t,
          typename UidType = uint64_t>
struct ORAM {
  using LinearORAM_ = LinearORAM::ORAM<T, UidType>;
  using ORAM_ = CircuitORAM::ORAM<T, 2, 20, PositionType, UidType, 4096, 2>;
  LinearORAM_* linearOram = NULL;
  ORAM_* treeOram = NULL;
  UidType nextUid = 0;  // uid 0 is reserved for dummy
  bool isLinear = false;
  static constexpr PositionType linear_oram_threshold = 500;

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
   *
   * @tparam Reader The type of the reader
   * @tparam Writer The type of the writer
   * @param reader must be sorted and has uid [0, reader.size())
   * @param writer write out the position map for tree oram as (uid, pos) pairs.
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
   * @brief Read a block from the ORAM and assign it to a new position
   *
   * @param pos The current position of the block
   * @param uid The unique id of the block
   * @param out The output block
   * @param newPos The new position of the block
   * @return PositionType The new position of the block
   */
  PositionType Read(PositionType pos, const UidType& uid, T& out,
                    PositionType newPos) {
    if (isLinear) {
      linearOram->Read(uid, out);
      return 0;
    } else {
      return treeOram->Read(pos, uid, out, newPos);
    }
  }

  /**
   * @brief Write a new block to the ORAM
   *
   * @param uid The unique id of the block
   * @param in The new block to write
   * @param newPos The new position of the block
   * @return PositionType The new position of the block
   */
  PositionType Write(const UidType& uid, const T& in, PositionType newPos) {
    if (isLinear) {
      linearOram->Write(uid, in);
      return 0;
    } else {
      return treeOram->Write(uid, in, newPos);
    }
  }

  /**
   * @brief Update a block in the ORAM and assign it to a random new position
   *
   * @tparam Func The type of the update function. The update function should
   * take a reference to the block and return a bool. If the return value is
   * true, the block is kept, otherwise it is deleted.
   * @param pos The current position of the block
   * @param uid The unique id of the block
   * @param updateFunc The update function
   * @return PositionType The random new position of the block
   */
  template <class Func>
    requires UpdateOrRemoveFunction<Func, T>
  PositionType Update(PositionType pos, const UidType& uid,
                      const Func& updateFunc) {
    if (isLinear) {
      linearOram->Update(uid, updateFunc);
      return 0;
    } else {
      return treeOram->Update(pos, uid, updateFunc);
    }
  }

  /**
   * @brief Update a block in the ORAM and assign it to a random new position
   *
   * @tparam Func The type of the update function. The update function should
   * take a reference to the block and return a bool. If the return value is
   * true, the block is kept, otherwise it is deleted.
   * @param pos The current position of the block
   * @param uid The unique id of the block
   * @param updateFunc The update function
   * @param out Output the updated block
   * @return PositionType The random new position of the block
   */
  template <class Func>
    requires UpdateOrRemoveFunction<Func, T>
  PositionType Update(PositionType pos, const UidType& uid,
                      const Func& updateFunc, T& out) {
    if (isLinear) {
      linearOram->Update(uid, updateFunc, out);
      return 0;
    } else {
      return treeOram->Update(pos, uid, updateFunc, out);
    }
  }

  /**
   * @brief Update a block in the ORAM and assign it to a random new position.
   * Possibly change the uid.
   *
   * @tparam Func The type of the update function. The update function should
   * take a reference to the block and return a bool. If the return value is
   * true, the block is kept, otherwise it is deleted.
   * @param pos The current position of the block
   * @param uid The unique id of the block
   * @param updateFunc The update function
   * @param out Output the updated block
   * @param updatedUid The new uid of the block
   * @return PositionType The random new position of the block
   */
  template <class Func>
    requires UpdateOrRemoveFunction<Func, T>
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

  /**
   * @brief Update a block in the ORAM and assign it to a new position
   *
   * @tparam Func The type of the update function. The update function should
   * take a reference to the block and return a bool. If the return value is
   * true, the block is kept, otherwise it is deleted.
   * @param pos The current position of the block
   * @param uid The unique id of the block
   * @param newPos The new position of the block
   * @param updateFunc The update function
   * @return PositionType
   */
  template <class Func>
    requires UpdateOrRemoveFunction<Func, T>
  PositionType Update(PositionType pos, const UidType& uid, PositionType newPos,
                      const Func& updateFunc) {
    if (isLinear) {
      linearOram->Update(uid, updateFunc);
      return 0;
    } else {
      return treeOram->Update(pos, uid, newPos, updateFunc);
    }
  }

  /**
   * @brief Update a block in the ORAM and assign it to a new position
   *
   * @tparam Func The type of the update function. The update function should
   * take a reference to the block and return a bool. If the return value is
   * true, the block is kept, otherwise it is deleted.
   * @param pos The current position of the block
   * @param uid The unique id of the block
   * @param newPos The new position of the block
   * @param updateFunc The update function
   * @param out Output the updated block
   * @return PositionType The new position of the block
   */
  template <class Func>
    requires UpdateOrRemoveFunction<Func, T>
  PositionType Update(PositionType pos, const UidType& uid, PositionType newPos,
                      const Func& updateFunc, T& out) {
    if (isLinear) {
      linearOram->Update(uid, updateFunc, out);
      return 0;
    } else {
      return treeOram->Update(pos, uid, newPos, updateFunc, out);
    }
  }

  /**
   * @brief Update a block in the ORAM and assign it to a new position.
   * Possibly change the uid.
   *
   * @tparam Func The type of the update function. The update function should
   * take a reference to the block and return a bool. If the return value is
   * true, the block is kept, otherwise it is deleted.
   * @param pos The current position of the block
   * @param uid The unique id of the block
   * @param newPos The new position of the block
   * @param updateFunc The update function
   * @param out Output the updated block
   * @param updatedUid The new uid of the block
   * @return PositionType The random new position of the block
   */
  template <class Func>
    requires UpdateOrRemoveFunction<Func, T>
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

  /**
   * @brief Read a batch of blocks and remove them from the ORAM.
   *
   * @param batchSize The number of blocks
   * @param pos The positions of the blocks
   * @param uid The uids of the blocks
   * @param out The output blocks
   *
   */
  void BatchReadAndRemove(uint64_t batchSize, PositionType* pos,
                          const UidType* uid, T* out) {
    if (isLinear) {
      linearOram->BatchRead(batchSize, uid, out);
    } else {
      treeOram->BatchReadAndRemove(batchSize, pos, uid, out);
    }
  }

  /**
   * @brief Write back a batch of blocks to the ORAM. For duplicate blocks,
   * ignore all but the first block.
   *
   * @param batchSize The number of blocks
   * @param uid The uids of the blocks. Requires uid to be sorted.
   * @param newPos The new positions of the blocks
   * @param in The new blocks
   * @param writeBackFlags The write back flags. If writeBackFlags[i] is true,
   * the block is written back, otherwise it is not written back.
   */
  void BatchWriteBack(uint64_t batchSize, const UidType* uid,
                      const PositionType* newPos, const T* in,
                      const std::vector<bool>& writeBackFlags) {
    if (isLinear) {
      linearOram->BatchWriteBack(batchSize, uid, in, writeBackFlags);
    } else {
      treeOram->BatchWriteBack(batchSize, uid, newPos, in, writeBackFlags);
    }
  }

  /**
   * @brief Update a batch of blocks in the ORAM and assign them to new
   * positions
   *
   * @tparam Func The type of the update function
   * @param batchSize The number of blocks
   * @param pos The positions of the blocks. May be modified.
   * @param uid The uids of the blocks. Requires uid to be sorted.
   * @param newPos The new positions of the blocks
   * @param updateFunc The update function. The update function should take a
   * integer batchSize and a pointer to an array of batchSize blocks, and return
   * a std::vector<bool> of batchSize, indicating whether each block should be
   * written back.
   * @param out Pointer to an array of output blocks
   */
  template <class Func>
    requires BatchUpdateOrRemoveFunction<Func, T>
  void BatchUpdate(uint64_t batchSize, PositionType* pos, const UidType* uid,
                   const PositionType* newPos, const Func& updateFunc, T* out) {
    if (isLinear) {
      linearOram->BatchUpdate(batchSize, uid, updateFunc, out);
    } else {
      treeOram->BatchUpdate(batchSize, pos, uid, newPos, updateFunc, out);
    }
  }

  /**
   * @brief Update a batch of blocks in the ORAM and assign them to new
   * positions
   *
   * @tparam Func The type of the update function
   * @param batchSize The number of blocks
   * @param pos The positions of the blocks. May be modified.
   * @param uid The uids of the blocks. Requires uid to be sorted.
   * @param newPos The new positions of the blocks
   * @param updateFunc The update function. The update function should take a
   * integer batchSize and a pointer to an array of batchSize blocks, and return
   * a std::vector<bool> of batchSize, indicating whether each block should be
   * written back.
   */
  template <class Func>
    requires BatchUpdateOrRemoveFunction<Func, T>
  void BatchUpdate(uint64_t batchSize, PositionType* pos, const UidType* uid,
                   const PositionType* newPos, const Func& updateFunc) {
    if (isLinear) {
      linearOram->BatchUpdate(batchSize, uid, updateFunc);
    } else {
      treeOram->BatchUpdate(batchSize, pos, uid, newPos, updateFunc);
    }
  }

  /**
   * @brief Update a batch of blocks in the ORAM and assign them to random new
   * positions
   *
   * @tparam Func The type of the update function
   * @param batchSize The number of blocks
   * @param pos The positions of the blocks. May be modified.
   * @param uid The uids of the blocks. Requires uid to be sorted.
   * @param updateFunc The update function. The update function should take a
   * integer batchSize and a pointer to an array of batchSize blocks, and return
   * a std::vector<bool> of batchSize, indicating whether each block should be
   * written back.
   * @return const std::vector<PositionType> The random new positions of the
   * blocks.
   */
  template <class Func>
    requires BatchUpdateOrRemoveFunction<Func, T>
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

}  // namespace ODSL::AdaptiveORAM