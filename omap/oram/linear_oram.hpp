#pragma once
#include <omp.h>

#include <functional>

#include "common/cpp_extended.hpp"
#include "external_memory/stdvector.hpp"
#include "oram/block.hpp"

/**
 * @brief A simple linear oram efficient at small size.
 *
 */
namespace ODSL::LinearORAM {
template <typename T, typename PositionType = uint64_t,
          typename UidType = uint64_t>
struct LinearORAM {
  using UidBlock_ = UidBlock<T, UidType>;
  StdVector<UidBlock_> data;
  PositionType _size;

  LinearORAM() : _size(0), data(0) {}
  LinearORAM(PositionType size) : _size(size), data(size) {}

  template <typename Reader>
  void InitFromReader(Reader reader) {
    UidType initSize = reader.size();
    for (UidType i = 0; i < initSize; i++) {
      data[i] = UidBlock_(reader.read(), i);
    }
    for (UidType i = initSize; i < _size; i++) {
      data[i].uid = DUMMY<UidType>();
    }
  }

  template <typename Reader, typename Writer>
  void InitFromReader(Reader reader, Writer writer) {
    InitFromReader(reader);
  }

  size_t size() const { return _size; }

  void SetSize(size_t size) {
    _size = size;
    data.resize(size);
  }

  static size_t GetMemoryUsage(size_t size) { return sizeof(UidBlock_) * size; }

  size_t GetMemoryUsage() const { return GetMemoryUsage(_size); }

  PositionType Read(PositionType pos, const UidType& uid, T& out) {
    for (UidBlock_& block : data) {
      obliMove(block.uid == uid, out, block.data);
    }
    return pos;
  }

  PositionType Write(const UidType& uid, const T& in) {
    bool inserted = false;
    for (UidBlock_& block : data) {
      bool emptyFlag = block.uid == DUMMY<UidType>();
      UidBlock_ newBlock(in, uid);
      obliMove(emptyFlag & !inserted, block, newBlock);
      inserted |= emptyFlag;
    }
    return 0;
  }

  template <class Func>
  PositionType Update(PositionType pos, const UidType& uid,
                      const Func& updateFunc) {
    T out;
    return Update(pos, uid, updateFunc, out);
  }

  template <class Func>
  PositionType Update(PositionType pos, const UidType& uid,
                      const Func& updateFunc, T& out) {
    return Update(pos, uid, updateFunc, out, uid);
  }

  template <class Func>
  PositionType Update(PositionType pos, const UidType& uid,
                      const Func& updateFunc, T& out,
                      const UidType& updatedUid) {
    for (const UidBlock_& block : data) {
      bool match = block.uid == uid;
      obliMove(match, out, block.data);
    }
    bool keepFlag = updateFunc(out);
    UidType newUid = DUMMY<UidType>();
    obliMove(keepFlag, newUid, updatedUid);
    UidBlock_ newBlock(out, newUid);
    for (UidBlock_& block : data) {
      bool match = block.uid == uid;
      obliMove(match, block, newBlock);
    }
    return pos;
  }

  void BatchReadAndRemove(uint64_t batchSize, const UidType* uid,
                          T* out) {
    for (uint64_t i = 0; i < batchSize; i++) {
      Read(0, uid[i], out[i]);
    }
    // don't actually remove since we will write back
  }

  void BatchWriteBack(uint64_t batchSize, const UidType* uid, const T* in,
                      const std::vector<bool>& keepFlag) {
    for (size_t i = 0; i < batchSize; ++i) {
      UidType searchUid = uid[i];
      if (i != 0) {
        obliMove(uid[i] == uid[i - 1], searchUid, DUMMY<UidType>());
        // don't update if same uid
      }
      UidBlock_ blockToWrite(in[i], searchUid);
      obliMove(!keepFlag[i], blockToWrite.uid, DUMMY<UidType>());
      for (UidBlock_& block : data) {
        bool match = block.uid == searchUid;
        obliMove(match, block, blockToWrite);
      }
    }
  }

  template <class Func>
  void BatchUpdate(uint64_t batchSize, const UidType* uid,
                          const Func& updateFunc, T* out) {
    BatchReadAndRemove(batchSize, uid, out);

    std::vector<bool> keepFlag = updateFunc(batchSize, out);

    BatchWriteBack(batchSize, uid, out, keepFlag);
    // deduplicate uids
  }

  template <class Func>
  void BatchUpdate(uint64_t batchSize, const UidType* uid,
                          const Func& updateFunc) {
    std::vector<T> out(batchSize);
    BatchUpdate(batchSize, uid, updateFunc, &out[0]);
  }
};
};  // namespace ODSL::LinearORAM
    // namespace ODSL::LinearORAM