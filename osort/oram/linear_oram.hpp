#pragma once
#include <functional>

#include "common/cpp_extended.hpp"
#include "external_memory/stdvector.hpp"
#include "oram/block.hpp"

namespace ORAM::LinearORAM {
template <typename T, typename PositionType = uint64_t,
          typename UidType = uint64_t>
struct LinearORAM {
  using UidBlock_ = UidBlock<T, UidType>;
  StdVector<UidBlock_> data;
  PositionType _size;

  LinearORAM() : _size(0) {}
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

  PositionType Update(PositionType pos, const UidType& uid,
                      std::function<void(T&)> updateFunc) {
    T out;
    return Update(pos, uid, updateFunc, out);
  }

  PositionType Update(PositionType pos, const UidType& uid,
                      std::function<void(T&)> updateFunc, T& out) {
    return Update(pos, uid, updateFunc, out, uid);
  }

  PositionType Update(PositionType pos, const UidType& uid,
                      std::function<void(T&)> updateFunc, T& out,
                      const UidType& updatedUid) {
    for (const UidBlock_& block : data) {
      bool match = block.uid == uid;
      obliMove(match, out, block.data);
    }
    updateFunc(out);
    UidBlock_ newBlock(out, updatedUid);
    for (UidBlock_& block : data) {
      bool match = block.uid == uid;
      obliMove(match, block, newBlock);
    }
    return pos;
  }
};
};  // namespace ORAM::LinearORAM
    // namespace ORAM::LinearORAM