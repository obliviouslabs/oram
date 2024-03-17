#pragma once
#include <omp.h>

#include <functional>

#include "common/cpp_extended.hpp"
#include "external_memory/stdvector.hpp"
#include "oram/block.hpp"

/**
 * @brief A simple linear oram efficient at small size.
 * Requires uid < size.
 *
 */
namespace ODSL::LinearORAM {
template <typename T, typename UidType = uint64_t>
struct LinearORAM {
  std::vector<T> data;
  UidType _size;

  LinearORAM() : _size(0), data(0) {}
  LinearORAM(UidType size) : _size(size), data(size) {}

  template <typename Reader>
  void InitFromReader(Reader reader) {
    UidType initSize = reader.size();
    for (UidType i = 0; i < initSize; i++) {
      data[i] = reader.read();
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

  static size_t GetMemoryUsage(size_t size) { return sizeof(T) * size; }

  size_t GetMemoryUsage() const { return GetMemoryUsage(_size); }

  void Read(const UidType& uid, T& out) {
    for (UidType i = 0; i < _size; i++) {
      obliMove(i == uid, out, data[i]);
    }
  }

  void Write(const UidType& uid, const T& in) {
    for (UidType i = 0; i < _size; i++) {
      obliMove(i == uid, data[i], in);
    }
  }

  template <class Func>
  void Update(const UidType& uid, const Func& updateFunc) {
    T out;
    Update(uid, updateFunc, out);
  }

  template <class Func>
  void Update(const UidType& uid, const Func& updateFunc, T& out) {
    Update(uid, updateFunc, out, uid);
  }

  template <class Func>
  void Update(const UidType& uid, const Func& updateFunc, T& out,
              const UidType& updatedUid) {
    Read(uid, out);
    bool keepFlag = updateFunc(out);
    UidType newUid = DUMMY<UidType>();
    obliMove(keepFlag, newUid, updatedUid);
    Write(newUid, out);
  }

  void BatchReadAndRemove(uint64_t batchSize, const UidType* uid, T* out) {
    for (uint64_t i = 0; i < batchSize; i++) {
      Read(uid[i], out[i]);
    }
    // don't actually remove since we will write back
  }

  void BatchWriteBack(uint64_t batchSize, const UidType* uid, const T* in,
                      const std::vector<bool>& keepFlag) {
    for (size_t i = 0; i < batchSize; ++i) {
      UidType newUid = uid[i];
      bool isDup = i > 0 && uid[i] == uid[i - 1];
      obliMove(!keepFlag[i] | isDup, newUid, DUMMY<UidType>());
      Write(newUid, in[i]);
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