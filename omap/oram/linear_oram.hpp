#pragma once
#include <omp.h>

#include <functional>

#include "common/cpp_extended.hpp"
#include "external_memory/algorithm/or_compact_shuffle.hpp"
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

  void BatchRead(uint64_t batchSize, const UidType* uid, T* out) {
    if (std::min(batchSize, _size) > 20) {
      BatchReadViaCompaction(batchSize, uid, out);
    } else {
      BatchReadNaive(batchSize, uid, out);
    }
  }

  void BatchReadNaive(uint64_t batchSize, const UidType* uid, T* out) {
    for (uint64_t i = 0; i < batchSize; i++) {
      Read(uid[i], out[i]);
    }
  }

  void BatchReadViaCompaction(uint64_t batchSize, const UidType* uid, T* out) {
    std::vector<UidType> uidCopy(std::max((uint64_t)_size + 1, batchSize),
                                 DUMMY<UidType>());
    std::vector<UidType> prefixSum(batchSize + 1);
    std::vector<T> oramCopy(std::max((uint64_t)_size, batchSize));
    std::copy(uid, uid + batchSize, uidCopy.begin());
    std::copy(data.begin(), data.end(), oramCopy.begin());
    prefixSum[0] = 0;
    prefixSum[1] = 1;  // the first uid must be kept
    for (uint64_t i = 1; i < batchSize; i++) {
      bool isDup = uid[i] == uid[i - 1];
      uidCopy[i] -= prefixSum[i];
      obliMove(isDup, uidCopy[i], DUMMY<UidType>());
      prefixSum[i + 1] = prefixSum[i] + !isDup;
    }
    EM::Algorithm::OrCompactSeparateMark(
        uidCopy.begin(), uidCopy.begin() + batchSize, prefixSum.begin());

    int distributeLevel = GetLogBaseTwo(_size - 1);
    for (uint64_t mask = 1UL << distributeLevel; mask != 0; mask >>= 1) {
      for (int64_t i = _size - mask - 1; i >= 0; --i) {
        bool movFlag =
            (!!(uidCopy[i] & mask)) & (uidCopy[i] != DUMMY<UidType>());
        obliSwap(movFlag, uidCopy[i + mask], uidCopy[i]);
      }
    }

    UidType oramPrefixSum = 0;
    for (uint64_t i = 0; i < _size; ++i) {
      bool isReal = uidCopy[i] != DUMMY<UidType>();
      uidCopy[i] = oramPrefixSum;
      oramPrefixSum += isReal;
    }
    uidCopy[_size] = oramPrefixSum;
    EM::Algorithm::OrCompactSeparateMark(
        oramCopy.begin(), oramCopy.begin() + _size, uidCopy.begin());

    EM::Algorithm::OrDistributeSeparateMark(
        oramCopy.begin(), oramCopy.begin() + batchSize, prefixSum.begin());
    for (uint64_t i = 0; i < batchSize; i++) {
      if (i > 0) {
        bool isDup = uid[i] == uid[i - 1];
        obliMove(isDup, oramCopy[i],
                 oramCopy[i - 1]);  // copy data for duplicate value
      }
      out[i] = oramCopy[i];
    }
  }

  void BatchWriteBack(uint64_t batchSize, const UidType* uid, const T* in,
                      const std::vector<bool>& keepFlag) {
    if (std::min(batchSize, _size) > 20) {
      BatchWriteBackViaCompaction(batchSize, uid, in, keepFlag);
    } else {
      BatchWriteBackNaive(batchSize, uid, in, keepFlag);
    }
  }

  void BatchWriteBackNaive(uint64_t batchSize, const UidType* uid, const T* in,
                           const std::vector<bool>& keepFlag) {
    for (size_t i = 0; i < batchSize; ++i) {
      UidType newUid = uid[i];
      bool isDup = i > 0 && uid[i] == uid[i - 1];
      obliMove(!keepFlag[i] | isDup, newUid, DUMMY<UidType>());
      Write(newUid, in[i]);
    }
  }

  void BatchWriteBackViaCompaction(uint64_t batchSize, const UidType* uid,
                                   const T* in,
                                   const std::vector<bool>& keepFlag) {
    std::vector<UidType> uidCopy(std::max((uint64_t)_size + 1, batchSize),
                                 DUMMY<UidType>());
    std::vector<UidType> prefixSum(batchSize + 1);
    std::vector<T> inCopy(std::max((uint64_t)_size, batchSize));
    std::copy(in, in + batchSize, inCopy.begin());
    prefixSum[0] = 0;
    prefixSum[1] = 1;  // the first uid must be kept
    for (size_t i = 0; i < batchSize; ++i) {
      UidType newUid = uid[i];
      bool isDup = (i > 0 && uid[i] == uid[i - 1]) | !keepFlag[i];
      newUid -= prefixSum[i];
      obliMove(isDup, newUid, DUMMY<UidType>());
      uidCopy[i] = newUid;
      prefixSum[i + 1] = prefixSum[i] + !isDup;
    }

    EM::Algorithm::OrCompactSeparateMark(
        uidCopy.begin(), uidCopy.begin() + batchSize, prefixSum.begin());

    int distributeLevel = GetLogBaseTwo(_size - 1);
    for (uint64_t mask = 1UL << distributeLevel; mask != 0; mask >>= 1) {
      for (int64_t i = _size - mask - 1; i >= 0; --i) {
        bool movFlag =
            (!!(uidCopy[i] & mask)) & (uidCopy[i] != DUMMY<UidType>());
        obliSwap(movFlag, uidCopy[i + mask], uidCopy[i]);
      }
    }

    UidType oramPrefixSum = 0;
    for (uint64_t i = 0; i < _size; ++i) {
      bool isReal = uidCopy[i] != DUMMY<UidType>();
      uidCopy[i] = oramPrefixSum;
      oramPrefixSum += isReal;
    }
    uidCopy[_size] = oramPrefixSum;
    // EM::Algorithm::OrCompactSeparateMark(
    //     oramCopy.begin(), oramCopy.begin() + _size, uidCopy.begin());
    EM::Algorithm::OrCompactSeparateMark(
        inCopy.begin(), inCopy.begin() + batchSize, prefixSum.begin());

    EM::Algorithm::OrDistributeSeparateMark(
        inCopy.begin(), inCopy.begin() + _size, uidCopy.begin());
    for (uint64_t i = 0; i < _size; i++) {
      bool isReal = uidCopy[i] != uidCopy[i + 1];
      obliMove(isReal, data[i], inCopy[i]);
    }
  }

  template <class Func>
  void BatchUpdate(uint64_t batchSize, const UidType* uid,
                   const Func& updateFunc, T* out) {
    BatchRead(batchSize, uid, out);

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