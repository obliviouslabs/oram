#pragma once
#include "algorithm/or_compact_shuffle.hpp"
#include "common/cpp_extended.hpp"
#include "odsl/block.hpp"

/**
 * @brief A simple linear oram efficient at small size.
 * Requires uid < size.
 *
 */
namespace ODSL::LinearORAM {
template <typename T, typename UidType = uint64_t>
struct ORAM {
 private:
  UidType _size;        // oram size
  std::vector<T> data;  // oram data, indexed by uid

 public:
  ORAM() : _size(0), data(0) {}
  explicit ORAM(UidType size) : _size(size), data(size) {}

  template <typename Reader>
    requires Readable<Reader, T>
  void InitFromReader(Reader reader) {
    UidType initSize = reader.size();
    for (UidType i = 0; i < initSize; i++) {
      data[i] = reader.read();
    }
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
    requires UpdateOrRemoveFunction<Func, T>
  void Update(const UidType& uid, const Func& updateFunc) {
    T out = T();
    Update(uid, updateFunc, out);
  }

  template <class Func>
    requires UpdateOrRemoveFunction<Func, T>
  void Update(const UidType& uid, const Func& updateFunc, T& out) {
    Update(uid, updateFunc, out, uid);
  }

  template <class Func>
    requires UpdateOrRemoveFunction<Func, T>
  void Update(const UidType& uid, const Func& updateFunc, T& out,
              const UidType& updatedUid) {
    Read(uid, out);
    bool keepFlag = updateFunc(out);
    UidType newUid = DUMMY<UidType>();
    obliMove(keepFlag, newUid, updatedUid);
    Write(newUid, out);
  }

  /**
   * @brief Read multiple uids from the oram, requires uid to be sorted.
   *
   * @param batchSize the number of uids to read
   * @param uid pointer to array of the uids to read
   * @param out pointer to array of the output
   */
  void BatchRead(uint64_t batchSize, const UidType* uid, T* out) {
    if (std::min(batchSize, (uint64_t)_size) > 40) {
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

  /**
   * @brief When both the batch size and the oram size are large, we use
   * compaction and distribution algorithms to read the oram. The time
   * complexity is O(N log N + B log B), where N is the oram size and B is the
   * batch size.
   */
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
      // if not dummy, let uidCopy[i] be the number of distinct uids before it
      uidCopy[i] -= prefixSum[i];
      obliMove(isDup, uidCopy[i], DUMMY<UidType>());
      prefixSum[i + 1] = prefixSum[i] + !isDup;
    }
    // compact the unique uids
    Algorithm::OrCompactSeparateMark(
        uidCopy.begin(), uidCopy.begin() + batchSize, prefixSum.begin());

    // prepare marks for distributing uids to there corresponding location
    int distributeLevel = GetLogBaseTwo(_size - 1);
    for (uint64_t mask = 1UL << distributeLevel; mask != 0; mask >>= 1) {
      for (int64_t i = _size - mask - 1; i >= 0; --i) {
        bool movFlag =
            (!!(uidCopy[i] & mask)) & (uidCopy[i] != DUMMY<UidType>());
        obliSwap(movFlag, uidCopy[i + mask], uidCopy[i]);
      }
    }

    // distribute the uids to their corresponding location
    // let uid store the prefix sum of the number of real uids before it
    UidType oramPrefixSum = 0;
    for (uint64_t i = 0; i < _size; ++i) {
      bool isReal = uidCopy[i] != DUMMY<UidType>();
      uidCopy[i] = oramPrefixSum;
      oramPrefixSum += isReal;
    }
    uidCopy[_size] = oramPrefixSum;

    // compact the oram entries we want to query
    Algorithm::OrCompactSeparateMark(oramCopy.begin(), oramCopy.begin() + _size,
                                     uidCopy.begin());

    // distribute these oram entries to align with the first unique uids in the
    // batch
    Algorithm::OrDistributeSeparateMark(
        oramCopy.begin(), oramCopy.begin() + batchSize, prefixSum.begin());

    // propagate results for duplicates, and copy the results to the output
    for (uint64_t i = 0; i < batchSize; i++) {
      if (i > 0) {
        bool isDup = uid[i] == uid[i - 1];
        obliMove(isDup, oramCopy[i],
                 oramCopy[i - 1]);  // copy data for duplicate value
      }
      out[i] = oramCopy[i];
    }
  }

  /**
   * @brief Write multiple uids to the oram, requires uid to be sorted.
   *
   * @param batchSize the number of uids to write
   * @param uid pointer to array of the uids to write
   * @param in pointer to array of the input
   * @param keepFlag pointer to array of the keep flags, if keepFlag[i] is
   * false, the i-th uid will be a dummy write
   */
  void BatchWriteBack(uint64_t batchSize, const UidType* uid, const T* in,
                      const std::vector<bool>& keepFlag) {
    if (std::min(batchSize, (uint64_t)_size) > 40) {
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
      obliMove((!keepFlag[i]) | isDup, newUid, DUMMY<UidType>());
      Write(newUid, in[i]);
    }
  }

  /**
   * @brief When both the batch size and the oram size are large, we use
   * compaction and distribution algorithms to write the oram. The time
   * complexity is O(N log N + B log B), where N is the oram size and B is the
   * batch size. Similar to the BatchReadViaCompaction function.
   */
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
      bool isDup = (i > 0 && uid[i] == uid[i - 1]) | (!keepFlag[i]);
      newUid -= prefixSum[i];
      obliMove(isDup, newUid, DUMMY<UidType>());
      uidCopy[i] = newUid;
      prefixSum[i + 1] = prefixSum[i] + !isDup;
    }

    Algorithm::OrCompactSeparateMark(
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
    Algorithm::OrCompactSeparateMark(inCopy.begin(), inCopy.begin() + batchSize,
                                     prefixSum.begin());

    Algorithm::OrDistributeSeparateMark(inCopy.begin(), inCopy.begin() + _size,
                                        uidCopy.begin());
    for (uint64_t i = 0; i < _size; i++) {
      bool isReal = uidCopy[i] != uidCopy[i + 1];
      obliMove(isReal, data[i], inCopy[i]);
    }
  }

  /**
   * @brief Update multiple uids in the oram, requires uid to be sorted.
   *
   * @tparam Func the type of the update function
   * @param batchSize the number of uids to update
   * @param uid pointer to array of the uids to update
   * @param updateFunc the update function
   * @param out pointer to array of the output
   */
  template <class Func>
    requires BatchUpdateOrRemoveFunction<Func, T>
  void BatchUpdate(uint64_t batchSize, const UidType* uid,
                   const Func& updateFunc, T* out) {
    BatchRead(batchSize, uid, out);

    std::vector<bool> keepFlag = updateFunc(batchSize, out);

    BatchWriteBack(batchSize, uid, out, keepFlag);
  }

  template <class Func>
    requires BatchUpdateOrRemoveFunction<Func, T>
  void BatchUpdate(uint64_t batchSize, const UidType* uid,
                   const Func& updateFunc) {
    std::vector<T> out(batchSize);
    BatchUpdate(batchSize, uid, updateFunc, &out[0]);
  }
};
};  // namespace ODSL::LinearORAM