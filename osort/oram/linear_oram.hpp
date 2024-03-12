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

  std::vector<PositionType> BatchUpdate(
      const std::vector<UidType>& uid,
      std::function<std::vector<bool>(std::vector<T>&)> updateFunc) {
    std::vector<T> out(uid.size());
    BatchUpdate(uid, updateFunc, out);
    return std::vector<PositionType>(uid.size(), 0);
  }

  void BatchUpdate(const std::vector<UidType>& uid,
                   std::function<std::vector<bool>(std::vector<T>&)> updateFunc,
                   std::vector<T>& out) {
    uint64_t batchSize = uid.size();
    Assert(batchSize == out.size());

    for (uint64_t i = 0; i < batchSize; i++) {
      Read(0, uid[i], out[i]);
    }

    const std::vector<bool>& writeBackFlags = updateFunc(out);

    for (uint64_t i = 0; i < batchSize; i++) {
      for (UidBlock_& block : data) {
        bool match = block.uid == uid[i];
        obliMove(match, block.data, out[i]);
        obliMove(match & !writeBackFlags[i], block.uid, DUMMY<UidType>());
      }
    }
  }

  void BatchReadAndRemove(const std::vector<UidType>& uid,
                          std::vector<T>& out) {
    for (uint64_t i = 0; i < uid.size(); i++) {
      Read(0, uid[i], out[i]);
    }
    // don't actually remove since we will write back
  }

  void BatchWriteBack(const std::vector<UidType>& uid, const std::vector<T>& in,
                      const std::vector<bool>& keepFlag) {
    for (size_t i = 0; i < uid.size(); ++i) {
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
  void BatchUpdateWithDup(const std::vector<UidType>& uid,
                          const Func& updateFunc, std::vector<T>& out) {
    // TODO: Implement hash table index / compaction
    // printf("Using linear oram\n");
    BatchReadAndRemove(uid, out);

    std::vector<bool> keepFlag = updateFunc(out);

    BatchWriteBack(uid, out, keepFlag);
    // deduplicate uids
  }

  // require uid to be sorted, can change uid
  // require oram to be sorted and contain keys 0 to size - 1
  template <class Func>
  void ParBatchUpdate(std::vector<UidType>& uid, const Func& updateFunc,
                      std::vector<T>& out, int numThreads = 0) {
    // TODO: Implement hash table index / compaction
    if (numThreads == 0) {
      numThreads = omp_get_max_threads();
    }
    if (true || size() < 2 * GetLogBaseTwo(uid.size())) {
      // printf("Using linear oram\n");
#pragma omp parallel for num_threads(numThreads)
      for (uint64_t i = 0; i < uid.size(); i++) {
        Read(0, uid[i], out[i]);
      }

      std::vector<bool> keepFlag = updateFunc(out);
      // printf("Update done %lu, numthreads = %d\n", uid.size(), numThreads);
      // deduplicate uids
      for (size_t i = uid.size() - 1; i > 0; --i) {
        obliMove(uid[i] == uid[i - 1], uid[i], DUMMY<UidType>());
      }
      int chunkSize = divRoundUp(uid.size(), numThreads);

// avoid write contention and false sharing
#pragma omp parallel for num_threads(numThreads) schedule(static, chunkSize)
      for (UidBlock_& block : data) {
        for (uint64_t i = 0; i < uid.size(); i++) {
          bool match = block.uid == uid[i];
          obliMove(match, block.data, out[i]);
          obliMove(match & !keepFlag[i], block.uid, DUMMY<UidType>());
        }
      }
      // printf("Write done\n");
      return;
    }
    // std::vector<uint32_t> queryPrefixSum(uid.size() + 1, 0);
    // queryPrefixSum[0] = 0;
    // for (uint64_t i = 0; i < uid.size(); i++) {
    //   queryPrefixSum[i + 1] = queryPrefixSum[i] + (uid[i] !=
    //   DUMMY<UidType>());
    // }
  }

  template <class Func>
  void ParBatchUpdate(std::vector<UidType>& uid, const Func& updateFunc,
                      int numThreads = 0) {
    std::vector<T> out(uid.size());
    ParBatchUpdate(uid, updateFunc, out, numThreads);
  }

  template <class Func>
  void BatchUpdateWithDup(const std::vector<UidType>& uid,
                          const Func& updateFunc) {
    std::vector<T> out(uid.size());
    BatchUpdateWithDup(uid, updateFunc, out);
  }
};
};  // namespace ODSL::LinearORAM
    // namespace ODSL::LinearORAM