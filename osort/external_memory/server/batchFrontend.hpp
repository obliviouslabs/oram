#pragma once
#include "external_memory/server/serverFrontend.hpp"
namespace EM {
namespace MemoryServer {
template <typename T, typename _BackendType = ::EM::Backend::MemServerBackend,
          bool ENCRYPTED = true, bool AUTH = true, bool LATE_INIT = true>
struct NonCachedBatchServerFrontend
    : NonCachedServerFrontendInstance<T, _BackendType, ENCRYPTED, AUTH,
                                      LATE_INIT> {
  typedef NonCachedServerFrontendInstance<T, _BackendType, ENCRYPTED, AUTH,
                                          LATE_INIT>
      _Base;
  using IndexType = typename _Base::IndexType;
  using BackendType = _BackendType;
  static inline constexpr auto sizeOfT = _Base::sizeOfT;

  static const uint16_t bufferCapacity = 16;
  // flag to determine whether the current queue is for read or write
  bool isRead = false;
  uint16_t bufferIdx = 0;
  // the counter for the last done job
  uint64_t doneCounter = 0;
  // the counter for the last received job
  uint64_t jobCounter = 0;
  typename T::Encrypted_t encs[ENCRYPTED ? bufferCapacity : 0];
  // output address of query
  T* outAddrs[bufferCapacity];
  // nounces for query
  typename _Base::nounce_t nounces[bufferCapacity];
  PageSlotState pageStates[bufferCapacity] = {EMPTY_PAGE};

  NonCachedBatchServerFrontend(NonCachedBatchServerFrontend& other)
      : _Base(other) {}

  NonCachedBatchServerFrontend(BackendType& _backend, uint64_t initialSize,
                               const T& _defaultVal)
      : _Base(_backend, initialSize, _defaultVal) {}

  NonCachedBatchServerFrontend(BackendType& _backend, uint64_t initialSize)
      : _Base(_backend, initialSize) {}

  // read index i to out, but does not need result until flush read is called
  // return a counter that can be used to check if the page is read
  uint64_t ReadLazy(const IndexType i, T& out, uint32_t auth_counter) {
    PERFCTR_INCREMENT(readCount);
    if constexpr (LATE_INIT) {
      if (!_Base::modified[i]) {
        out = _Base::defaultVal;
        return doneCounter;
      }
    }
    if (!isRead) {
      flushWrite();  // don't read and write at the same time
    }
    isRead = true;
    if (pageStates[bufferIdx] == DONE_PAGE) {
      if constexpr (ENCRYPTED) {
        if constexpr (AUTH) {
          encs[bufferIdx].Decrypt(*(outAddrs[bufferIdx]),
                                  nounces[bufferIdx].bytes);
        } else {
          encs[bufferIdx].Decrypt(*(outAddrs[bufferIdx]));
        }
      }
      pageStates[bufferIdx] = EMPTY_PAGE;
    } else if (pageStates[bufferIdx] == PENDING_PAGE) {
      flushRead();
      // not enough space, flush everything
    }
    Assert(pageStates[bufferIdx] == EMPTY_PAGE);
    pageStates[bufferIdx] = PENDING_PAGE;
    outAddrs[bufferIdx] = &out;
    if constexpr (AUTH) {
      // needs to record the index to authenticate freshness
      typename _Base::nounce_t& nounceCopy = nounces[bufferIdx];
      nounceCopy = _Base::nounce;
      nounceCopy.identifiers.indexPart ^= i;
      nounceCopy.identifiers.counterPart ^= auth_counter;
    }
    if constexpr (ENCRYPTED) {
      typename T::Encrypted_t& inEnc = encs[bufferIdx];
      _Base::backend.ReadLazy(_Base::slot.base + i * sizeOfT, sizeOfT,
                              reinterpret_cast<uint8_t*>(&inEnc),
                              pageStates[bufferIdx]);

    } else {
      _Base::backend.ReadLazy(_Base::slot.base + i * sizeOfT, sizeOfT,
                              reinterpret_cast<uint8_t*>(&out),
                              pageStates[bufferIdx]);
    }

    bufferIdx = (bufferIdx + 1) % bufferCapacity;
    return ++jobCounter;
  }

  // write in to index i, but does not need result until flush write is called
  // return a counter that can be used to check if the page is written
  uint64_t WriteLazy(const IndexType i, T& in, uint32_t auth_counter) {
    PERFCTR_INCREMENT(writeCount);
    if constexpr (LATE_INIT) {
      _Base::modified[i] = true;
    }
    if (isRead) {
      flushRead();
    }
    isRead = false;
    if (pageStates[bufferIdx] == DONE_PAGE) {
      pageStates[bufferIdx] = EMPTY_PAGE;
      ++doneCounter;
    } else if (pageStates[bufferIdx] == PENDING_PAGE) {
      flushWrite();
      // not enough space, flush everything
    }
    Assert(pageStates[bufferIdx] == EMPTY_PAGE);
    pageStates[bufferIdx] = PENDING_PAGE;
    if constexpr (ENCRYPTED) {
      typename T::Encrypted_t& inEnc = encs[bufferIdx];
      if constexpr (AUTH) {
        typename _Base::nounce_t nounceCopy = _Base::nounce;
        nounceCopy.identifiers.indexPart ^= i;
        nounceCopy.identifiers.counterPart ^= auth_counter;
        inEnc.Encrypt(in, nounceCopy.bytes);
      } else {
        inEnc.Encrypt(in);
      }
      _Base::backend.WriteLazy(_Base::slot.base + i * sizeOfT, sizeOfT,
                               reinterpret_cast<const uint8_t*>(&inEnc),
                               pageStates[bufferIdx]);
    } else {
      _Base::backend.WriteLazy(_Base::slot.base + i * sizeOfT, sizeOfT,
                               reinterpret_cast<const uint8_t*>(&in),
                               pageStates[bufferIdx]);
    }

    bufferIdx = (bufferIdx + 1) % bufferCapacity;
    return ++jobCounter;
  }

  void flushRead() {
    for (uint16_t idx = 0; idx < bufferCapacity; ++idx) {
      if (pageStates[idx] == PENDING_PAGE) {
        _Base::backend.FlushRead();
        Assert(pageStates[idx] == DONE_PAGE);
      }
      Assert(pageStates[idx] != PENDING_PAGE);
      if (pageStates[idx] == DONE_PAGE) {
        if constexpr (AUTH) {
          encs[idx].Decrypt(*(outAddrs[idx]), nounces[idx].bytes);
        } else {
          encs[idx].Decrypt(*outAddrs[idx]);
        }
        pageStates[idx] = EMPTY_PAGE;
      }
    }
    doneCounter = jobCounter;
    // all jobs should be done
  }

  void flushWrite() {
    for (uint16_t idx = 0; idx < bufferCapacity; ++idx) {
      if (pageStates[idx] == PENDING_PAGE) {
        _Base::backend.FlushWrite();
        Assert(pageStates[idx] == DONE_PAGE);
      }
      Assert(pageStates[idx] != PENDING_PAGE);
      if (pageStates[idx] == DONE_PAGE) {
        pageStates[idx] = EMPTY_PAGE;
      }
    }
    doneCounter = jobCounter;
    // all jobs should be done
  }

  bool isJobDone(uint64_t _jobCounter) {
    Assert((int64_t)(jobCounter - _jobCounter) >= 0L);
    return (int64_t)(doneCounter - _jobCounter) >= 0L;
  }
};
}  // namespace MemoryServer
}  // namespace EM