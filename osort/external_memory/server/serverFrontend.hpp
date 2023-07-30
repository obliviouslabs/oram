#pragma once
#include <cinttypes>
#include <vector>

#include "common/dmcache.hpp"
#include "common/encrypted.hpp"
#include "common/lrucache.hpp"
#include "common/tracing/tracer.hpp"
#include "common/utils.hpp"
#include "external_memory/server/cached.hpp"
#include "external_memory/server/serverAllocator.hpp"
#include "external_memory/server/serverBackend.hpp"

namespace EM {
namespace MemoryServer {

template <typename T, typename _BackendType = ::EM::Backend::MemServerBackend,
          bool ENCRYPTED = true, bool AUTH = true, bool LATE_INIT = false>
// #ifndef ENCLAVE_MODE
// requires EM::Backend::BackendServer<BackendType>
// #endif
struct NonCachedServerFrontendInstance {
  // Just forwards reads and writes to the server. Call the allocator during
  // construction and resizes.
  //

  using AllocatorSlot = typename EM::LargeBlockAllocator::AllocatorSlot;
  using BackendType = _BackendType;

  typedef uint64_t IndexType;
  BackendType& backend;

  static inline constexpr auto sizeOfT =
      ENCRYPTED ? sizeof(typename T::Encrypted_t) : sizeof(T);

  T defaultVal;

  std::vector<bool> modified;

  AllocatorSlot slot;

  typedef union {
    uint8_t bytes[IV_SIZE];
    struct {
      IndexType indexPart;
      uint32_t counterPart;  // in case an index is written multiple times
      uint32_t padding;
    } identifiers;
  } nounce_t;

  nounce_t nounce;
  // prevent the slot to be freed, the slot should be saved elsewhere to avoid
  // memory leak
  void preventFree() { slot.base = -1; }

  NonCachedServerFrontendInstance(NonCachedServerFrontendInstance& other)
      : backend(other.backend) {
    defaultVal = other.defaultVal;
    slot = other.slot;
    if constexpr (AUTH) {
      nounce = other.nounce;
    }
    other.slot.base = -1;
    if constexpr (LATE_INIT) {
      std::swap(modified, other.modified);
    }
  }

  NonCachedServerFrontendInstance(BackendType& _backend, uint64_t initialSize,
                                  const T& _defaultVal)
      : backend(_backend) {
    IndexType requiredSize = initialSize * sizeOfT;
    slot = backend.Allocate(requiredSize);
    // std::cout << "Alloc: " << slot.base << "--" << slot.base + slot.size <<
    // std::endl;
    if constexpr (AUTH) {
      nounce.identifiers.indexPart = UniformRandom();
      nounce.identifiers.counterPart = UniformRandom32();
    }
    if constexpr (LATE_INIT) {
      this->defaultVal = _defaultVal;
      modified.resize(initialSize, false);
    } else {
      for (uint64_t i = 0; i < initialSize; i++) {
        Write(i, _defaultVal);
      }
    }
  }

  NonCachedServerFrontendInstance(BackendType& _backend, uint64_t initialSize)
      : backend(_backend) {
    IndexType requiredSize = initialSize * sizeOfT;
    slot = backend.Allocate(requiredSize);
    if constexpr (AUTH) {
      nounce.identifiers.indexPart = UniformRandom();
      nounce.identifiers.counterPart = UniformRandom32();
    }
  }

  ~NonCachedServerFrontendInstance() {
    if (slot.base == -1) {
      return;
    }
    // std::cout << "Freed: " << slot.base << "--" << slot.base + slot.size <<
    // std::endl;
    backend.Free(slot);
    slot.size = 0;
  }

  void Write(const IndexType i, const T& in) {
    if constexpr (AUTH) {
      Write(i, in, 0);
    } else {
      if constexpr (LATE_INIT) {
        modified[i] = true;
      }
      PERFCTR_INCREMENT(writeCount);
      if constexpr (ENCRYPTED) {
        typename T::Encrypted_t inEnc;
        inEnc.Encrypt(in);

        backend.Write(slot.base + i * sizeOfT, sizeOfT,
                      reinterpret_cast<const uint8_t*>(&inEnc));
      } else {
        backend.Write(slot.base + i * sizeOfT, sizeOfT,
                      reinterpret_cast<const uint8_t*>(&in));
      }
    }
  }

  void Write(const IndexType i, const T& in, uint32_t counter) {
    static_assert(AUTH);
    if constexpr (LATE_INIT) {
      modified[i] = true;
    }
    PERFCTR_INCREMENT(writeCount);

    typename T::Encrypted_t inEnc;
    nounce_t nounceCopy = nounce;
    nounceCopy.identifiers.indexPart ^= i;
    nounceCopy.identifiers.counterPart ^= counter;
    inEnc.Encrypt(in, nounceCopy.bytes);
    backend.Write(slot.base + i * sizeOfT, sizeOfT,
                  reinterpret_cast<const uint8_t*>(&inEnc));
  }

  // dummy implementation to be overloaded
  uint64_t WriteLazy(const IndexType i, const T& in) {
    Write(i, in);
    return 0;
  }

  uint64_t WriteLazy(const IndexType i, const T& in, uint32_t counter) {
    Write(i, in, counter);
    return 0;
  }

  void Read(const IndexType i, T& out) {
    if constexpr (AUTH) {
      Read(i, out, 0);
    } else {
      if constexpr (LATE_INIT) {
        if (!modified[i]) {
          out = defaultVal;
          return;
        }
      }
      PERFCTR_INCREMENT(readCount);

      if constexpr (ENCRYPTED) {
        typename T::Encrypted_t inEnc;
        backend.Read(slot.base + i * sizeOfT, sizeOfT,
                     reinterpret_cast<uint8_t*>(&inEnc));
        inEnc.Decrypt(out);
      } else {
        backend.Read(slot.base + i * sizeOfT, sizeOfT,
                     reinterpret_cast<uint8_t*>(&out));
      }
    }
  }

  void Read(const IndexType i, T& out, uint32_t counter) {
    static_assert(AUTH);
    if constexpr (LATE_INIT) {
      if (!modified[i]) {
        out = defaultVal;
        return;
      }
    }
    PERFCTR_INCREMENT(readCount);
    nounce_t nounceCopy = nounce;
    nounceCopy.identifiers.indexPart ^= i;
    nounceCopy.identifiers.counterPart ^= counter;
    typename T::Encrypted_t inEnc;
    backend.Read(slot.base + i * sizeOfT, sizeOfT,
                 reinterpret_cast<uint8_t*>(&inEnc));
    inEnc.Decrypt(out, nounceCopy.bytes);
  }

  uint64_t ReadLazy(const IndexType i, T& out) {
    Read(i, out);
    return 0;
  }

  uint64_t ReadLazy(const IndexType i, T& out, uint32_t counter) {
    Read(i, out, counter);
    return 0;
  }

  void flushRead() {}

  void flushRead(uint32_t counter) {}

  void flushWrite() {}

  void flushWrite(uint32_t counter) {}
};

template <typename T, typename BackendType = ::EM::Backend::MemServerBackend,
          bool ENCRYPTED = true, bool AUTH = true,
          uint64_t CACHE_SIZE = SERVER__CACHE_SIZE, uint64_t TLB_SIZE = 2>
using ServerFrontendInstance = CACHE::Cached<
    T, NonCachedServerFrontendInstance<T, BackendType, ENCRYPTED, AUTH, false>,
    CACHE_SIZE, TLB_SIZE>;

}  // namespace MemoryServer
}  // namespace EM