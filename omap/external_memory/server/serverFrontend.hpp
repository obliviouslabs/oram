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

// enumerator for different types of encryption level
enum class EncryptType {
  NONE,              // no encryption
  ENCRYPT,           // encrypt only, no authentication
  ENCRYPT_AND_AUTH,  // encrypt and authenticate with the index of the page, but
                     // no freshness check
  ENCRYPT_AND_AUTH_FRESH  // encrypt and authenticate with the index of the page
                          // and a freshness check. Note that the freshness
                          // check only works for 2^32 writes.
};

template <typename T, typename _BackendType = ::EM::Backend::MemServerBackend,
          const EncryptType enc_type = EncryptType::ENCRYPT_AND_AUTH,
          bool LATE_INIT = false>
struct NonCachedServerFrontendInstance {
  // Just forwards reads and writes to the server. Call the allocator during
  // construction and resizes.
  //

  using AllocatorSlot = typename EM::LargeBlockAllocator::AllocatorSlot;
  using BackendType = _BackendType;

  static constexpr bool ENCRYPTED = enc_type >= EncryptType::ENCRYPT;
  static constexpr bool AUTH = enc_type >= EncryptType::ENCRYPT_AND_AUTH;
  static constexpr bool FRESH_CHECK =
      enc_type >= EncryptType::ENCRYPT_AND_AUTH_FRESH;

  typedef uint64_t IndexType;
  BackendType& backend;

  static inline constexpr auto sizeOfT =
      ENCRYPTED ? sizeof(typename T::Encrypted_t) : sizeof(T);

  T defaultVal;

  std::vector<bool> modified;

  // counters for each page
  std::vector<uint32_t> counters;

  AllocatorSlot slot;

  typedef union {
    uint8_t bytes[IV_SIZE];
    struct {
      uint64_t index;
      uint32_t counter;  // in case an index is written multiple times
      uint32_t padding;  // pad to ensure counter part is right after index part
    };
  } nounce_t;

  nounce_t nounce;

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
    if constexpr (FRESH_CHECK) {
      std::swap(counters, other.counters);
    }
  }

  NonCachedServerFrontendInstance(BackendType& _backend, uint64_t initialSize,
                                  const T& _defaultVal)
      : backend(_backend) {
    if (initialSize == 0) {
      slot.base = -1;
      return;
    }
    IndexType requiredSize = initialSize * sizeOfT;
    slot = backend.Allocate(requiredSize);
    // std::cout << "Alloc: " << slot.base << "--" << slot.base + slot.size <<
    // std::endl;
    if constexpr (AUTH) {
      nounce.index = UniformRandom();
      nounce.counter = UniformRandom32();
    }
    if constexpr (FRESH_CHECK) {
      counters.resize(initialSize, 0);
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
      : NonCachedServerFrontendInstance(_backend, initialSize, T()) {}

  ~NonCachedServerFrontendInstance() {
    if (slot.base == (EM::LargeBlockAllocator::Size_t)-1) {
      return;
    }
    backend.Free(slot);
    slot.size = 0;
  }

  void Write(const IndexType i, const T& in) {
    if constexpr (LATE_INIT) {
      modified[i] = true;
    }
    PERFCTR_INCREMENT(writeCount);
    if constexpr (ENCRYPTED) {
      typename T::Encrypted_t inEnc;
      if constexpr (AUTH) {
        nounce_t nounceCopy = nounce;
        nounceCopy.index ^= i;
        if constexpr (FRESH_CHECK) {
          nounceCopy.counter ^= ++counters[i];
        }
        inEnc.Encrypt(in, nounceCopy.bytes);
      } else {
        inEnc.Encrypt(in);
      }

      backend.Write(slot.base + i * sizeOfT, sizeOfT,
                    reinterpret_cast<const uint8_t*>(&inEnc));
    } else {
      backend.Write(slot.base + i * sizeOfT, sizeOfT,
                    reinterpret_cast<const uint8_t*>(&in));
    }
  }

  // dummy implementation to be overloaded
  uint64_t WriteLazy(const IndexType i, const T& in) {
    Write(i, in);
    return 0;
  }

  void Read(const IndexType i, T& out) {
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
      if constexpr (AUTH) {
        nounce_t nounceCopy = nounce;
        nounceCopy.index ^= i;
        if constexpr (FRESH_CHECK) {
          nounceCopy.counter ^= counters[i];
        }

        inEnc.Decrypt(out, nounceCopy.bytes);
      } else {
        inEnc.Decrypt(out);
      }
    } else {
      backend.Read(slot.base + i * sizeOfT, sizeOfT,
                   reinterpret_cast<uint8_t*>(&out));
    }
  }

  uint64_t ReadLazy(const IndexType i, T& out) {
    Read(i, out);
    return 0;
  }

  void flushRead() {}

  void flushRead(uint32_t counter) {}

  void flushWrite() {}

  void flushWrite(uint32_t counter) {}
};

template <typename T, typename BackendType = ::EM::Backend::MemServerBackend,
          const EncryptType enc_type = EncryptType::ENCRYPT_AND_AUTH,
          uint64_t CACHE_SIZE = SERVER__CACHE_SIZE, uint64_t TLB_SIZE = 2>
using ServerFrontendInstance = CACHE::Cached<
    T, NonCachedServerFrontendInstance<T, BackendType, enc_type, false>,
    CACHE_SIZE, TLB_SIZE>;

}  // namespace MemoryServer
}  // namespace EM