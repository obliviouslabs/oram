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
  requires EM::Backend::BackendServer<_BackendType>
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
  using Encrypted_t = typename T::Encrypted_t;
  static inline constexpr auto sizeOfT =
      ENCRYPTED ? sizeof(Encrypted_t) : sizeof(T);

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
      : backend(other.backend), defaultVal(other.defaultVal), slot(other.slot) {
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
    SetSize(initialSize, _defaultVal);
  }

  NonCachedServerFrontendInstance(BackendType& _backend, uint64_t initialSize)
      : NonCachedServerFrontendInstance(_backend) {
    SetSize(initialSize);
  }

  NonCachedServerFrontendInstance(BackendType& _backend)
      : backend(_backend), slot({-1UL, 0UL}) {
    if constexpr (AUTH) {
      nounce.index = 0;
      nounce.counter = 0;
    }
  }

  void SetSize(uint64_t initialSize, const T& _defaultVal) {
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

  void SetSize(uint64_t initialSize) { SetSize(initialSize, T()); }

  ~NonCachedServerFrontendInstance() {
    if (slot.base == (EM::LargeBlockAllocator::Size_t)-1) {
      return;
    }
    backend.Free(slot);
    slot.size = 0;
  }

  void Encrypt(const T& in, Encrypted_t& inEnc, const IndexType i) {
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
  }

  void Decrypt(Encrypted_t& inEnc, T& out, const IndexType i) {
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
  }

  void Write(const IndexType i, const T& in) {
    if constexpr (LATE_INIT) {
      modified[i] = true;
    }
    PERFCTR_INCREMENT(writeCount);
    if constexpr (ENCRYPTED) {
      Encrypted_t inEnc;
      Encrypt(in, inEnc, i);
      backend.Write(slot.base + i * sizeOfT, sizeOfT,
                    reinterpret_cast<const uint8_t*>(&inEnc));
    } else {
      backend.Write(slot.base + i * sizeOfT, sizeOfT,
                    reinterpret_cast<const uint8_t*>(&in));
    }
  }

  std::vector<Encrypted_t> writeBackBuffer;
  std::vector<uint64_t> writeBackBufferOffsets;

  // dummy implementation to be overloaded
  void WriteLazy(const IndexType i, const T& in) {
    // #ifdef DISK_IO
    writeBackBuffer.emplace_back();
    Encrypt(in, writeBackBuffer.back(), i);
    writeBackBufferOffsets.push_back(slot.base + i * sizeOfT);
    // #else
    //     Write(i, in);
    // #endif
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
      Encrypted_t inEnc;
      backend.Read(slot.base + i * sizeOfT, sizeOfT,
                   reinterpret_cast<uint8_t*>(&inEnc));
      Decrypt(inEnc, out, i);
    } else {
      backend.Read(slot.base + i * sizeOfT, sizeOfT,
                   reinterpret_cast<uint8_t*>(&out));
    }
  }

  std::vector<Encrypted_t> readBuffer;
  std::vector<uint64_t> readBufferOffsets;
  std::vector<T*> readBufferOuts;

  void ReadLazy(const IndexType i, T& out) {
    // #ifdef DISK_IO
    readBufferOffsets.push_back(slot.base + i * sizeOfT);
    readBufferOuts.push_back(&out);
    // #else
    //     Read(i, out);
    // #endif
  }

  void flushRead() {
    uint64_t batchSize = readBufferOffsets.size();
    if (batchSize == 0) {
      return;
    }
    readBuffer.resize(batchSize);
    backend.ReadBatch(readBufferOffsets.size(), sizeOfT,
                      readBufferOffsets.data(),
                      reinterpret_cast<uint8_t*>(readBuffer.data()));
    for (uint64_t i = 0; i < batchSize; i++) {
      Decrypt(readBuffer[i], *readBufferOuts[i],
              (IndexType)((readBufferOffsets[i] - slot.base) / sizeOfT));
    }
    readBufferOffsets.clear();
    readBufferOuts.clear();
    readBuffer.clear();
  }

  void flushWrite() {
    uint64_t batchSize = writeBackBufferOffsets.size();
    if (batchSize == 0) {
      return;
    }
    backend.WriteBatch(writeBackBufferOffsets.size(), sizeOfT,
                       writeBackBufferOffsets.data(),
                       reinterpret_cast<uint8_t*>(writeBackBuffer.data()));
    writeBackBufferOffsets.clear();
    writeBackBuffer.clear();
  }
};

template <typename T, typename BackendType = ::EM::Backend::MemServerBackend,
          const EncryptType enc_type = EncryptType::ENCRYPT_AND_AUTH,
          uint64_t cache_size = SERVER__CACHE_SIZE, uint64_t tlb_size = 2>
using ServerFrontendInstance = CACHE::Cached<
    T, NonCachedServerFrontendInstance<T, BackendType, enc_type, false>,
    cache_size, tlb_size>;

}  // namespace MemoryServer
}  // namespace EM