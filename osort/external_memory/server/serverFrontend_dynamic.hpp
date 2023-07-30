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

template <typename _BackendType = ::EM::Backend::MemServerBackend,
          bool ENCRYPTED = true, bool AUTH = true>
// #ifndef ENCLAVE_MODE
// requires EM::Backend::BackendServer<BackendType>
// #endif
struct DynamicServerFrontendInstance {
  // Just forwards reads and writes to the server. Call the allocator during
  // construction and resizes.
  //

  using AllocatorSlot = typename EM::LargeBlockAllocator::AllocatorSlot;
  using BackendType = _BackendType;

  typedef uint64_t IndexType;
  BackendType& backend;

  size_t sizeOfT;

  uint8_t* bufferEncrypted;

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

  DynamicServerFrontendInstance(DynamicServerFrontendInstance& other)
      : backend(other.backend), sizeOfT(other.sizeOfT) {
    slot = other.slot;
    if constexpr (AUTH) {
      nounce = other.nounce;
    }
    other.slot.base = -1;
    bufferEncrypted = (uint8_t*)malloc(sizeOfT + AES_BLOCK_SIZE);
    if (!bufferEncrypted) {
      printf("not enough space to launch frontend encrypted buffer page\n");
      abort();
    }
  }

  DynamicServerFrontendInstance(BackendType& _backend, uint64_t initialSize,
                                size_t _sizeOfT)
      : sizeOfT(_sizeOfT), backend(_backend) {
    IndexType requiredSize = initialSize * (sizeOfT + AES_BLOCK_SIZE);
    slot = backend.Allocate(requiredSize);

    // std::cout << "Alloc: " << slot.base << "--" << slot.base + slot.size <<
    // std::endl;
    if constexpr (AUTH) {
      nounce.identifiers.indexPart = UniformRandom();
      nounce.identifiers.counterPart = UniformRandom32();
    }

    bufferEncrypted = (uint8_t*)malloc(sizeOfT + AES_BLOCK_SIZE);
    if (!bufferEncrypted) {
      printf("not enough space to launch frontend encrypted buffer page\n");
      abort();
    }
  }

  ~DynamicServerFrontendInstance() {
    free(bufferEncrypted);
    if (slot.base == -1) {
      return;
    }
    // std::cout << "Freed: " << slot.base << "--" << slot.base + slot.size <<
    // std::endl;
    backend.Free(slot);
    slot.size = 0;
  }

  void Write(const IndexType i, const void* in) { Write(i, in, 0); }

  void Write(const IndexType i, const void* in, uint32_t counter) {
    // static_assert(AUTH);
    PERFCTR_INCREMENT(writeCount);

    nounce_t nounceCopy = nounce;
    nounceCopy.identifiers.indexPart ^= i;
    nounceCopy.identifiers.counterPart ^= counter;
    for (int r = 0; r < IO_ROUND; ++r) {
      aes_256_gcm_encrypt(sizeOfT, (uint8_t*)in, KEY, nounceCopy.bytes,
                          bufferEncrypted, bufferEncrypted + AES_BLOCK_SIZE);
    }
    backend.Write(slot.base + i * (sizeOfT + AES_BLOCK_SIZE),
                  sizeOfT + AES_BLOCK_SIZE, bufferEncrypted);
  }

  void Read(const IndexType i, void* out) { Read(i, out, 0); }

  void Read(const IndexType i, void* out, uint32_t counter) {
    PERFCTR_INCREMENT(readCount);
    nounce_t nounceCopy = nounce;
    nounceCopy.identifiers.indexPart ^= i;
    nounceCopy.identifiers.counterPart ^= counter;
    backend.Read(slot.base + i * (sizeOfT + AES_BLOCK_SIZE),
                 sizeOfT + AES_BLOCK_SIZE, bufferEncrypted);
    for (int r = 0; r < IO_ROUND; ++r) {
      bool success =
          aes_256_gcm_decrypt(sizeOfT, bufferEncrypted + AES_BLOCK_SIZE, KEY,
                              nounceCopy.bytes, bufferEncrypted, (uint8_t*)out);
      if (!success) {
        printf("authentication failure\n");
        abort();
      }
    }
  }
};

}  // namespace MemoryServer
}  // namespace EM