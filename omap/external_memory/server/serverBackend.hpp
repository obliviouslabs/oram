#pragma once
#include <cinttypes>
#include <concepts>
#include <cstring>
#include <vector>

#include "common/lrucache.hpp"
#include "common/tracing/tracer.hpp"
#include "external_memory/server/serverAllocator.hpp"

namespace EM {
enum PageSlotState { EMPTY_PAGE, PENDING_PAGE, DONE_PAGE };
namespace Backend {

template <class FS>
concept BackendServer = requires(
    FS fs, uint64_t indexType, const EM::LargeBlockAllocator::AllocatorSlot cas,
    uint8_t* u8p, const uint8_t* cu8p) {
  { FS(indexType) };
  {
    fs.Allocate(indexType)
  } -> std::same_as<const EM::LargeBlockAllocator::AllocatorSlot>;
  { fs.Free(cas) };
  { fs.Write(indexType, indexType, cu8p) };
  { fs.Read(indexType, indexType, u8p) };
};

#define SERVER_SIZE (1 << 16)

struct ServerBackend : EM::LargeBlockAllocator {
  // Read and Write bytes, contains an allocator.
  //
  using typename EM::LargeBlockAllocator::AllocatorSlot;
  using typename EM::LargeBlockAllocator::Size_t;

  explicit ServerBackend(uint64_t _size) : EM::LargeBlockAllocator(_size) {}
};

struct MemServerBackend : ServerBackend {
  uint8_t* data;
#ifndef ENCLAVE_MODE
  void ocall_InitServer(uint8_t** data_ptr, uint64_t pageSize,
                        uint64_t numPage) {
    *data_ptr = new uint8_t[pageSize * numPage];
  }

  void ocall_DeleteServer() { delete[] (data); }
#endif

#ifndef DISK_IO
  void ocall_Read(uint64_t offset, uint64_t sz, uint8_t* tmp) {
    std::memcpy(tmp, &data[offset], sz);
  }
  void ocall_Write(uint64_t offset, uint64_t sz, const uint8_t* tmp) {
    std::memcpy(&data[offset], tmp, sz);
  }
#endif

  explicit MemServerBackend(uint64_t _size) : ServerBackend(_size) {
    ocall_InitServer(&data, 4096, (size - 1) / 4096 + 1);
  }

  ~MemServerBackend() { ocall_DeleteServer(); }

  void Read(uint64_t offset, uint64_t sz, uint8_t* to) {
    ocall_Read(offset, sz, to);
  }

  void Write(uint64_t offset, uint64_t sz, const uint8_t* from) {
    ocall_Write(offset, sz, from);
  }

#ifndef DISK_IO
  void ocall_Read_Batch(uint64_t* offsets, uint64_t size, uint8_t* tmp,
                        uint64_t chunkNum, uint64_t totalSize) {
    uint8_t* pos = tmp;
    for (uint64_t i = 0; i < chunkNum; ++i) {
      uint64_t offset = *(offsets + i);
      ocall_Read(offset, size, pos);
      pos += size;
    }
  }
  void ocall_Write_Batch(uint64_t* offsets, uint64_t size, uint8_t* tmp,
                         uint64_t chunkNum, uint64_t totalSize) {
    uint8_t* pos = tmp;
    for (uint64_t i = 0; i < chunkNum; ++i) {
      uint64_t offset = *(offsets + i);
      ocall_Write(offset, size, pos);
      pos += size;
    }
  }
#endif
  void ReadBatch(uint64_t chunkNum, uint64_t pageSize, uint64_t* offsets,
                 uint8_t* out) {
    ocall_Read_Batch(offsets, pageSize, out, chunkNum, pageSize * chunkNum);
  }

  void WriteBatch(uint64_t chunkNum, uint64_t pageSize, uint64_t* offsets,
                  uint8_t* in) {
    ocall_Write_Batch(offsets, pageSize, in, chunkNum, pageSize * chunkNum);
  }
};
extern MemServerBackend* g_DefaultBackend;
}  // namespace Backend
}  // namespace EM