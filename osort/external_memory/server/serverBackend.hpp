#pragma once
#include <cinttypes>
#include <cstring>
#include <vector>

#include "common/lrucache.hpp"
#include "common/tracing/tracer.hpp"
#include "external_memory/server/serverAllocator.hpp"

#ifndef ENCLAVE_MODE
#include <concepts>
#endif

#ifndef IO_ROUND
#define IO_ROUND 1
#endif

namespace EM {
enum PageSlotState { EMPTY_PAGE, PENDING_PAGE, DONE_PAGE };
namespace Backend {

#ifndef ENCLAVE_MODE
template <class FS>
concept BackendServer = requires(
    FS fs, uint64_t indexType, const EM::LargeBlockAllocator::AllocatorSlot cas,
    uint8_t* u8p, const uint8_t* cu8p) {
  { FS(indexType) };
  {
    fs.Allocate(indexType)
  } -> std::same_as<EM::LargeBlockAllocator::AllocatorSlot>;
  { fs.Free(cas) };
  { fs.Write(indexType, indexType, cu8p) };
  { fs.Read(indexType, indexType, u8p) };
};
#endif

#define SERVER_SIZE (1 << 16)

struct ServerBackend : EM::LargeBlockAllocator {
  // Read and Write bytes, contains an allocator.
  //
  using typename EM::LargeBlockAllocator::AllocatorSlot;
  using typename EM::LargeBlockAllocator::Size_t;

  ServerBackend(uint64_t _size) : EM::LargeBlockAllocator(_size) {}
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

  MemServerBackend(uint64_t _size) : ServerBackend(_size) {
    ocall_InitServer(&data, 4096, (size - 1) / 4096 + 1);
  }

  ~MemServerBackend() { ocall_DeleteServer(); }

  void Read(uint64_t offset, uint64_t sz, uint8_t* to) {
    for (int r = 0; r < IO_ROUND; ++r) {
      ocall_Read(offset, sz, to);
    }
  }

  void Write(uint64_t offset, uint64_t sz, const uint8_t* from) {
    for (int r = 0; r < IO_ROUND; ++r) {
      ocall_Write(offset, sz, from);
    }
  }

  const static uint64_t maxChunkNum = 16;
  struct Chunks {
    uint8_t tmp[4096 * maxChunkNum];
    uint8_t* addrs[maxChunkNum];
    uint64_t offsets[maxChunkNum];
    uint64_t sizes[maxChunkNum];
    PageSlotState* states[maxChunkNum];
    uint64_t chunkNum = 0;
    uint64_t tmpOffset = 0;
  };

  Chunks readChunks, writeChunks;

#ifndef DISK_IO
  void ocall_Read_Batch(uint64_t* offsets, uint64_t* sizes, uint8_t* tmp,
                        uint64_t chunkNum, uint64_t totalSize) {
    uint8_t* pos = tmp;
    for (uint64_t i = 0; i < chunkNum; ++i) {
      uint64_t offset = *(offsets + i);
      uint64_t size = *(sizes + i);
      ocall_Read(offset, size, pos);
      pos += size;
    }
  }
  void ocall_Write_Batch(uint64_t* offsets, uint64_t* sizes, uint8_t* tmp,
                         uint64_t chunkNum, uint64_t totalSize) {
    uint8_t* pos = tmp;
    for (uint64_t i = 0; i < chunkNum; ++i) {
      uint64_t offset = *(offsets + i);
      uint64_t size = *(sizes + i);
      ocall_Write(offset, size, pos);
      pos += size;
    }
  }
#endif

  void ReadLazy(uint64_t offset, uint64_t sz, uint8_t* to,
                PageSlotState& state) {
    Assert(state == PENDING_PAGE);
    Assert(sz <= 4096);
    uint64_t& chunkNum = readChunks.chunkNum;
    uint64_t& tmpOffset = readChunks.tmpOffset;
    if (tmpOffset + sz > sizeof(readChunks.tmp) || chunkNum + 1 > maxChunkNum) {
      FlushRead();
      Assert(chunkNum == 0);
      Assert(tmpOffset == 0);
    }

    readChunks.states[chunkNum] = &state;
    readChunks.addrs[chunkNum] = to;
    readChunks.offsets[chunkNum] = offset;
    readChunks.sizes[chunkNum] = sz;
    ++chunkNum;
    tmpOffset += sz;
  }

  void FlushRead() {
    uint64_t& chunkNum = readChunks.chunkNum;
    uint64_t& tmpOffset = readChunks.tmpOffset;
    uint8_t* tmp = readChunks.tmp;
    uint64_t* offsets = readChunks.offsets;
    uint64_t* sizes = readChunks.sizes;
    uint8_t** addrs = readChunks.addrs;
    if (chunkNum == 0) {
      Assert(tmpOffset == 0);
      return;
    }
    uint64_t totalSize = 0;
    for (size_t i = 0; i < chunkNum; ++i) {
      totalSize += sizes[i];
    }
    for (int r = 0; r < IO_ROUND; ++r) {
      ocall_Read_Batch(offsets, sizes, tmp, chunkNum, totalSize);
    }
    uint8_t* pos = tmp;
    for (uint64_t i = 0; i < chunkNum; ++i) {
      for (int r = 0; r < IO_ROUND; ++r) {
        std::memcpy(addrs[i], pos, sizes[i]);
      }
      pos += sizes[i];
      *readChunks.states[i] = DONE_PAGE;
    }
    chunkNum = 0;
    tmpOffset = 0;
  }

  void WriteLazy(uint64_t offset, uint64_t sz, const uint8_t* from,
                 PageSlotState& state) {
    Assert(state == PENDING_PAGE);
    Assert(sz <= 4096);
    uint64_t& chunkNum = writeChunks.chunkNum;
    uint64_t& tmpOffset = writeChunks.tmpOffset;
    if (tmpOffset + sz > sizeof(writeChunks.tmp) ||
        chunkNum + 1 > maxChunkNum) {
      FlushWrite();
      Assert(chunkNum == 0);
      Assert(tmpOffset == 0);
    }

    writeChunks.states[chunkNum] = &state;
    writeChunks.offsets[chunkNum] = offset;
    writeChunks.sizes[chunkNum] = sz;
    for (int r = 0; r < IO_ROUND; ++r) {
      std::memcpy(writeChunks.tmp + tmpOffset, from, sz);
    }
    ++chunkNum;
    tmpOffset += sz;
  }

  void FlushWrite() {
    uint64_t& chunkNum = writeChunks.chunkNum;
    uint64_t& tmpOffset = writeChunks.tmpOffset;
    uint8_t* tmp = writeChunks.tmp;
    uint64_t* offsets = writeChunks.offsets;
    uint64_t* sizes = writeChunks.sizes;
    if (chunkNum == 0) {
      Assert(tmpOffset == 0);
      return;
    }
    uint64_t totalSize = 0;
    for (size_t i = 0; i < chunkNum; ++i) {
      totalSize += sizes[i];
    }
    for (int r = 0; r < IO_ROUND; ++r) {
      ocall_Write_Batch(offsets, sizes, tmp, chunkNum, totalSize);
    }
    for (uint64_t i = 0; i < chunkNum; ++i) {
      *writeChunks.states[i] = DONE_PAGE;
    }
    chunkNum = 0;
    tmpOffset = 0;
  }
};
extern MemServerBackend* g_DefaultBackend;
}  // namespace Backend
}  // namespace EM