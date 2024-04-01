#pragma once
#include <cinttypes>
#include <fstream>

#include "common/cpp_extended.hpp"

#define DEFAULT_FILENAME "/ssdmount/storage.bin"

static uint64_t N;
static std::unique_ptr<std::fstream> lbios(nullptr);

uint8_t* ocall_InitServer(uint64_t sizeOfT, uint64_t N_) {
  lbios.reset(new std::fstream());
  N = N_;
  lbios->open(DEFAULT_FILENAME, std::fstream::in | std::fstream::out |
                                    std::fstream::binary | std::fstream::trunc);
  Assert(lbios->is_open());
  return NULL;
}

void ocall_Read(size_t pos, uint64_t length, uint8_t* page) {
  std::ifstream file(DEFAULT_FILENAME, std::ios::binary);
  std::streampos filePos = pos;
  file.seekg(filePos, std::ios::beg);
  file.read((char*)page, length);
  if (!file) {
    throw std::runtime_error("read failed");
  }
}

void ocall_Write(uint64_t pos, uint64_t length, const uint8_t* page) {
  std::ofstream file(DEFAULT_FILENAME,
                     std::ios::binary | std::ios::in | std::ios::out);
  std::streampos filePos = pos;
  file.seekp(filePos, std::ios::beg);
  file.write((char*)page, length);
  if (!file) {
    throw std::runtime_error("write failed");
  }
}

void ocall_Read_Batch(uint64_t batchSize, uint64_t pageBytes,
                      uint64_t totalBytes, uint64_t* offsets, uint8_t* buffer) {
  // printf("batch read %ld chunks of total size %ld\n", batchSize, totalBytes);
  uint8_t* pos = buffer;
  for (uint64_t i = 0; i < batchSize; ++i) {
    uint64_t offset = *(offsets + i);
    ocall_Read(offset, pageBytes, pos);
    pos += pageBytes;
  }
}
void ocall_Write_Batch(uint64_t batchSize, uint64_t pageBytes,
                       uint64_t totalBytes, uint64_t* offsets,
                       uint8_t* buffer) {
  //  printf("batch write %ld chunks of total size %ld\n", batchSize,
  //  totalBytes);
  uint8_t* pos = buffer;
  for (uint64_t i = 0; i < batchSize; ++i) {
    uint64_t offset = *(offsets + i);
    ocall_Write(offset, pageBytes, pos);
    pos += pageBytes;
  }
}

void ocall_DeleteServer() { lbios->close(); }