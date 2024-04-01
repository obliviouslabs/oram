#pragma once
#include <cinttypes>
#include <cstring>
#ifndef ENCLAVE_MODE
// #warning Defining enclave mode
#define ENCLAVE_MODE
// #define NDEBUG
#endif

#include "common/cpp_extended.hpp"

static uint64_t N;
uint8_t* data;

uint8_t* ocall_InitServer(uint64_t sizeOfT, uint64_t N_) {
  // printf("init mem server\n");
  data = new uint8_t[sizeOfT * N_];
  N = N_;
  return data;
}

void ocall_Read(size_t pos, uint64_t length, uint8_t* page) {
  std::memcpy(page, &data[pos], length);
}

void ocall_Write(uint64_t pos, uint64_t length, const uint8_t* page) {
  std::memcpy(&data[pos], page, length);
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

void ocall_DeleteServer() {
  if (data) {
    delete[] data;
  }
  data = NULL;
}