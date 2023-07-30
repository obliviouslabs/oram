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

void ocall_DeleteServer() {
  if (data) {
    delete[] data;
  }
  data = NULL;
}