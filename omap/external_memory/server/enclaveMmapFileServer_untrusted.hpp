#pragma once
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cinttypes>
#include <cstring>
#include <fstream>
#include <iostream>

#include "common/cpp_extended.hpp"

#define DEFAULT_FILENAME "/ssdmount/storage.bin"

static uint64_t N;
int fd;
uint8_t* data;
uint8_t* ocall_InitServer(uint64_t sizeOfT, uint64_t N_) {
  N = N_;
  int fd = open(DEFAULT_FILENAME, O_RDWR | O_CREAT | O_TRUNC, 0644);
  int allocRes = fallocate64(fd, 0, 0, N * sizeOfT);
  if (allocRes != 0) {
    std::cerr << "Error allocating space for the file." << std::endl;
    close(fd);
    return NULL;
  }
  if (fd == -1) {
    std::cerr << "Error opening file for writing." << std::endl;
    return NULL;
  }
  data = static_cast<uint8_t*>(
      mmap(nullptr, N * sizeOfT, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0));
  if (data == MAP_FAILED) {
    std::cerr << "Error memory mapping the file." << std::endl;
    close(fd);
    return NULL;
  }
  return NULL;
}

void ocall_Read(size_t pos, uint64_t length, uint8_t* page) {
  std::memcpy(page, data + pos, length);
}

void ocall_Write(uint64_t pos, uint64_t length, const uint8_t* page) {
  std::memcpy(data + pos, page, length);
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
  if (munmap(data, N) == -1) {
    std::cerr << "Error unmapping the file." << std::endl;
  }
  close(fd);
}