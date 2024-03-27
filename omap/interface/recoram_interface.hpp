#pragma once
#include <cinttypes>

typedef void* oram_t;

struct ORAMBindingSingleton {
  oram_t oram;
  ORAMBindingSingleton();
  void InitORAM(uint32_t size);
  void Write(uint32_t addr, uint64_t val);
  uint64_t Read(uint32_t addr);
  ~ORAMBindingSingleton();
};
