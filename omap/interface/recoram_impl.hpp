#pragma once
#include "recoram_interface.hpp"
#include "odsl/recursive_oram.hpp"

ORAMBindingSingleton::ORAMBindingSingleton() { oram = nullptr; }

using ORAMType = ODSL::RecursiveORAM<T, uint32_t>;

void ORAMBindingSingleton::InitORAM(uint32_t size) {
  Assert(oram == nullptr);
  oram = (void*)(new ORAMType(size));
  ((ORAMType*)oram)->InitDefault(T());
}

void ORAMBindingSingleton::InitORAMExternal(uint32_t size,
                                            uint64_t cacheBytes) {
  Assert(oram == nullptr);
  oram = (void*)(new ORAMType(size, cacheBytes));
  ((ORAMType*)oram)->InitDefault(T());
}

void ORAMBindingSingleton::Write(uint32_t addr, T val) {
  ((ORAMType*)oram)->Write(addr, val);
}

T ORAMBindingSingleton::Read(uint32_t addr) {
  T ret;
  ((ORAMType*)oram)->Read(addr, ret);
  return ret;
}

ORAMBindingSingleton::~ORAMBindingSingleton() {
  if (oram) {
    delete (ORAMType*)oram;
  }
}