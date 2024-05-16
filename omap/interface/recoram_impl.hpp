#pragma once
#include "odsl/recursive_oram.hpp"
#include "recoram_interface.hpp"

ORAMBinding::ORAMBinding() { oram = nullptr; }

using ORAMType = ODSL::RecursiveORAM<T, uint32_t>;

void ORAMBinding::InitORAM(uint32_t size) {
  Assert(oram == nullptr);
  oram = (void*)(new ORAMType(size));
  ((ORAMType*)oram)->InitDefault(T());
}

void ORAMBinding::InitORAMExternal(uint32_t size, uint64_t cacheBytes) {
  Assert(oram == nullptr);
  oram = (void*)(new ORAMType(size, cacheBytes));
  ((ORAMType*)oram)->InitDefault(T());
}

void ORAMBinding::Write(uint32_t addr, T val) {
  ((ORAMType*)oram)->Write(addr, val);
}

T ORAMBinding::Read(uint32_t addr) {
  T ret;
  ((ORAMType*)oram)->Read(addr, ret);
  return ret;
}

ORAMBinding::~ORAMBinding() {
  if (oram) {
    delete (ORAMType*)oram;
  }
}