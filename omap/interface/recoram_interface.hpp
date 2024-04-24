#pragma once
#include <stdint.h>

typedef void* oram_t;
typedef uint64_t T;

struct ORAMBindingSingleton {
  oram_t oram;  // pointer to the ORAM object
  ORAMBindingSingleton();

  /**
   * @brief Initializes an empty ORAM. Assumes sufficiently large enclave size.
   *
   * @param size The maximum number of elements the ORAM can hold.
   */
  void InitORAM(uint32_t size);

  /**
   * @brief Initializes an empty ORAM using at most cacheBytes bytes of EPC
   * memory.
   *
   * @param size The maximum number of elements the ORAM can hold.
   * @param cacheBytes The maximum number of bytes to use for caching.
   */
  void InitORAMExternal(uint32_t size, uint64_t cacheBytes);

  /**
   * @brief Writes a value to the ORAM.
   *
   * @param addr The address to write to.
   * @param val The value to write.
   */
  void Write(uint32_t addr, T val);

  /**
   * @brief Reads a value from the ORAM.
   *
   * @param addr The address to read from.
   * @return The value at the address.
   */
  T Read(uint32_t addr);

  ~ORAMBindingSingleton();
};

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