#pragma once
#include <stdint.h>

#include "omap_interface.hpp"

struct OMapBindingSingleton {
  OMapGenericBinding<8, 8> omap;
  OMapBindingSingleton();
  void InitEmpty(uint32_t size);
  void InitEmptyExternal(uint32_t size, uint64_t cacheBytes);
  void StartInit(uint32_t size);
  void StartInitExternal(uint32_t size, uint64_t cacheBytes);
  void FinishInit();
  bool Insert(const void* keyPtr, const void* valPtr);
  bool OInsert(const void* keyPtr, const void* valPtr);
  bool Erase(const void* keyPtr);
  bool OErase(const void* keyPtr);
  bool Find(const void* keyPtr, void* valPtr);
  ~OMapBindingSingleton();
};