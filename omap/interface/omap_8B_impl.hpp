#pragma once
#include "omap_8B_interface.hpp"
#include "omap_impl.hpp"

OMapBindingSingleton::OMapBindingSingleton() {}

void OMapBindingSingleton::InitEmpty(uint32_t size) { omap.InitEmpty(size); }

void OMapBindingSingleton::InitEmptyExternal(uint32_t size,
                                             uint64_t cacheBytes) {
  omap.InitEmptyExternal(size, cacheBytes);
}

void OMapBindingSingleton::StartInit(uint32_t size) { omap.StartInit(size); }

void OMapBindingSingleton::StartInitExternal(uint32_t size,
                                             uint64_t cacheBytes) {
  omap.StartInitExternal(size, cacheBytes);
}

void OMapBindingSingleton::FinishInit() { omap.FinishInit(); }

bool OMapBindingSingleton::Insert(const void* keyPtr, const void* valPtr) {
  return omap.Insert(keyPtr, valPtr);
}

bool OMapBindingSingleton::OInsert(const void* keyPtr, const void* valPtr) {
  return omap.OInsert(keyPtr, valPtr);
}

bool OMapBindingSingleton::Erase(const void* keyPtr) {
  return omap.Erase(keyPtr);
}

bool OMapBindingSingleton::OErase(const void* keyPtr) {
  return omap.OErase(keyPtr);
}

bool OMapBindingSingleton::Find(const void* keyPtr, void* valPtr) {
  return omap.Find(keyPtr, valPtr);
}

OMapBindingSingleton::~OMapBindingSingleton() {}
