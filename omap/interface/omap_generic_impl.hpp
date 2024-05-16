#pragma once
#include "odsl/omap.hpp"
#include "omap_generic_interface.hpp"

template <uint64_t key_size, uint64_t val_size>
using MapGenericType = ODSL::OMap<Bytes<key_size>, Bytes<val_size>, uint32_t>;

template <uint64_t key_size, uint64_t val_size>
using InitializerGenericType =
    typename MapGenericType<key_size, val_size>::InitContext;

#define MapType MapGenericType<key_size, val_size>
#define InitializerType InitializerGenericType<key_size, val_size>

template <uint64_t key_size, uint64_t val_size>
OMapGenericBinding<key_size, val_size>::OMapGenericBinding() {
  omap = nullptr;
  initializer = nullptr;
}

template <uint64_t key_size, uint64_t val_size>
void OMapGenericBinding<key_size, val_size>::InitEmpty(uint32_t size) {
  omap = (void*)(new MapType(size));
  ((MapType*)omap)->Init();
}

template <uint64_t key_size, uint64_t val_size>
void OMapGenericBinding<key_size, val_size>::InitEmptyExternal(
    uint32_t size, uint64_t cacheBytes) {
  omap = (void*)(new MapType(size, cacheBytes));
  ((MapType*)omap)->Init();
}

template <uint64_t key_size, uint64_t val_size>
void OMapGenericBinding<key_size, val_size>::StartInit(uint32_t size) {
  omap = (void*)(new MapType(size));
  initializer = (void*)(((MapType*)omap)->NewInitContext());
}

template <uint64_t key_size, uint64_t val_size>
void OMapGenericBinding<key_size, val_size>::StartInitExternal(
    uint32_t size, uint64_t cacheBytes) {
  omap = (void*)(new MapType(size, cacheBytes));
  initializer = (void*)(((MapType*)omap)->NewInitContext());
}

template <uint64_t key_size, uint64_t val_size>
void OMapGenericBinding<key_size, val_size>::FinishInit() {
  Assert(initializer, "FinishInit without StartInit");
  ((InitializerType*)initializer)->Finalize();
  delete (InitializerType*)initializer;
  initializer = nullptr;
}

template <uint64_t key_size, uint64_t val_size>
bool OMapGenericBinding<key_size, val_size>::Insert(const void* keyPtr,
                                                    const void* valPtr) {
  Bytes<key_size> key(keyPtr);
  Bytes<val_size> val(valPtr);
  if (initializer) {
    ((InitializerType*)initializer)->Insert(key, val);
    return false;
  }
  return ((MapType*)omap)->Insert(key, val);
}

template <uint64_t key_size, uint64_t val_size>
bool OMapGenericBinding<key_size, val_size>::OInsert(const void* keyPtr,
                                                     const void* valPtr) {
  Bytes<key_size> key(keyPtr);
  Bytes<val_size> val(valPtr);
  if (initializer) {
    ((InitializerType*)initializer)->Insert(key, val);
    return false;
  }
  return ((MapType*)omap)->OInsert(key, val);
}

template <uint64_t key_size, uint64_t val_size>
bool OMapGenericBinding<key_size, val_size>::Find(const void* keyPtr,
                                                  void* valPtr) {
  Bytes<key_size> key(keyPtr);
  Bytes<val_size> val(valPtr);
  Assert(!initializer, "Find during initialization");
  bool found = ((MapType*)omap)->Find(key, val);
  memcpy(valPtr, val.GetData(), val_size);
  return found;
}

template <uint64_t key_size, uint64_t val_size>
bool OMapGenericBinding<key_size, val_size>::Erase(const void* keyPtr) {
  Assert(!initializer, "Erase during initialization");
  Bytes<key_size> key(keyPtr);
  return ((MapType*)omap)->Erase(key);
}

template <uint64_t key_size, uint64_t val_size>
bool OMapGenericBinding<key_size, val_size>::OErase(const void* keyPtr) {
  Assert(!initializer, "Erase during initialization");
  Bytes<key_size> key(keyPtr);
  return ((MapType*)omap)->OErase(key);
}

template <uint64_t key_size, uint64_t val_size>
OMapGenericBinding<key_size, val_size>::~OMapGenericBinding() {
  if (omap) {
    delete (MapType*)omap;
  }
  if (initializer) {
    delete (InitializerType*)initializer;
  }
}
