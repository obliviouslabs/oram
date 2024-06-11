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
  initializer = nullptr;
  ((MapType*)omap)->Init();
}

template <uint64_t key_size, uint64_t val_size>
void OMapGenericBinding<key_size, val_size>::InitEmptyExternal(
    uint32_t size, uint64_t cacheBytes) {
  omap = (void*)(new MapType(size, cacheBytes));
  initializer = nullptr;
  ((MapType*)omap)->Init();
}

template <uint64_t key_size, uint64_t val_size>
void OMapGenericBinding<key_size, val_size>::StartInit(uint32_t size) {
  try {
    omap = (void*)(new MapType(size));
    initializer = (void*)(((MapType*)omap)->NewInitContext());
  } catch (const std::exception& e) {
    printf("Caught exception in StartInit: %s\n", e.what());
  }
}

template <uint64_t key_size, uint64_t val_size>
void OMapGenericBinding<key_size, val_size>::StartInitExternal(
    uint32_t size, uint64_t cacheBytes) {
  try {
    omap = (void*)(new MapType(size, cacheBytes));
    initializer = (void*)(((MapType*)omap)->NewInitContext());
  } catch (const std::exception& e) {
    printf("Caught exception in StartInitExternal: %s\n", e.what());
  }
}

template <uint64_t key_size, uint64_t val_size>
void OMapGenericBinding<key_size, val_size>::FinishInit() {
  Assert(initializer, "FinishInit without StartInit");
  try {
    ((InitializerType*)initializer)->Finalize();
    delete (InitializerType*)initializer;
    initializer = nullptr;
  } catch (const std::exception& e) {
    printf("Caught exception in FinishInit: %s\n", e.what());
  }
}

template <uint64_t key_size, uint64_t val_size>
bool OMapGenericBinding<key_size, val_size>::Insert(const void* keyPtr,
                                                    const void* valPtr) {
  Bytes<key_size> key(keyPtr);
  Bytes<val_size> val(valPtr);
  if (initializer) {
    try {
      ((InitializerType*)initializer)->Insert(key, val);
    } catch (const std::exception& e) {
      printf("Caught exception in Insert (during initialization): %s\n",
             e.what());
    }
    return false;
  }
  if (!omap) {
    printf("Insert: omap is null\n");
    return false;
  }
  try {
    return ((MapType*)omap)->Insert(key, val);
  } catch (const std::exception& e) {
    printf("Caught exception in Insert: %s\n", e.what());
  }
  return false;
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
  if (!omap) {
    printf("OInsert: omap is null\n");
    return false;
  }
  try {
    return ((MapType*)omap)->OInsert(key, val);
  } catch (const std::exception& e) {
    printf("Caught exception in OInsert: %s\n", e.what());
  }
  return false;
}

template <uint64_t key_size, uint64_t val_size>
bool OMapGenericBinding<key_size, val_size>::Find(const void* keyPtr,
                                                  void* valPtr) {
  Bytes<key_size> key(keyPtr);
  Bytes<val_size> val(valPtr);
  Assert(!initializer, "Find during initialization");
  if (!omap) {
    printf("Find: omap is null\n");
    return false;
  }
  try {
    bool found = ((MapType*)omap)->Find(key, val);
    memcpy(valPtr, val.GetData(), val_size);
    return found;
  } catch (const std::exception& e) {
    printf("Caught exception in Find: %s\n", e.what());
  }
  return false;
}

template <uint64_t key_size, uint64_t val_size>
bool OMapGenericBinding<key_size, val_size>::Erase(const void* keyPtr) {
  Assert(!initializer, "Erase during initialization");
  Bytes<key_size> key(keyPtr);
  if (!omap) {
    printf("Erase: omap is null\n");
    return false;
  }
  try {
    return ((MapType*)omap)->Erase(key);
  } catch (const std::exception& e) {
    printf("Caught exception in Erase: %s\n", e.what());
    return false;
  }
}

template <uint64_t key_size, uint64_t val_size>
bool OMapGenericBinding<key_size, val_size>::OErase(const void* keyPtr) {
  Assert(!initializer, "Erase during initialization");
  if (!omap) {
    printf("OErase: omap is null\n");
    return false;
  }
  try {
    Bytes<key_size> key(keyPtr);
    return ((MapType*)omap)->OErase(key);
  } catch (const std::exception& e) {
    printf("Caught exception in OErase: %s\n", e.what());
    return false;
  }
}

template <uint64_t key_size, uint64_t val_size>
void OMapGenericBinding<key_size, val_size>::Destroy() {
  if (omap) {
    delete (MapType*)omap;
  }
  if (initializer) {
    delete (InitializerType*)initializer;
  }
}

template <uint64_t key_size, uint64_t val_size>
OMapGenericBinding<key_size, val_size>::~OMapGenericBinding() {
  Destroy();
}
