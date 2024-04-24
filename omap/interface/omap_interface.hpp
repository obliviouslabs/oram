#pragma once
#include <stdint.h>

typedef void* omap_t;
typedef void* initializer_t;

/**
 * @brief Binds to OMap implementation.
 *
 */
template <uint64_t key_size, uint64_t val_size>
struct OMapGenericBinding {
  omap_t omap;                // pointer to the OMap object
  initializer_t initializer;  // pointer to the OMap initializer
  OMapGenericBinding();

  /**
   * @brief Initializes an empty OMap. Assumes sufficiently large enclave size.
   *
   * @param size The maximum number of key-value pairs the map can hold.
   */
  void InitEmpty(uint32_t size);

  /**
   * @brief Initializes an empty OMap. Uses at most cacheBytes bytes of EPC
   * memory.
   *
   * @param size The maximum number of key-value pairs the map can hold.
   * @param cacheBytes The maximum number of bytes to use for caching.
   */
  void InitEmptyExternal(uint32_t size, uint64_t cacheBytes);

  /**
   * @brief Inserts a key-value pair into the map. The method can be called
   * either during initialization or after initialization. The method leaks
   * whether the key already exists.
   *
   * @param key The key to insert.
   * @param val The value to insert.
   * @return Always false if called during initialization. Otherwise, returns
   * true if the key already exists, false otherwise.
   */
  bool Insert(const void* keyPtr, const void* valPtr);

  /**
   * @brief Inserts a key-value pair into the map obliviously. The method can be
   * called either during initialization or after initialization. If the method
   * is called during initialization, it is same as Insert. Otherwise, the
   * insertion is fully oblivious, and is about 3x slower than Insert.
   *
   * @param key
   * @param val
   * @return true
   * @return false
   */
  bool OInsert(const void* keyPtr, const void* valPtr);

  /**
   * @brief Finds a key in the map. The method is oblivious.
   *
   * @param key The key to find.
   * @param val Outputs the value corresponding to the key.
   * @return true if the key is found, false otherwise.
   */
  bool Find(const void* keyPtr, void* valPtr);

  /**
   * @brief Erases a key from the map. The method is not oblivious, and can leak
   * if the key exists.
   *
   * @param key The key to erase.
   * @return true if the key existed in the map, false otherwise.
   */
  bool Erase(const void* keyPtr);

  /**
   * @brief Erases a key from the map obliviously. 1~2x slower than Erase.
   *
   * @param key The key to erase.
   * @return true if the key existed in the map, false otherwise.
   */
  bool OErase(const void* keyPtr);

  /**
   * @brief Start initializing the OMap. Assumes sufficiently large EPC. Cannot
   * be called together with InitEmpty or InitEmptyExternal. Only Insert or
   * OInsert can be called during initialization.
   *
   * @param size The maximum number of key-value pairs the map can hold.
   */
  void StartInit(uint32_t size);

  /**
   * @brief Start initializing the OMap using at most cacheBytes bytes of EPC
   * memory. Cannot be called together with InitEmpty or InitEmptyExternal. Only
   * Insert or OInsert can be called during initialization.
   *
   * @param size The maximum number of key-value pairs the map can hold.
   */
  void StartInitExternal(uint32_t size, uint64_t cacheBytes);

  /**
   * @brief Finish the initialization. Must be called after StartInit or
   * StartInitExternal.
   *
   */
  void FinishInit();

  ~OMapGenericBinding();
};

#include "odsl/omap.hpp"

template <uint64_t key_size, uint64_t val_size>
using MapGenericType = ODSL::OHashMap<Bytes<key_size>, Bytes<val_size>, ODSL::FULL_OBLIVIOUS, uint32_t>;

template <uint64_t key_size, uint64_t val_size>
using InitializerGenericType = typename MapGenericType<key_size, val_size>::InitContext;

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
void OMapGenericBinding<key_size, val_size>::InitEmptyExternal(uint32_t size,
                                                     uint64_t cacheBytes) {
  omap = (void*)(new MapType(size, cacheBytes));
  ((MapType*)omap)->Init();
}

template <uint64_t key_size, uint64_t val_size>
void OMapGenericBinding<key_size, val_size>::StartInit(uint32_t size) {
  omap = (void*)(new MapType(size));
  initializer = (void*)(((MapType*)omap)->NewInitContext());
}

template <uint64_t key_size, uint64_t val_size>
void OMapGenericBinding<key_size, val_size>::StartInitExternal(uint32_t size,
                                                     uint64_t cacheBytes) {
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
bool OMapGenericBinding<key_size, val_size>::Insert(const void* keyPtr, const void* valPtr) {
  Bytes<key_size> key(keyPtr);
  Bytes<val_size> val(valPtr);
  if (initializer) {
    ((InitializerType*)initializer)->Insert(key, val);
    return false;
  }
  return ((MapType*)omap)->Insert(key, val);
}

template <uint64_t key_size, uint64_t val_size>
bool OMapGenericBinding<key_size, val_size>::OInsert(const void* keyPtr, const void* valPtr) {
  Bytes<key_size> key(keyPtr);
  Bytes<val_size> val(valPtr);
  if (initializer) {
    ((InitializerType*)initializer)->Insert(key, val);
    return false;
  }
  return ((MapType*)omap)->OInsert(key, val);
}

template <uint64_t key_size, uint64_t val_size>
bool OMapGenericBinding<key_size, val_size>::Find(const void* keyPtr, void* valPtr) {
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
