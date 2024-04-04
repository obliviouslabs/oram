#pragma once
#include <stdint.h>

typedef void* omap_t;
typedef void* initializer_t;
typedef uint64_t K;
typedef uint64_t V;
/**
 * @brief Binds to OMap implementation.
 *
 */
struct OMapBindingSingleton {
  omap_t omap;                // pointer to the OMap object
  initializer_t initializer;  // pointer to the OMap initializer
  OMapBindingSingleton();

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
  bool Insert(K key, V val);

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
  bool OInsert(K key, V val);

  /**
   * @brief Finds a key in the map. The method is oblivious.
   *
   * @param key The key to find.
   * @param val Outputs the value corresponding to the key.
   * @return true if the key is found, false otherwise.
   */
  bool Find(K key, V& val);

  /**
   * @brief Erases a key from the map. The method is not oblivious, and can leak
   * if the key exists.
   *
   * @param key The key to erase.
   * @return true if the key existed in the map, false otherwise.
   */
  bool Erase(K key);

  /**
   * @brief Erases a key from the map obliviously. 1~2x slower than Erase.
   *
   * @param key The key to erase.
   * @return true if the key existed in the map, false otherwise.
   */
  bool OErase(K key);

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

  ~OMapBindingSingleton();
};
