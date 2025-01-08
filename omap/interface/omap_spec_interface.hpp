#pragma once
#include <stdint.h>

#include "omap_generic_interface.hpp"

#define OMapBinding(OMAP_KEY_SIZE, OMAP_VAL_SIZE) \
  OMapBinding##_##OMAP_KEY_SIZE##_##OMAP_VAL_SIZE

#define DECLARE_OMAP_BINDING(OMAP_KEY_SIZE, OMAP_VAL_SIZE)      \
  struct OMapBinding(OMAP_KEY_SIZE, OMAP_VAL_SIZE) {            \
    OMapGenericBinding<OMAP_KEY_SIZE, OMAP_VAL_SIZE> omap;      \
    OMapBinding(OMAP_KEY_SIZE, OMAP_VAL_SIZE)();                \
    void InitEmpty(uint32_t size);                              \
    void InitEmptyExternal(uint32_t size, uint64_t cacheBytes); \
    void StartInit(uint32_t size);                              \
    void StartInitExternal(uint32_t size, uint64_t cacheBytes); \
    void FinishInit();                                          \
    bool Insert(const void* keyPtr, const void* valPtr);        \
    bool OInsert(const void* keyPtr, const void* valPtr);       \
    bool Erase(const void* keyPtr);                             \
    bool OErase(const void* keyPtr);                            \
    bool Find(const void* keyPtr, void* valPtr);                \
    void Destroy();                                             \
    ~OMapBinding(OMAP_KEY_SIZE, OMAP_VAL_SIZE)();               \
  };

#define DECLARE_OMAP(KS, VS) DECLARE_OMAP_BINDING(KS, VS)
#include "omap_declare.cfg"
#undef DECLARE_OMAP