#ifndef OMAP_C_H
#define OMAP_C_H

#include <stdint.h>
// A pure C header file for C binding
#ifdef __cplusplus
extern "C" {
#endif

#define OMAP_TYPE(KS, VS) OMapType##_##KS##_##VS

enum OMapEnum {
#undef DECLARE_OMAP
#define DECLARE_OMAP(KS, VS) OMAP_TYPE(KS, VS) = KS * 65536 + VS,
#include "omap_declare.cfg"
#undef DECLARE_OMAP
};

struct OblivMap {
  void* omap;
  enum OMapEnum type;
};

typedef struct OblivMap OblivMap;

OblivMap NewOblivMap(int keySize, int valSize);
void InitEmpty(OblivMap omap, uint32_t size);
void InitEmptyExternal(OblivMap omap, uint32_t size, uint64_t cacheBytes);
void StartInit(OblivMap omap, uint32_t size);
void StartInitExternal(OblivMap omap, uint32_t size, uint64_t cacheBytes);
void FinishInit(OblivMap omap);
uint8_t Insert(OblivMap omap, const void* keyPtr, const void* valPtr);
uint8_t OInsert(OblivMap omap, const void* keyPtr, const void* valPtr);
uint8_t Erase(OblivMap omap, const void* keyPtr);
uint8_t OErase(OblivMap omap, const void* keyPtr);
uint8_t Find(OblivMap omap, const void* keyPtr, void* valPtr);
void Destroy(OblivMap omap);

#ifdef __cplusplus
}
#endif

#endif