#include "omap_c.h"

#include <stdio.h>

#include "omap_spec_interface.hpp"

static OMapEnum getOMapType(int keySize, int valSize) {
  OMapEnum omapEnum = static_cast<OMapEnum>(keySize * 65536 + valSize);
  return omapEnum;
}

extern "C" OblivMap NewOblivMap(int keySize, int valSize) {
  OblivMap omap;
  omap.type = getOMapType(keySize, valSize);
  switch (omap.type) {
#define DECLARE_OMAP(KS, VS)             \
  case OMAP_TYPE(KS, VS):                \
    omap.omap = new OMapBinding(KS, VS); \
    break;
#include "omap_declare.cfg"
#undef DECLARE_OMAP
    default:
      printf("Invalid omap type\n");
      omap.omap = nullptr;
      break;
  }
  return omap;
}
extern "C" void InitEmpty(OblivMap omap, uint32_t size) {
  switch (omap.type) {
#define DECLARE_OMAP(KS, VS)                                       \
  case OMAP_TYPE(KS, VS):                                          \
    static_cast<OMapBinding(KS, VS)*>(omap.omap)->InitEmpty(size); \
    break;

#include "omap_declare.cfg"
#undef DECLARE_OMAP
    default:
      printf("Invalid omap type\n");
      break;
  }
}
extern "C" void InitEmptyExternal(OblivMap omap, uint32_t size,
                                  uint64_t cacheBytes) {
  // omap->InitEmptyExternal(size, cacheBytes);
  switch (omap.type) {
#define DECLARE_OMAP(KS, VS)                                         \
  case OMAP_TYPE(KS, VS):                                            \
    static_cast<OMapBinding(KS, VS)*>(omap.omap)->InitEmptyExternal( \
        size, cacheBytes);                                           \
    break;
#include "omap_declare.cfg"
#undef DECLARE_OMAP
    default:
      printf("Invalid omap type\n");
      break;
  }
}
extern "C" void StartInit(OblivMap omap, uint32_t size) {
  switch (omap.type) {
#define DECLARE_OMAP(KS, VS)                                       \
  case OMAP_TYPE(KS, VS):                                          \
    static_cast<OMapBinding(KS, VS)*>(omap.omap)->StartInit(size); \
    break;
#include "omap_declare.cfg"
#undef DECLARE_OMAP
    default:
      printf("Invalid omap type\n");
      break;
  }
}
extern "C" void StartInitExternal(OblivMap omap, uint32_t size,
                                  uint64_t cacheBytes) {
  switch (omap.type) {
#define DECLARE_OMAP(KS, VS)                                         \
  case OMAP_TYPE(KS, VS):                                            \
    static_cast<OMapBinding(KS, VS)*>(omap.omap)->StartInitExternal( \
        size, cacheBytes);                                           \
    break;

#include "omap_declare.cfg"
#undef DECLARE_OMAP
    default:
      printf("Invalid omap type\n");
      break;
  }
  // omap->StartInitExternal(size, cacheBytes);
}
extern "C" void FinishInit(OblivMap omap) {
  switch (omap.type) {
#define DECLARE_OMAP(KS, VS)                                    \
  case OMAP_TYPE(KS, VS):                                       \
    static_cast<OMapBinding(KS, VS)*>(omap.omap)->FinishInit(); \
    break;

#include "omap_declare.cfg"
#undef DECLARE_OMAP
    default:
      printf("Invalid omap type\n");
      break;
  }
}

extern "C" uint8_t Insert(OblivMap omap, const void* keyPtr,
                          const void* valPtr) {
  switch (omap.type) {
#define DECLARE_OMAP(KS, VS)                                             \
  case OMAP_TYPE(KS, VS):                                                \
    return static_cast<OMapBinding(KS, VS)*>(omap.omap)->Insert(keyPtr,  \
                                                                valPtr); \
    break;
#include "omap_declare.cfg"
#undef DECLARE_OMAP
    default:
      printf("Invalid omap type\n");
      break;
  }
  return 0;
}
extern "C" uint8_t OInsert(OblivMap omap, const void* keyPtr,
                           const void* valPtr) {
  switch (omap.type) {
#define DECLARE_OMAP(KS, VS)                                              \
  case OMAP_TYPE(KS, VS):                                                 \
    return static_cast<OMapBinding(KS, VS)*>(omap.omap)->OInsert(keyPtr,  \
                                                                 valPtr); \
    break;
#include "omap_declare.cfg"
#undef DECLARE_OMAP
    default:
      printf("Invalid omap type\n");
      break;
  }
  return 0;
}
extern "C" uint8_t Erase(OblivMap omap, const void* keyPtr) {
  switch (omap.type) {
#define DECLARE_OMAP(KS, VS)                                            \
  case OMAP_TYPE(KS, VS):                                               \
    return static_cast<OMapBinding(KS, VS)*>(omap.omap)->Erase(keyPtr); \
    break;
#include "omap_declare.cfg"
#undef DECLARE_OMAP
    default:
      printf("Invalid omap type\n");
      break;
  }
  return 0;
}
extern "C" uint8_t OErase(OblivMap omap, const void* keyPtr) {
  switch (omap.type) {
#define DECLARE_OMAP(KS, VS)                                             \
  case OMAP_TYPE(KS, VS):                                                \
    return static_cast<OMapBinding(KS, VS)*>(omap.omap)->OErase(keyPtr); \
    break;
#include "omap_declare.cfg"
#undef DECLARE_OMAP
    default:
      printf("Invalid omap type\n");
      break;
  }
  return 0;
}
extern "C" uint8_t Find(OblivMap omap, const void* keyPtr, void* valPtr) {
  switch (omap.type) {
#define DECLARE_OMAP(KS, VS)                                                   \
  case OMAP_TYPE(KS, VS):                                                      \
    return static_cast<OMapBinding(KS, VS)*>(omap.omap)->Find(keyPtr, valPtr); \
    break;
#include "omap_declare.cfg"
#undef DECLARE_OMAP
    default:
      printf("Invalid omap type\n");
      break;
  }
  return 0;
}
extern "C" void Destroy(OblivMap omap) {
  switch (omap.type) {
#define DECLARE_OMAP(KS, VS)                                 \
  case OMAP_TYPE(KS, VS):                                    \
    static_cast<OMapBinding(KS, VS)*>(omap.omap)->Destroy(); \
    break;
#include "omap_declare.cfg"
#undef DECLARE_OMAP
    default:
      printf("Invalid omap type\n");
      break;
  }
}