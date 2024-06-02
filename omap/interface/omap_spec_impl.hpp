#pragma once
#include "omap_generic_impl.hpp"
#include "omap_spec_interface.hpp"

#define DECLARE_OMAP_BINDING_IMPL(KS, VS)                                      \
  OMapBinding(KS, VS)::OMapBinding(KS, VS)() {}                                \
  void OMapBinding(KS, VS)::InitEmpty(uint32_t size) { omap.InitEmpty(size); } \
  void OMapBinding(KS, VS)::InitEmptyExternal(uint32_t size,                   \
                                              uint64_t cacheBytes) {           \
    omap.InitEmptyExternal(size, cacheBytes);                                  \
  }                                                                            \
  void OMapBinding(KS, VS)::StartInit(uint32_t size) { omap.StartInit(size); } \
  void OMapBinding(KS, VS)::StartInitExternal(uint32_t size,                   \
                                              uint64_t cacheBytes) {           \
    omap.StartInitExternal(size, cacheBytes);                                  \
  }                                                                            \
  void OMapBinding(KS, VS)::FinishInit() { omap.FinishInit(); }                \
  bool OMapBinding(KS, VS)::Insert(const void* keyPtr, const void* valPtr) {   \
    return omap.Insert(keyPtr, valPtr);                                        \
  }                                                                            \
  bool OMapBinding(KS, VS)::OInsert(const void* keyPtr, const void* valPtr) {  \
    return omap.OInsert(keyPtr, valPtr);                                       \
  }                                                                            \
  bool OMapBinding(KS, VS)::Erase(const void* keyPtr) {                        \
    return omap.Erase(keyPtr);                                                 \
  }                                                                            \
  bool OMapBinding(KS, VS)::OErase(const void* keyPtr) {                       \
    return omap.OErase(keyPtr);                                                \
  }                                                                            \
  bool OMapBinding(KS, VS)::Find(const void* keyPtr, void* valPtr) {           \
    return omap.Find(keyPtr, valPtr);                                          \
  }                                                                            \
  void OMapBinding(KS, VS)::Destroy() { omap.Destroy(); }                      \
  OMapBinding(KS, VS)::~OMapBinding(KS, VS)() {}

#undef DECLARE_OMAP
#define DECLARE_OMAP(KS, VS) DECLARE_OMAP_BINDING_IMPL(KS, VS)
#include "omap_declare.cfg"