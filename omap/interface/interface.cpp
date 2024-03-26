#include "interface/recoram_interface.hpp"

#include "odsl/recursive_oram.hpp"

using namespace ODSL;

ORAMBindingSingleton::ORAMBindingSingleton() {
  uint64_t BackendSize = 1e10;
  EM::Backend::g_DefaultBackend =
      new EM::Backend::MemServerBackend(BackendSize);
  oram = nullptr;
}

void ORAMBindingSingleton::InitORAM(uint32_t size) {
  // ASSERT(oram == nullptr);
  oram = (void*)(new RecursiveORAM<uint64_t, uint32_t>(size));
  ((RecursiveORAM<uint64_t, uint32_t>*)oram)->InitDefault(0);
}

void ORAMBindingSingleton::Write(uint32_t addr, uint64_t val) {
  ((RecursiveORAM<uint64_t, uint32_t>*)oram)->Write(addr, val);
}

uint64_t ORAMBindingSingleton::Read(uint32_t addr) {
  uint64_t ret;
  ((RecursiveORAM<uint64_t, uint32_t>*)oram)->Read(addr, ret);
  return ret;
}

#include "interface/omap_interface.hpp"
#include "odsl/omap.hpp"

using namespace ODSL;

using K = uint64_t;
using V = uint64_t;
using MapType = OHashMap<K, V, true, uint32_t>;

OMapBindingSingleton::OMapBindingSingleton() {
  uint64_t BackendSize = 1e10;
  EM::Backend::g_DefaultBackend =
      new EM::Backend::MemServerBackend(BackendSize);
  omap = nullptr;
}

void OMapBindingSingleton::InitOMap(uint32_t size) {
  // ASSERT(omap == nullptr);
  omap = (void*)(new MapType(size));
  ((MapType*)omap)->Init();
}

bool OMapBindingSingleton::Insert(K key, V val) {
  return ((MapType*)omap)->Insert(key, val);
}

bool OMapBindingSingleton::OInsert(K key, V val) {
  return ((MapType*)omap)->OInsert(key, val);
}

bool OMapBindingSingleton::Find(K key, V& val) {
  return ((MapType*)omap)->Find(key, val);
}

bool OMapBindingSingleton::Erase(K key) { return ((MapType*)omap)->Erase(key); }

bool OMapBindingSingleton::OErase(K key) {
  return ((MapType*)omap)->OErase(key);
}
