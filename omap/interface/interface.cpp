#include "interface/recoram_interface.hpp"
#include "odsl/recursive_oram.hpp"

using namespace ODSL;

ORAMBindingSingleton::ORAMBindingSingleton() { oram = nullptr; }

using ORAMType = RecursiveORAM<T, uint32_t>;

void ORAMBindingSingleton::InitORAM(uint32_t size) {
  Assert(oram == nullptr);
  oram = (void*)(new ORAMType(size));
  ((ORAMType*)oram)->InitDefault(T());
}

void ORAMBindingSingleton::InitORAMExternal(uint32_t size,
                                            uint64_t cacheBytes) {
  Assert(oram == nullptr);
  oram = (void*)(new ORAMType(size, cacheBytes));
  ((ORAMType*)oram)->InitDefault(T());
}

void ORAMBindingSingleton::Write(uint32_t addr, T val) {
  ((ORAMType*)oram)->Write(addr, val);
}

T ORAMBindingSingleton::Read(uint32_t addr) {
  T ret;
  ((ORAMType*)oram)->Read(addr, ret);
  return ret;
}

ORAMBindingSingleton::~ORAMBindingSingleton() {
  if (oram) {
    delete (ORAMType*)oram;
  }
}

#include "interface/omap_interface.hpp"
#include "odsl/omap.hpp"

using MapType = OHashMap<K, V, true, uint32_t>;
using InitializerType = typename MapType::InitContext;

OMapBindingSingleton::OMapBindingSingleton() {
  omap = nullptr;
  initializer = nullptr;
}

void OMapBindingSingleton::InitEmpty(uint32_t size) {
  omap = (void*)(new MapType(size));
  ((MapType*)omap)->Init();
}

void OMapBindingSingleton::InitEmptyExternal(uint32_t size,
                                             uint64_t cacheBytes) {
  omap = (void*)(new MapType(size, cacheBytes));
  ((MapType*)omap)->Init();
}

void OMapBindingSingleton::StartInit(uint32_t size) {
  omap = (void*)(new MapType(size));
  initializer = (void*)(((MapType*)omap)->NewInitContext());
}

void OMapBindingSingleton::StartInitExternal(uint32_t size,
                                             uint64_t cacheBytes) {
  omap = (void*)(new MapType(size, cacheBytes));
  initializer = (void*)(((MapType*)omap)->NewInitContext());
}

void OMapBindingSingleton::FinishInit() {
  Assert(initializer, "FinishInit without StartInit");
  ((InitializerType*)initializer)->Finalize();
  delete (InitializerType*)initializer;
  initializer = nullptr;
}

bool OMapBindingSingleton::Insert(K key, V val) {
  if (initializer) {
    ((InitializerType*)initializer)->Insert(key, val);
    return false;
  }
  return ((MapType*)omap)->Insert(key, val);
}

bool OMapBindingSingleton::OInsert(K key, V val) {
  if (initializer) {
    ((InitializerType*)initializer)->Insert(key, val);
    return false;
  }
  return ((MapType*)omap)->OInsert(key, val);
}

bool OMapBindingSingleton::Find(K key, V& val) {
  Assert(!initializer, "Find during initialization");
  return ((MapType*)omap)->Find(key, val);
}

bool OMapBindingSingleton::Erase(K key) {
  Assert(!initializer, "Erase during initialization");
  return ((MapType*)omap)->Erase(key);
}

bool OMapBindingSingleton::OErase(K key) {
  Assert(!initializer, "Erase during initialization");
  return ((MapType*)omap)->OErase(key);
}

OMapBindingSingleton::~OMapBindingSingleton() {
  if (omap) {
    delete (MapType*)omap;
  }
  if (initializer) {
    delete (InitializerType*)initializer;
  }
}

#include "interface/par_omap_interface.hpp"
#include "odsl/par_omap.hpp"

using ParMapType = ParOMap<K, V, uint32_t>;
using ParInitializerType = typename ParMapType::InitContext;

ParOMapBindingSingleton::ParOMapBindingSingleton() {
  omap = nullptr;
  initializer = nullptr;
}

void ParOMapBindingSingleton::InitEmpty(uint32_t size, uint32_t numCores) {
  uint32_t shardCount = ParMapType::GetSuitableShardCount(numCores, true);
  omap = (void*)(new ParMapType(size, shardCount));
  ((ParMapType*)omap)->Init();
}

void ParOMapBindingSingleton::InitEmptyExternal(uint32_t size,
                                                uint32_t numCores,
                                                uint64_t cacheBytes) {
  uint32_t shardCount = ParMapType::GetSuitableShardCount(numCores, true);
  omap = (void*)(new ParMapType(size, shardCount));
  ((ParMapType*)omap)->Init(cacheBytes);
}

void ParOMapBindingSingleton::StartInit(uint32_t size, uint32_t initSize,
                                        uint32_t numCores) {
  uint32_t shardCount = ParMapType::GetSuitableShardCount(numCores, false);
  omap = (void*)(new ParMapType(size, shardCount));
  initializer = (void*)(((ParMapType*)omap)->NewInitContext(initSize));
}

void ParOMapBindingSingleton::StartInitExternal(uint32_t size,
                                                uint32_t initSize,
                                                uint32_t numCores,
                                                uint64_t cacheBytes) {
  uint32_t shardCount = ParMapType::GetSuitableShardCount(numCores, false);
  omap = (void*)(new ParMapType(size, shardCount));
  initializer =
      (void*)(((ParMapType*)omap)->NewInitContext(initSize, cacheBytes));
}

void ParOMapBindingSingleton::FinishInit() {
  Assert(initializer, "FinishInit without StartInit");
  ((ParInitializerType*)initializer)->Finalize();
  delete (ParInitializerType*)initializer;
  initializer = nullptr;
}

void ParOMapBindingSingleton::InsertBatch(uint32_t batchSize, const K* keys,
                                          const V* vals, bool* existFlags) {
  if (initializer) {
    for (uint32_t i = 0; i < batchSize; ++i) {
      ((ParInitializerType*)initializer)->Insert(keys[i], vals[i]);
    }
    return;
  }
  std::vector<uint8_t> resFlags =
      ((ParMapType*)omap)->InsertBatch(keys, keys + batchSize, vals);
  for (uint32_t i = 0; i < batchSize; ++i) {
    existFlags[i] = resFlags[i];
  }
}

void ParOMapBindingSingleton::FindBatch(uint32_t batchSize, const K* keys,
                                        V* vals, bool* existFlags) {
  Assert(!initializer, "Find during initialization");
  std::vector<uint8_t> resFlags =
      ((ParMapType*)omap)->FindBatch(keys, keys + batchSize, vals);
  for (uint32_t i = 0; i < batchSize; ++i) {
    existFlags[i] = resFlags[i];
  }
}

void ParOMapBindingSingleton::FindBatchDeferMaintain(uint32_t batchSize,
                                                     const K* keys, V* vals,
                                                     bool* existFlags) {
  Assert(!initializer, "Find during initialization");
  std::vector<uint8_t> resFlags =
      ((ParMapType*)omap)
          ->FindBatchDeferWriteBack(keys, keys + batchSize, vals);
  for (uint32_t i = 0; i < batchSize; ++i) {
    existFlags[i] = resFlags[i];
  }
}

void ParOMapBindingSingleton::FindBatchMaintain() {
  Assert(!initializer, "Find during initialization");
  ((ParMapType*)omap)->WriteBack();
}

void ParOMapBindingSingleton::EraseBatch(uint32_t batchSize, const K* keys,
                                         bool* existFlags) {
  Assert(!initializer, "Erase during initialization");
  std::vector<uint8_t> resFlags =
      ((ParMapType*)omap)->EraseBatch(keys, keys + batchSize);
  for (uint32_t i = 0; i < batchSize; ++i) {
    existFlags[i] = resFlags[i];
  }
}

#include "interface/common_interface.hpp"
void ResetBackend(uint64_t size) {
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
  EM::Backend::g_DefaultBackend = new EM::Backend::MemServerBackend(size);
}

void DeleteBackend() {
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
}

void HelloWorld(uint32_t num) { printf("Hello, world %u!\n", num); }