#pragma once
#include "odsl/par_omap.hpp"
#include "par_omap_interface.hpp"

using ParMapType = ODSL::ParOMap<K, V, uint32_t>;
using ParInitializerType = typename ParMapType::InitContext;

ParOMapBinding::ParOMapBinding() {
  omap = nullptr;
  initializer = nullptr;
}

void ParOMapBinding::InitEmpty(uint32_t size, uint32_t numCores) {
  uint32_t shardCount = ParMapType::GetSuitableShardCount(numCores, true);
  omap = (void*)(new ParMapType(size, shardCount));
  ((ParMapType*)omap)->Init();
}

void ParOMapBinding::InitEmptyExternal(uint32_t size, uint32_t numCores,
                                       uint64_t cacheBytes) {
  uint32_t shardCount = ParMapType::GetSuitableShardCount(numCores, true);
  omap = (void*)(new ParMapType(size, shardCount));
  ((ParMapType*)omap)->Init(cacheBytes);
}

void ParOMapBinding::StartInit(uint32_t size, uint32_t initSize,
                               uint32_t numCores) {
  uint32_t shardCount = ParMapType::GetSuitableShardCount(numCores, false);
  omap = (void*)(new ParMapType(size, shardCount));
  initializer = (void*)(((ParMapType*)omap)->NewInitContext(initSize));
}

void ParOMapBinding::StartInitExternal(uint32_t size, uint32_t initSize,
                                       uint32_t numCores, uint64_t cacheBytes) {
  uint32_t shardCount = ParMapType::GetSuitableShardCount(numCores, false);
  omap = (void*)(new ParMapType(size, shardCount));
  initializer =
      (void*)(((ParMapType*)omap)->NewInitContext(initSize, cacheBytes));
}

void ParOMapBinding::FinishInit() {
  Assert(initializer, "FinishInit without StartInit");
  ((ParInitializerType*)initializer)->Finalize();
  delete (ParInitializerType*)initializer;
  initializer = nullptr;
}

void ParOMapBinding::InsertBatch(uint32_t batchSize, const K* keys,
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

void ParOMapBinding::FindBatch(uint32_t batchSize, const K* keys, V* vals,
                               bool* existFlags) {
  Assert(!initializer, "Find during initialization");
  std::vector<uint8_t> resFlags =
      ((ParMapType*)omap)->FindBatch(keys, keys + batchSize, vals);
  for (uint32_t i = 0; i < batchSize; ++i) {
    existFlags[i] = resFlags[i];
  }
}

void ParOMapBinding::FindBatchDeferMaintain(uint32_t batchSize, const K* keys,
                                            V* vals, bool* existFlags) {
  Assert(!initializer, "Find during initialization");
  std::vector<uint8_t> resFlags =
      ((ParMapType*)omap)
          ->FindBatchDeferWriteBack(keys, keys + batchSize, vals);
  for (uint32_t i = 0; i < batchSize; ++i) {
    existFlags[i] = resFlags[i];
  }
}

void ParOMapBinding::FindBatchMaintain() {
  Assert(!initializer, "Find during initialization");
  ((ParMapType*)omap)->WriteBack();
}

void ParOMapBinding::EraseBatch(uint32_t batchSize, const K* keys,
                                bool* existFlags) {
  Assert(!initializer, "Erase during initialization");
  std::vector<uint8_t> resFlags =
      ((ParMapType*)omap)->EraseBatch(keys, keys + batchSize);
  for (uint32_t i = 0; i < batchSize; ++i) {
    existFlags[i] = resFlags[i];
  }
}
