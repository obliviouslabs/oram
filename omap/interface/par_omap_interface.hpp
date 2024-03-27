#pragma once
#include <cinttypes>

typedef void* par_omap_t;
typedef void* initializer_t;

struct ParOMapBindingSingleton {
  par_omap_t omap;
  initializer_t initializer;
  ParOMapBindingSingleton();
  void InitEmpty(uint32_t size, uint32_t shardCount);
  void InitEmptyExternal(uint32_t size, uint32_t shardCount,
                         uint64_t cacheBytes);
  void InsertBatch(uint32_t batchSize, const uint64_t* keys,
                   const uint64_t* vals, bool* existFlags);
  void FindBatch(uint32_t batchSize, const uint64_t* keys, uint64_t* vals,
                 bool* existFlags);
  void FindBatchDeferMaintain(uint32_t batchSize, const uint64_t* keys,
                              uint64_t* vals, bool* existFlags);
  void FindBatchMaintain();

  void EraseBatch(uint32_t batchSize, const uint64_t* keys, bool* existFlags);
  void StartInit(uint32_t size, uint32_t initSize, uint32_t shardCount);
  void StartInitExternal(uint32_t size, uint32_t initSize, uint32_t shardCount,
                         uint64_t cacheBytes);
  void FinishInit();
  ~ParOMapBindingSingleton();
};
