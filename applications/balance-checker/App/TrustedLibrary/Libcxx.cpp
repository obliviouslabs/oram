#include <stdio.h>

#include <thread>

#include "../../common.hpp"
#include "../App.h"
#include "Enclave_u.h"
#ifdef DISK_IO
#include "external_memory/server/enclaveFileServer_untrusted.hpp"
#else
#include "external_memory/server/enclaveMemServer_untrusted.hpp"
#endif

void ActualMain(void) {
  sgx_status_t ret = SGX_ERROR_UNEXPECTED;
  size_t mapSize = 1000;
  ret = ecall_omap_init(global_eid, mapSize);

  if (ret != SGX_SUCCESS) abort();

  key_type key = 20;
  val_type val = 0;
  int findFlag = 0;
  ret = ecall_omap_find(global_eid, &findFlag, (uint8_t*)&key, (uint8_t*)&val,
                        sizeof(key_type), sizeof(val_type));
  if (ret != SGX_SUCCESS) abort();
  if (findFlag == 0) {
    printf("Did not find %lu\n", key);
    abort();
  }
  printf("find %lu %lu\n", key, val);
  // printf("Did not abort\n");
}
