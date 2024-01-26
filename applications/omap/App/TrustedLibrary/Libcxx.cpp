#include <stdio.h>

#include <thread>

#include "../App.h"
#include "Enclave_u.h"
#ifdef DISK_IO
#include "external_memory/server/enclaveFileServer_untrusted.hpp"
#else
#include "external_memory/server/enclaveMemServer_untrusted.hpp"
#endif

void ActualMain(void) {
  sgx_status_t ret = SGX_ERROR_UNEXPECTED;

  ret = ecall_omap_perf(global_eid);

  if (ret != SGX_SUCCESS) abort();

  // printf("Did not abort\n");
}
