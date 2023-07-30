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

  ret = ecall_sort_perf(global_eid);
  // ret = ecall_test_sanity(global_eid);
  // ret = ecall_pageswap_perf(global_eid);
  // ret = ecall_pageswap_with_crypt_perf(global_eid);
  // ret = ecall_mergesplit_perf(global_eid);
  //   ret = ecall_mergesplit_compare(global_eid);

  if (ret != SGX_SUCCESS) abort();

  // printf("Did not abort\n");
}
