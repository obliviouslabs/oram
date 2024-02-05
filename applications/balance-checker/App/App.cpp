#include "App.h"

#include <assert.h>
#include <openssl/sha.h>
#include <pwd.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "Enclave_u.h"
#include "sgx_dcap_ql_wrapper.h"
#include "sgx_eid.h"   /* sgx_enclave_id_t */
#include "sgx_error.h" /* sgx_status_t */
#include "sgx_quote_3.h"
#include "sgx_urts.h"

/* Global EID shared by multiple threads */
sgx_enclave_id_t global_eid = 0;

typedef struct _sgx_errlist_t {
  sgx_status_t err;
  const char *msg;
  const char *sug; /* Suggestion */
} sgx_errlist_t;

/* Error code returned by sgx_create_enclave */
static sgx_errlist_t sgx_errlist[] = {
    {SGX_ERROR_UNEXPECTED, "Unexpected error occurred.", NULL},
    {SGX_ERROR_INVALID_PARAMETER, "Invalid parameter.", NULL},
    {SGX_ERROR_OUT_OF_MEMORY, "Out of memory.", NULL},
    {SGX_ERROR_ENCLAVE_LOST, "Power transition occurred.",
     "Please refer to the sample \"PowerTransition\" for details."},
    {SGX_ERROR_INVALID_ENCLAVE, "Invalid enclave image.", NULL},
    {SGX_ERROR_INVALID_ENCLAVE_ID, "Invalid enclave identification.", NULL},
    {SGX_ERROR_INVALID_SIGNATURE, "Invalid enclave signature.", NULL},
    {SGX_ERROR_OUT_OF_EPC, "Out of EPC memory.", NULL},
    {SGX_ERROR_NO_DEVICE, "Invalid SGX device.",
     "Please make sure SGX module is enabled in the BIOS, and install SGX "
     "driver afterwards."},
    {SGX_ERROR_MEMORY_MAP_CONFLICT, "Memory map conflicted.", NULL},
    {SGX_ERROR_INVALID_METADATA, "Invalid enclave metadata.", NULL},
    {SGX_ERROR_DEVICE_BUSY, "SGX device was busy.", NULL},
    {SGX_ERROR_INVALID_VERSION, "Enclave version was invalid.", NULL},
    {SGX_ERROR_INVALID_ATTRIBUTE, "Enclave was not authorized.", NULL},
    {SGX_ERROR_ENCLAVE_FILE_ACCESS, "Can't open enclave file.", NULL},
    {SGX_ERROR_NDEBUG_ENCLAVE,
     "The enclave is signed as product enclave, and can not be created as "
     "debuggable enclave.",
     NULL},
};

/* Check error conditions for loading enclave */
void print_error_message(sgx_status_t ret) {
  size_t idx = 0;
  size_t ttl = sizeof sgx_errlist / sizeof sgx_errlist[0];

  for (idx = 0; idx < ttl; idx++) {
    if (ret == sgx_errlist[idx].err) {
      if (NULL != sgx_errlist[idx].sug)
        printf("Info: %s\n", sgx_errlist[idx].sug);
      printf("Error: %s\n", sgx_errlist[idx].msg);
      break;
    }
  }

  if (idx == ttl) printf("Error: Unexpected error occurred.\n");
}

/* Initialize the enclave:
 *   Call sgx_create_enclave to initialize an enclave instance
 */
int initialize_enclave(void) {
  sgx_status_t ret = SGX_ERROR_UNEXPECTED;

  /* Call sgx_create_enclave to initialize an enclave instance */
  /* Debug Support: set 2nd parameter to 1 */
  ret = sgx_create_enclave(ENCLAVE_FILENAME, SGX_DEBUG_FLAG, NULL, NULL,
                           &global_eid, NULL);
  if (ret != SGX_SUCCESS) {
    print_error_message(ret);
    return -1;
  }

  return 0;
}

/* OCall functions */
void ocall_print_string(const char *str) {
  /* Proxy/Bridge will check the length and null-terminate
   * the input string to prevent buffer overflow.
   */
  printf("%s", str);
}

void ocall_FAIL() {
  printf("Called FAIL oCall\n");
  exit(-1);
}

#include <chrono>

uint64_t ocall_measure_time(void) {
  // returns linux-epoch-time in nanoseconds.
  auto now_c = std::chrono::system_clock::now();
  uint64_t ret = now_c.time_since_epoch().count();
  return ret;
}

void sha256sum(const uint8_t *data, uint32_t data_size, uint8_t *hash) {
  SHA256_CTX sha256;
  SHA256_Init(&sha256);
  SHA256_Update(&sha256, data, data_size);
  SHA256_Final(hash, &sha256);
}

void printh(uint8_t *hash, size_t sz) {
  std::stringstream ss;
  for (int i = 0; i < sz; i++) {
    ss << std::hex << std::setw(2) << std::setfill('0') << (int)hash[i];
  }
  std::cout << ss.str() << std::endl;
}

bool create_app_enclave_report(const char *enclave_path,
                               sgx_target_info_t qe_target_info,
                               sgx_report_t *app_report,
                               const sgx_report_data_t *p_data);

const char *format_hex_buffer(char *buffer, uint maxSize, uint8_t *data,
                              size_t size);
const char *uint16_to_buffer(char *buffer, uint maxSize, uint16_t data,
                             size_t size);

/* Application entry */
int SGX_CDECL main(int argc, char *argv[]) {
  (void)(argc);
  (void)(argv);

  /* Initialize the enclave */
  if (initialize_enclave() < 0) {
    printf("Enter a character before exit ...\n");
    getchar();
    return -1;
  }

  /* Utilize trusted libraries */
  ActualMain();

  /* Destroy the enclave */
  sgx_destroy_enclave(global_eid);

  return 0;
}

const char *uint16_to_buffer(char *buffer, uint maxSize, uint16_t n,
                             size_t size) {
  if (size * 2 >= maxSize || size < 2) return "DEADBEEF";
  sprintf(&buffer[0], "%02X", uint8_t(n));
  sprintf(&buffer[2], "%02X", uint8_t(n >> 8));

  for (int i = 2; i < size; i++) {
    sprintf(&buffer[i * 2], "%02X", 0);
  }
  buffer[size * 2 + 1] = '\0';
  return buffer;
}

const char *format_hex_buffer(char *buffer, uint maxSize, uint8_t *data,
                              size_t size) {
  if (size * 2 >= maxSize) return "DEADBEEF";

  for (int i = 0; i < size; i++) {
    sprintf(&buffer[i * 2], "%02X", data[i]);
  }
  buffer[size * 2 + 1] = '\0';
  return buffer;
}

bool create_app_enclave_report(const char *enclave_path,
                               sgx_target_info_t qe_target_info,
                               sgx_report_t *app_report,
                               const sgx_report_data_t *p_data) {
  bool ret = true;
  uint32_t retval = 0;
  sgx_status_t sgx_status = SGX_SUCCESS;
  sgx_enclave_id_t eid = 0;
  int launch_token_updated = 0;
  sgx_launch_token_t launch_token = {0};

  sgx_status = sgx_create_enclave(enclave_path, SGX_DEBUG_FLAG, &launch_token,
                                  &launch_token_updated, &eid, NULL);
  if (SGX_SUCCESS != sgx_status) {
    printf("Error, call sgx_create_enclave fail [%s], SGXError:%04x.\n",
           __FUNCTION__, sgx_status);
    ret = false;
    sgx_destroy_enclave(eid);
    return ret;
  }

  sgx_status =
      enclave_create_report(eid, &retval, &qe_target_info, p_data, app_report);
  if ((SGX_SUCCESS != sgx_status) || (0 != retval)) {
    printf("\nCall to get_app_enclave_report() failed\n");
    ret = false;
    sgx_destroy_enclave(eid);
    return ret;
  }

  sgx_destroy_enclave(eid);
  return ret;
}
