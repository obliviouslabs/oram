#pragma once

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include <string>

#include "sgx_eid.h"
#include "sgx_error.h"

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#if defined(__GNUC__)
#define TOKEN_FILENAME "enclave.token"
#define ENCLAVE_FILENAME "enclave.signed.so"
#endif

extern sgx_enclave_id_t global_eid; /* global enclave id */

#if defined(__cplusplus)
extern "C" {
#endif

std::string GetPublicKeyBase64(void);
void ActualMain(const char* dbPath);

#if defined(__cplusplus)
}
#endif
