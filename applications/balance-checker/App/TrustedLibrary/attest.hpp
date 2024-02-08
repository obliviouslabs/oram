/*
https://github.com/Azure-Samples/microsoft-azure-attestation/blob/master/sgx.attest.sample.intel.sdk/genquotes/host/host.cpp#L207
 * Copyright (C) 2011-2018 Intel Corporation. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in
 *     the documentation and/or other materials provided with the
 *     distribution.
 *   * Neither the name of Intel Corporation nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
#pragma once
#include <assert.h>
#include <openssl/sha.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "sgx_dcap_ql_wrapper.h"
#include "sgx_eid.h"   /* sgx_enclave_id_t */
#include "sgx_error.h" /* sgx_status_t */
#include "sgx_quote_3.h"
#include "sgx_urts.h"

void sha256sum(const uint8_t *data, uint32_t data_size, uint8_t *hash) {
  SHA256_CTX sha256;
  SHA256_Init(&sha256);
  SHA256_Update(&sha256, data, data_size);
  SHA256_Final(hash, &sha256);
}

bool create_app_enclave_report(const char *enclave_path,
                               sgx_target_info_t qe_target_info,
                               sgx_report_t *app_report,
                               const sgx_report_data_t *p_data);

const char *format_hex_buffer(char *buffer, uint maxSize, const uint8_t *data,
                              size_t size);
const char *uint16_to_buffer(char *buffer, uint maxSize, uint16_t data,
                             size_t size);

int get_target_info(sgx_target_info_t *target_info) {
  printf("\nStep1: Call sgx_qe_get_target_info: ");
  quote3_error_t qe3_ret = sgx_qe_get_target_info(target_info);
  if (SGX_QL_SUCCESS != qe3_ret) {
    printf("Error, call sgx_qe_get_target_info fail, ret = 0x%x\n", qe3_ret);
    return -1;
  }
  return 0;
}

void get_quote_size() {
  printf("\nStep3: Call sgx_qe_get_quote_size: ");
  uint32_t quote_size = 0;
  quote3_error_t qe3_ret = sgx_qe_get_quote_size(&quote_size);
  if (SGX_QL_SUCCESS != qe3_ret) {
    printf("Error in sgx_qe_get_quote_size. 0x%04x\n", qe3_ret);
    abort();
  }
  printf("succeed!\n");
}

uint8_t *get_quote(const sgx_report_t *app_report, uint32_t quote_size) {
  uint8_t *p_quote_buffer = (uint8_t *)malloc(quote_size);
  if (NULL == p_quote_buffer) {
    printf("\nCouldn't allocate quote_buffer\n");
    abort();
  }
  memset(p_quote_buffer, 0, quote_size);
  printf("\nStep4: Call sgx_qe_get_quote: ");
  quote3_error_t qe3_ret =
      sgx_qe_get_quote(app_report, quote_size, p_quote_buffer);
  if (SGX_QL_SUCCESS != qe3_ret) {
    printf("Error in sgx_qe_get_quote. 0x%04x\n", qe3_ret);
    abort();
  }
  printf("succeed!\n");
  return p_quote_buffer;
  //   sgx_quote3_t *p_quote = (_sgx_quote3_t *)p_quote_buffer;
  //   sgx_ql_ecdsa_sig_data_t *p_sig_data =
  //       (sgx_ql_ecdsa_sig_data_t *)p_quote->signature_data;
  //   sgx_ql_auth_data_t *p_auth_data =
  //       (sgx_ql_auth_data_t *)p_sig_data->auth_certification_data;
  //   sgx_ql_certification_data_t *p_cert_data =
  //       (sgx_ql_certification_data_t *)((uint8_t *)p_auth_data +
  //                                       sizeof(*p_auth_data) +
  //                                       p_auth_data->size);
}

void store_quote(const sgx_report_t *app_report, uint32_t quote_size,
                 uint8_t *p_quote_buffer, const uint8_t *enclave_held_data,
                 size_t enclave_held_data_size) {
  const int hex_buffer_size = 1024 * 64;
  char hex_buffer[hex_buffer_size];

  std::string output_dir("./out/");
  std::string cmd("mkdir -p " + output_dir);
  std::string file(output_dir + std::string("enclave.info.debug.json"));
  int result = system(cmd.c_str());
  printf("\nExecuted command '%s' with the result:%u", cmd.c_str(), result);
  printf(
      "\nStep5: Saving quote to JSON file "
      "name = %s\n",
      file.c_str());
  FILE *fp = fopen(file.c_str(), "w");
  fprintf(fp, "%s\n", "{");
  fprintf(fp, "  \"Type\": %d,\n", (int)2);
  // In open-enclave sdk enclave type 2 means OE_ENCLAVE_TYPE_SGX:
  // https://github.com/openenclave/openenclave/blob/3e15573418caed43f9094ff8aec36cdde4f278f7/include/openenclave/bits/types.h#L127
  fprintf(fp, "  \"MrEnclaveHex\": \"%s\",\n",
          format_hex_buffer(hex_buffer, hex_buffer_size,
                            app_report->body.mr_enclave.m, SGX_HASH_SIZE));
  fprintf(fp, "  \"MrSignerHex\": \"%s\",\n",
          format_hex_buffer(hex_buffer, hex_buffer_size,
                            app_report->body.mr_signer.m, SGX_HASH_SIZE));
  fprintf(fp, "  \"ProductIdHex\": \"%s\",\n",
          uint16_to_buffer(hex_buffer, hex_buffer_size,
                           (uint16_t)app_report->body.isv_prod_id, 16));
  fprintf(fp, "  \"SecurityVersion\": %u,\n", (int)app_report->body.isv_svn);
  fprintf(fp, "  \"Attributes\": %lu,\n",
          (uint64_t)app_report->body.attributes.flags);
  fprintf(fp, "  \"QuoteHex\": \"%s\",\n",
          format_hex_buffer(hex_buffer, hex_buffer_size, p_quote_buffer,
                            quote_size));
  fprintf(fp, "  \"EnclaveHeldDataHex\": \"%s\"\n",
          format_hex_buffer(hex_buffer, hex_buffer_size, enclave_held_data,
                            enclave_held_data_size));
  fprintf(fp, "%s\n", "}");
  fclose(fp);

  if (NULL != p_quote_buffer) {
    free(p_quote_buffer);
  }
}

sgx_report_data_t get_report_data_hash(const uint8_t *public_key_data,
                                       size_t public_key_data_size) {
  sgx_report_data_t hash;
  sha256sum(public_key_data, public_key_data_size, hash.d);
  printf("succeed!\n");
  return hash;
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

const char *format_hex_buffer(char *buffer, uint maxSize, const uint8_t *data,
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

  sgx_status = enclave_create_report(global_eid, &retval, &qe_target_info,
                                     p_data, app_report);
  if ((SGX_SUCCESS != sgx_status) || (0 != retval)) {
    printf("\nCall to get_app_enclave_report() failed\n");
    ret = false;
  }

  return ret;
}