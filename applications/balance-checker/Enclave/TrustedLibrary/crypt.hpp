#pragma once
#include <string>

#include "../../common.hpp"
#include "sgx_tcrypto.h"
#include "sgx_tseal.h"
sgx_status_t generate_key_pair(sgx_ec256_private_t* p_private,
                               sgx_ec256_public_t* p_public) {
  sgx_ecc_state_handle_t ecc_handle;
  sgx_status_t status = sgx_ecc256_open_context(&ecc_handle);
  if (status != SGX_SUCCESS) {
    return status;
  }

  status = sgx_ecc256_create_key_pair(p_private, p_public, ecc_handle);
  sgx_ecc256_close_context(ecc_handle);

  return status;
}

sgx_status_t compute_shared_key(const sgx_ec256_private_t* p_private_b,
                                const sgx_ec256_public_t* p_public_a,
                                sgx_ec256_dh_shared_t* p_shared_key) {
  sgx_ecc_state_handle_t ecc_handle;
  sgx_status_t status = sgx_ecc256_open_context(&ecc_handle);
  if (status != SGX_SUCCESS) {
    return status;
  }

  status = sgx_ecc256_compute_shared_dhkey(p_private_b, p_public_a,
                                           p_shared_key, ecc_handle);
  sgx_ecc256_close_context(ecc_handle);

  return status;
}

ec256_public_t convert_to_ec256_public_t(const sgx_ec256_public_t& public_key) {
  ec256_public_t res;
  for (int i = 0; i < 32; ++i) {
    res.gx[i] = public_key.gx[31 - i];
    res.gy[i] = public_key.gy[31 - i];
  }
  return res;
}

sgx_ec256_public_t convert_to_sgx_ec256_public_t(
    const ec256_public_t& public_key_big_endian) {
  sgx_ec256_public_t res;
  for (int i = 0; i < 32; ++i) {
    res.gx[i] = public_key_big_endian.gx[31 - i];
    res.gy[i] = public_key_big_endian.gy[31 - i];
  }
  return res;
}

typedef sgx_ec256_private_t private_key_t;

sgx_status_t seal_private_key(const private_key_t* private_key,
                              uint8_t** sealed_blob_out,
                              size_t* sealed_size_out) {
  sgx_status_t status;

  // Specify the key policy (e.g., MRENCLAVE or MRSIGNER) and the enclave's
  // execution attributes
  uint16_t key_policy = SGX_KEYPOLICY_MRENCLAVE;
  sgx_attributes_t attribute_mask = {SGX_FLAGS_INITTED | SGX_FLAGS_DEBUG, 0};
  sgx_misc_select_t misc_mask = 0;
  size_t private_key_size = sizeof(private_key_t);

  // Calculate the size of the sealed data
  uint32_t sealed_data_size = sgx_calc_sealed_data_size(0, private_key_size);
  if (sealed_data_size == UINT32_MAX) {
    return SGX_ERROR_UNEXPECTED;
  }

  // Allocate space for the sealed data
  sgx_sealed_data_t* sealed_data = (sgx_sealed_data_t*)malloc(sealed_data_size);
  if (!sealed_data) {
    return SGX_ERROR_OUT_OF_MEMORY;
  }

  // Seal the private key
  status = sgx_seal_data_ex(key_policy, attribute_mask, misc_mask, 0, NULL,
                            private_key_size, (uint8_t*)private_key,
                            sealed_data_size, sealed_data);
  if (status != SGX_SUCCESS) {
    free(sealed_data);
    return status;
  }

  // Set the output parameters
  *sealed_blob_out = (uint8_t*)sealed_data;
  *sealed_size_out = sealed_data_size;

  return SGX_SUCCESS;
}

sgx_status_t unseal_private_key(const uint8_t* sealed_data,
                                size_t sealed_data_size, private_key_t* key) {
  sgx_status_t sgx_status;
  const sgx_sealed_data_t* sealed_blob = (const sgx_sealed_data_t*)sealed_data;
  uint32_t decrypted_data_len = sgx_get_encrypt_txt_len(sealed_blob);

  if (decrypted_data_len != sizeof(private_key_t)) return SGX_ERROR_UNEXPECTED;

  // Unseal the private key
  sgx_status = sgx_unseal_data(sealed_blob, NULL, NULL, (uint8_t*)key,
                               &decrypted_data_len);
  if (SGX_SUCCESS != sgx_status) {
    return sgx_status;
  }

  return SGX_SUCCESS;
}