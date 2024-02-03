#pragma once
#include <string>

#include "sgx_tcrypto.h"
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