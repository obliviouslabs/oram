#pragma once
#include <string>

#include "../../common.hpp"
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