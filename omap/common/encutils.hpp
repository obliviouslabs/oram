#pragma once
#include <cassert>
#include <random>

#include "common/defs.hpp"
#ifdef ENCLAVE_MODE
#include "bearssl_rand.h"
#include "sgx_trts.h"
#endif
#ifndef AES_BLOCK_SIZE
#define AES_BLOCK_SIZE 32
#endif
void aes_init();
void aes_256_gcm_encrypt(uint64_t plaintextSize, uint8_t* plaintext,
                         const uint8_t key[AES_BLOCK_SIZE],
                         uint8_t iv[AES_BLOCK_SIZE],
                         uint8_t tag[AES_BLOCK_SIZE], uint8_t* ciphertext);
bool aes_256_gcm_decrypt(uint64_t ciphertextSize, uint8_t* ciphertext,
                         const uint8_t key[AES_BLOCK_SIZE],
                         uint8_t iv[AES_BLOCK_SIZE],
                         uint8_t tag[AES_BLOCK_SIZE], uint8_t* plaintext);

void aes_256_ctr_encrypt(uint64_t plaintextSize, uint8_t* plaintext,
                         const uint8_t key[AES_BLOCK_SIZE],
                         const uint8_t iv[AES_BLOCK_SIZE], uint8_t* ciphertext);
bool aes_256_ctr_decrypt(uint64_t ciphertextSize, uint8_t* ciphertext,
                         const uint8_t key[AES_BLOCK_SIZE],
                         const uint8_t iv[AES_BLOCK_SIZE], uint8_t* plaintext);

template <class result_type_param = uint64_t>
class RandomDevice : public std::random_device {
 public:
  using result_type = result_type_param;

  static constexpr result_type max() {
    return std::numeric_limits<result_type>::max();
  }
  static constexpr result_type min() {
    return std::numeric_limits<result_type>::min();
  }
#ifdef ENCLAVE_MODE
  result_type operator()() {
    result_type val;
    if (sgx_read_rand((unsigned char*)&val, sizeof(result_type)) ==
        SGX_SUCCESS) {
      return val;
    }
    throw 0;
  }
#endif
};

class RandGen {
  RandomDevice<> rd;
#ifndef ENCLAVE_MODE
  std::minstd_rand engine;
#endif
 public:
  RandGen();
  RandGen(uint64_t seed);
  uint64_t rand64();
  uint32_t rand32();
  uint8_t rand1();
};

#ifdef ENCLAVE_MODE
#include "common/encutils.cpp"
#endif

// "random" (:P) AES key:
inline constexpr uint8_t KEY[16] = {0x41, 0x41, 0x41, 0x41, 0x41, 0x41,
                                    0x41, 0x41, 0x41, 0x41, 0x41, 0x41,
                                    0x41, 0x41, 0x41, 0x41};

void read_rand(uint8_t* output, size_t size);
uint64_t secure_hash_with_salt(const uint8_t* data, size_t data_size,
                               const uint8_t (&salt)[16]);

template <typename T>
uint64_t secure_hash_with_salt(const T& data, const uint8_t (&salt)[16]);

template <typename T>
void secure_hash_with_salt(const T& data, const uint8_t (&salt)[16], void* res,
                           uint8_t resSize);
#ifndef ENCLAVE_MODE_ENCLAVE
#include <bearssl.h>
#include <cstdint>
#include <cstring>

template <typename T>
uint64_t secure_hash_with_salt(const T& data, const uint8_t (&salt)[16]) {
  return secure_hash_with_salt((const uint8_t*)&data, sizeof(T), salt);
}

template <typename T>
void secure_hash_with_salt(const T& data, const uint8_t (&salt)[16], uint8_t* res, size_t resSize) {
    // Initialize the hash context
    br_sha256_context ctx;
    br_sha256_init(&ctx);

    // Hash the salt
    br_sha256_update(&ctx, salt, sizeof(salt));

    // Hash the data
    br_sha256_update(&ctx, &data, sizeof(T));

    // Finalize the hash and get the result
    unsigned char hash[br_sha256_SIZE]; // br_sha256_SIZE is normally 32 for SHA-256
    br_sha256_out(&ctx, hash);

    // Copy the result to the output buffer, up to resSize bytes
    memcpy(res, hash, resSize < br_sha256_SIZE ? resSize : br_sha256_SIZE);
}

#else
#include <cstring>

#include "sgx_tcrypto.h"
#include "sgx_trts.h"
template <typename T>
uint64_t secure_hash_with_salt(const T& data, const uint8_t (&salt)[16]) {
  uint8_t data_with_salt[sizeof(T) + 16];
  memcpy(data_with_salt, &data, sizeof(T));
  memcpy(data_with_salt + sizeof(T), salt, 16);
  sgx_sha256_hash_t hash;
  sgx_sha256_msg((uint8_t*)&data_with_salt[0], sizeof(T) + 16, &hash);
  uint64_t result;
  memcpy(&result, &hash, sizeof(result));
  return result;
}

template <typename T>
void secure_hash_with_salt(const T& data, const uint8_t (&salt)[16], void* res,
                           uint8_t resSize) {
  Assert(resSize < sizeof(sgx_sha256_hash_t));
  uint8_t data_with_salt[sizeof(T) + 16];
  memcpy(data_with_salt, &data, sizeof(T));
  memcpy(data_with_salt + sizeof(T), salt, 16);
  sgx_sha256_hash_t hash;
  sgx_sha256_msg((uint8_t*)&data_with_salt[0], sizeof(T) + 16, &hash);
  memcpy(res, &hash, resSize);
}

uint64_t secure_hash_with_salt(const uint8_t* data, size_t data_size,
                               const uint8_t (&salt)[16]) {
  sgx_sha_state_handle_t sha_handle;
  sgx_sha256_hash_t hash;
  uint64_t result = 0;

  sgx_status_t status = sgx_sha256_init(&sha_handle);
  if (status != SGX_SUCCESS) {
    // Handle error
    return 0;
  }

  // Hash the salt
  status = sgx_sha256_update(salt, sizeof(salt), sha_handle);
  if (status != SGX_SUCCESS) {
    // Handle error
    sgx_sha256_close(sha_handle);
    return 0;
  }

  // Hash the data
  status = sgx_sha256_update(data, data_size, sha_handle);
  if (status != SGX_SUCCESS) {
    // Handle error
    sgx_sha256_close(sha_handle);
    return 0;
  }

  // Finalize the hash
  status = sgx_sha256_get_hash(sha_handle, &hash);
  sgx_sha256_close(sha_handle);
  if (status != SGX_SUCCESS) {
    // Handle error
    return 0;
  }

  // Use the first 8 bytes of the hash as the result
  memcpy(&result, hash, sizeof(result));

  return result;
}

#endif