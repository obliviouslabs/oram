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
                         uint8_t iv[AES_BLOCK_SIZE], uint8_t* ciphertext);
bool aes_256_ctr_decrypt(uint64_t ciphertextSize, uint8_t* ciphertext,
                         const uint8_t key[AES_BLOCK_SIZE],
                         uint8_t iv[AES_BLOCK_SIZE], uint8_t* plaintext);

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
#ifdef ENCLAVE_MODE
  br_aesctr_drbg_context ctx;
  size_t idx = 0;
  uint8_t buffer[4096];
#else
  std::minstd_rand engine;
#endif
 public:
  RandGen();
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