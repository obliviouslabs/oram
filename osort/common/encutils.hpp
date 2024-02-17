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
#ifdef ENCLAVE_MODE
  br_aesctr_drbg_context ctx;
  size_t idx = 0;
  uint8_t buffer[4096];
#else
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
#ifndef NOOPENSSL
#include <openssl/evp.h>
#include <openssl/rand.h>
template <typename T>
class PRP {
  unsigned char aes_key[EVP_MAX_KEY_LENGTH];
  unsigned char iv[EVP_MAX_IV_LENGTH];

  static void generate_aes_key_and_iv(unsigned char* aes_key,
                                      unsigned char* iv) {
    RAND_bytes(aes_key, EVP_CIPHER_key_length(EVP_aes_128_cbc()));
    RAND_bytes(iv, EVP_CIPHER_iv_length(EVP_aes_128_cbc()));
  }

  static void prp_encrypt(const unsigned char* plaintext, int plaintext_len,
                          const unsigned char* key, const unsigned char* iv,
                          unsigned char* ciphertext, int* ciphertext_len) {
    EVP_CIPHER_CTX* ctx;

    // Create and initialize the context
    ctx = EVP_CIPHER_CTX_new();
    EVP_EncryptInit_ex(ctx, EVP_aes_128_cbc(), NULL, key, iv);

    // Perform encryption
    EVP_EncryptUpdate(ctx, ciphertext, ciphertext_len, plaintext,
                      plaintext_len);
    EVP_EncryptFinal_ex(ctx, ciphertext + *ciphertext_len, ciphertext_len);

    // Clean up
    EVP_CIPHER_CTX_free(ctx);
  }

 public:
  // Generate a random AES key and IV
  PRP() { generate_aes_key_and_iv(aes_key, iv); }
  T prp(T val) const {
    uint8_t outBuffer[32];
    int buffer_len;
    prp_encrypt((unsigned char*)&val, sizeof(T), aes_key, iv, outBuffer,
                &buffer_len);
    return *(T*)outBuffer;
  }
};
#else
// TODO implement later
template <typename T>
class PRP {
  uint64_t shift;

 public:
  PRP() {
    RandGen randGen;
    shift = randGen.rand64();
  }
  T prp(T val) const { return val + shift; }
};
#endif