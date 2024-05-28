#include <stdexcept>
#include <utility>
#ifdef ENCLAVE_MODE_ENCLAVE
// #include "../Enclave.h"
#include "Enclave_t.h"
#endif

#include <x86intrin.h>

#ifndef ENCLAVE_MODE
#include "bearssl_aead.h"
#include "bearssl_hash.h"
#endif

#include "common/encutils.hpp"
#define memset_s(s, smax, c, n) memset(s, c, n);

RandGen default_rand;
// random AES key:
uint8_t KEY[AES_BLOCK_SIZE] = {0};
br_aes_x86ni_ctr_keys aes_ctx;

// We have two versions of bearssl. For the enclave we need the sgx prepacked
// one, for non enclave we need one installed in the OS. This file allows us to
// use SGX library functions instead of bearssl inside of the enclave whenever
// they are available.

#ifndef ENCLAVE_MODE_ENCLAVE
void handleErrors(void) { throw std::runtime_error("BearSSL error"); }
#endif

void __attribute__((noinline)) sgxsd_br_clear_stack() {
  uint8_t stack[4096];
  memset_s(&stack, sizeof(stack), 0, sizeof(stack));
  _mm256_zeroall();
}

#ifndef ENCLAVE_MODE_ENCLAVE
uint64_t secure_hash_with_salt(const uint8_t *data, size_t data_size,
                               const uint8_t (&salt)[16]) {
  uint64_t res;
  // Initialize the hash context
  br_sha256_context ctx;
  br_sha256_init(&ctx);

  // Hash the salt
  br_sha256_update(&ctx, salt, sizeof(salt));

  // Hash the data
  br_sha256_update(&ctx, data, data_size);

  // Finalize the hash and get the result
  unsigned char
      hash[br_sha256_SIZE];  // br_sha256_SIZE is normally 32 for SHA-256
  br_sha256_out(&ctx, hash);

  // Copy the result to the output buffer, up to resSize bytes
  memcpy(&res, hash, 8);

  return res;
}
#else
#include <cstring>

#include "sgx_tcrypto.h"
#include "sgx_trts.h"
uint64_t secure_hash_with_salt(const uint8_t *data, size_t data_size,
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

bool sgxsd_aes_gcm_run(bool encrypt, const void *p_src, uint32_t src_len,
                       void *p_dst, const uint8_t p_iv[SGXSD_AES_GCM_IV_SIZE],
                       const void *p_aad, uint32_t aad_len,
                       uint8_t p_mac[SGXSD_AES_GCM_MAC_SIZE],
                       br_aes_x86ni_ctr_keys *aes_ctx_ptr) {
  if (((p_src == NULL || p_dst == NULL) && src_len != 0) || p_iv == NULL ||
      (p_aad == NULL && aad_len != 0) || p_mac == NULL) {
    return 0;
  }
  br_gcm_context aes_gcm_ctx;
  br_gcm_init(&aes_gcm_ctx, &(aes_ctx_ptr->vtable), &br_ghash_pclmul);
  br_gcm_reset(&aes_gcm_ctx, p_iv, SGXSD_AES_GCM_IV_SIZE);
  if (aad_len != 0) {
    br_gcm_aad_inject(&aes_gcm_ctx, p_aad, aad_len);
  }
  br_gcm_flip(&aes_gcm_ctx);
  if (src_len != 0) {
    memmove(p_dst, p_src, src_len);
    br_gcm_run(&aes_gcm_ctx, encrypt, p_dst, src_len);
  }
  bool tag_res;
  if (encrypt) {
    br_gcm_get_tag(&aes_gcm_ctx, p_mac);
    tag_res = true;
  } else {
    tag_res = br_gcm_check_tag(&aes_gcm_ctx, p_mac);
  }
  sgxsd_br_clear_stack();
  memset_s(&aes_gcm_ctx, sizeof(aes_gcm_ctx), 0, sizeof(aes_gcm_ctx));
  if (tag_res) {
    return 1;
  } else {
    if (p_dst != NULL) {
      memset_s(p_dst, src_len, 0, src_len);
    }
    return 0;
  }
}

void aes_256_gcm_encrypt(uint64_t plaintextSize, uint8_t *plaintext,
                         const uint8_t iv[SGXSD_AES_GCM_IV_SIZE],
                         uint8_t tag[SGXSD_AES_GCM_MAC_SIZE],
                         uint8_t *ciphertext,
                         br_aes_x86ni_ctr_keys *aes_ctx_ptr) {
  sgxsd_aes_gcm_run(true, plaintext, plaintextSize, ciphertext, iv, nullptr, 0,
                    tag, aes_ctx_ptr);
}

bool aes_256_gcm_decrypt(uint64_t ciphertextSize, uint8_t *ciphertext,
                         const uint8_t iv[SGXSD_AES_GCM_IV_SIZE],
                         uint8_t tag[SGXSD_AES_GCM_MAC_SIZE],
                         uint8_t *plaintext,
                         br_aes_x86ni_ctr_keys *aes_ctx_ptr) {
  return sgxsd_aes_gcm_run(false, ciphertext, ciphertextSize, plaintext, iv,
                           nullptr, 0, tag, aes_ctx_ptr);
}

bool sgxsd_aes_ctr_run(bool encrypt, const void *p_src, uint32_t src_len,
                       void *p_dst, const uint8_t p_iv[SGXSD_AES_GCM_IV_SIZE],
                       br_aes_x86ni_ctr_keys *aes_ctx_ptr) {
  if (((p_src == NULL || p_dst == NULL) && src_len != 0) || p_iv == NULL) {
    return 0;
  }
  if (src_len != 0) {
    memmove(p_dst, p_src, src_len);
    br_aes_x86ni_ctr_run(aes_ctx_ptr, p_iv, 0, p_dst, src_len);
  }
  return 1;
}

void aes_256_ctr_encrypt(uint64_t plaintextSize, uint8_t *plaintext,
                         const uint8_t iv[SGXSD_AES_GCM_IV_SIZE],
                         uint8_t *ciphertext,
                         br_aes_x86ni_ctr_keys *aes_ctx_ptr) {
  sgxsd_aes_ctr_run(true, plaintext, plaintextSize, ciphertext, iv,
                    aes_ctx_ptr);
}

bool aes_256_ctr_decrypt(uint64_t ciphertextSize, uint8_t *ciphertext,
                         const uint8_t iv[SGXSD_AES_GCM_IV_SIZE],
                         uint8_t *plaintext,
                         br_aes_x86ni_ctr_keys *aes_ctx_ptr) {
  return sgxsd_aes_ctr_run(false, ciphertext, ciphertextSize, plaintext, iv,
                           aes_ctx_ptr);
}

#ifndef ENCLAVE_MODE_ENCLAVE
RandGen::RandGen() { new (this) RandGen(rd()); }
RandGen::RandGen(uint64_t seed) : engine(seed) {}
uint64_t RandGen::rand64() {
  std::uniform_int_distribution<uint64_t> d;
  return d(engine);
}
uint32_t RandGen::rand32() {
  std::uniform_int_distribution<uint32_t> d;
  return d(engine);
}
uint8_t RandGen::rand1() {
  std::uniform_int_distribution<short> d(0, 1);
  return d(engine);
}

void read_rand(uint8_t *output, size_t size) {
  FILE *fp = fopen("/dev/urandom", "rb");
  if (fp == NULL) {
    perror("Failed to open /dev/urandom");
    return;  // Failure
  }

  size_t read = fread(output, 1, size, fp);
  fclose(fp);

  if (read != size) {
    perror("Failed to read enough bytes");
    // Handle the error, not enough data was read
    return;  // Failure
  }

  return;  // Success
}

#else
RandGen::RandGen() {}

uint64_t RandGen::rand64() {
  uint64_t output;
  sgx_read_rand((uint8_t *)&output, sizeof(output));
  return output;
}

uint32_t RandGen::rand32() {
  uint32_t output;
  sgx_read_rand((uint8_t *)&output, sizeof(output));
  return output;
}
uint8_t RandGen::rand1() {
  uint8_t output;
  sgx_read_rand((uint8_t *)&output, sizeof(output));
  return output & 1;
}

void read_rand(uint8_t *output, size_t size) { sgx_read_rand(output, size); }

#endif
struct GlobalRandomKeySetter {
  GlobalRandomKeySetter() {
    read_rand(KEY, AES_BLOCK_SIZE);
    br_aes_x86ni_ctr_init(&aes_ctx, KEY, SGXSD_AES_GCM_KEY_SIZE);
  }
};

GlobalRandomKeySetter global_random_key_setter;