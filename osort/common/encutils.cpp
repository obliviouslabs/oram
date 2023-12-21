#include <utility>
#ifdef ENCLAVE_MODE
// #include "../Enclave.h"
#include "Enclave_t.h"
#endif

#ifndef NOOPENSSL
#include <openssl/aes.h>
#include <openssl/crypto.h>
#include <openssl/evp.h>
#include <openssl/rand.h>

#include "common/encutils.hpp"

RandGen default_rand;

void aes_init() {
  static int init = 0;
  if (init == 0) {
    OpenSSL_add_all_ciphers();
    // int rv = RAND_load_file("/dev/urandom", 32);
    init = 1;
  }
}

void aes_256_gcm_encrypt(uint64_t plaintextSize, uint8_t *plaintext,
                         const uint8_t key[AES_BLOCK_SIZE],
                         uint8_t iv[AES_BLOCK_SIZE],
                         uint8_t tag[AES_BLOCK_SIZE], uint8_t *ciphertext) {
  aes_init();

  size_t enc_length = ((plaintextSize + 15) / 16) * 16;

  int actual_size = 0, final_size = 0;
  EVP_CIPHER_CTX *e_ctx = EVP_CIPHER_CTX_new();
  EVP_CIPHER_CTX_ctrl(e_ctx, EVP_CTRL_GCM_SET_IVLEN, 16, NULL);
  EVP_CIPHER_CTX_set_padding(e_ctx, 0);
  EVP_EncryptInit(e_ctx, EVP_aes_256_gcm(), key, iv);

  EVP_EncryptUpdate(e_ctx, ciphertext, &actual_size, plaintext, plaintextSize);
  int ok = EVP_EncryptFinal(e_ctx, &ciphertext[actual_size], &final_size);
  Assert(final_size <= enc_length);
  EVP_CIPHER_CTX_ctrl(e_ctx, EVP_CTRL_GCM_GET_TAG, 16, tag);
  EVP_CIPHER_CTX_free(e_ctx);
  Assert(ok == 1);
  IGNORE_UNUSED(enc_length);
  IGNORE_UNUSED(ok);
}

bool aes_256_gcm_decrypt(uint64_t ciphertextSize, uint8_t *ciphertext,
                         const uint8_t key[AES_BLOCK_SIZE],
                         uint8_t iv[AES_BLOCK_SIZE],
                         uint8_t tag[AES_BLOCK_SIZE], uint8_t *plaintext) {
  aes_init();

  int actual_size = 0, final_size = 0;
  EVP_CIPHER_CTX *d_ctx = EVP_CIPHER_CTX_new();
  EVP_CIPHER_CTX_ctrl(d_ctx, EVP_CTRL_GCM_SET_IVLEN, 16, NULL);
  EVP_CIPHER_CTX_set_padding(d_ctx, 0);
  EVP_DecryptInit(d_ctx, EVP_aes_256_gcm(), key, iv);
  EVP_DecryptUpdate(d_ctx, plaintext, &actual_size, ciphertext, ciphertextSize);
  EVP_CIPHER_CTX_ctrl(d_ctx, EVP_CTRL_GCM_SET_TAG, 16, tag);
  int ok;
  ok = EVP_DecryptFinal(d_ctx, &plaintext[actual_size], &final_size);
  EVP_CIPHER_CTX_free(d_ctx);

  Assert(ok == 1);
  return ok == 1;
}

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

#else
#include <x86intrin.h>

#include "bearssl_aead.h"
#include "bearssl_hash.h"
#include "common/encutils.hpp"
#define memset_s(s, smax, c, n) memset(s, c, n);
#define MAX_ENC_SIZE 8192

#define SGXSD_AES_GCM_IV_SIZE 12
#define SGXSD_AES_GCM_MAC_SIZE 16
#define SGXSD_AES_GCM_KEY_SIZE 32
#define SGXSD_CURVE25519_KEY_SIZE 32
#define SGXSD_SHA256_HASH_SIZE 32

RandGen default_rand;

// In enclave mode we can't use openssl, so we use libsodium, which should be
// packed in the enclave directly.
void aes_init() {
  static int init = 0;
  if (init == 0) {
    // Assert(crypto_aead_aes256gcm_is_available());
    // int rv = RAND_load_file("/dev/urandom", 32);
    init = 1;
  }
}

void __attribute__((noinline)) sgxsd_br_clear_stack() {
  uint8_t stack[4096];
  memset_s(&stack, sizeof(stack), 0, sizeof(stack));
  _mm256_zeroall();
}

sgx_status_t sgxsd_aes_gcm_run(bool encrypt,
                               const uint8_t p_key[SGXSD_AES_GCM_KEY_SIZE],
                               const void *p_src, uint32_t src_len, void *p_dst,
                               const uint8_t p_iv[SGXSD_AES_GCM_KEY_SIZE],
                               const void *p_aad, uint32_t aad_len,
                               uint8_t p_mac[SGXSD_AES_GCM_KEY_SIZE]) {
  if (p_key == NULL || ((p_src == NULL || p_dst == NULL) && src_len != 0) ||
      p_iv == NULL || (p_aad == NULL && aad_len != 0) || p_mac == NULL) {
    return SGX_ERROR_INVALID_PARAMETER;
  }
  br_aes_x86ni_ctr_keys aes_ctx;
  br_aes_x86ni_ctr_init(&aes_ctx, p_key, SGXSD_AES_GCM_KEY_SIZE);
  br_gcm_context aes_gcm_ctx;
  br_gcm_init(&aes_gcm_ctx, &aes_ctx.vtable, &br_ghash_pclmul);
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
  memset_s(&aes_ctx, sizeof(aes_ctx), 0, sizeof(aes_ctx));
  memset_s(&aes_gcm_ctx, sizeof(aes_gcm_ctx), 0, sizeof(aes_gcm_ctx));
  if (tag_res) {
    return SGX_SUCCESS;
  } else {
    if (p_dst != NULL) {
      memset_s(p_dst, src_len, 0, src_len);
    }
    return SGX_ERROR_MAC_MISMATCH;
  }
}

void aes_256_gcm_encrypt(uint64_t plaintextSize, uint8_t *plaintext,
                         const uint8_t key[AES_BLOCK_SIZE],
                         uint8_t iv[AES_BLOCK_SIZE],
                         uint8_t tag[AES_BLOCK_SIZE], uint8_t *ciphertext) {
  aes_init();
  sgxsd_aes_gcm_run(true, key, plaintext, plaintextSize, ciphertext, iv,
                    nullptr, 0, tag);
}

bool aes_256_gcm_decrypt(uint64_t ciphertextSize, uint8_t *ciphertext,
                         const uint8_t key[AES_BLOCK_SIZE],
                         uint8_t iv[AES_BLOCK_SIZE],
                         uint8_t tag[AES_BLOCK_SIZE], uint8_t *plaintext) {
  aes_init();
  return 0 == sgxsd_aes_gcm_run(false, key, ciphertext, ciphertextSize,
                                plaintext, iv, nullptr, 0, tag);
}

sgx_status_t sgxsd_aes_ctr_run(bool encrypt,
                               const uint8_t p_key[SGXSD_AES_GCM_KEY_SIZE],
                               const void *p_src, uint32_t src_len, void *p_dst,
                               const uint8_t p_iv[SGXSD_AES_GCM_KEY_SIZE]) {
  if (p_key == NULL || ((p_src == NULL || p_dst == NULL) && src_len != 0) ||
      p_iv == NULL) {
    return SGX_ERROR_INVALID_PARAMETER;
  }
  br_aes_x86ni_ctr_keys aes_ctx;
  br_aes_x86ni_ctr_init(&aes_ctx, p_key, SGXSD_AES_GCM_KEY_SIZE);
  if (src_len != 0) {
    memmove(p_dst, p_src, src_len);
    br_aes_x86ni_ctr_run(&aes_ctx, p_iv, 0, p_dst, src_len);
  }

  // sgxsd_br_clear_stack();
  // memset_s(&aes_ctx, sizeof(aes_ctx), 0, sizeof(aes_ctx));
  return SGX_SUCCESS;
}

void aes_256_ctr_encrypt(uint64_t plaintextSize, uint8_t *plaintext,
                         const uint8_t key[AES_BLOCK_SIZE],
                         const uint8_t iv[AES_BLOCK_SIZE],
                         uint8_t *ciphertext) {
  aes_init();
  sgxsd_aes_ctr_run(true, key, plaintext, plaintextSize, ciphertext, iv);
}

bool aes_256_ctr_decrypt(uint64_t ciphertextSize, uint8_t *ciphertext,
                         const uint8_t key[AES_BLOCK_SIZE],
                         const uint8_t iv[AES_BLOCK_SIZE], uint8_t *plaintext) {
  aes_init();
  return 0 == sgxsd_aes_ctr_run(false, key, ciphertext, ciphertextSize,
                                plaintext, iv);
}

RandGen::RandGen() { new (this) RandGen(rd()); }

RandGen::RandGen(uint64_t seed) {
  const br_block_ctr_class *aes_vtable =
      &br_aes_x86ni_ctr_vtable;  // in theory we should use
                                 // br_aes_x86ni_ctr_get_vtable() here, but it
                                 // seems bearssl does not recognize the cpuid
                                 // and always return NUll
  // if (aes_vtable == NULL) {
  //   printf("vtable fail\n");
  //   throw 1;
  // }
  br_aesctr_drbg_init(&ctx, aes_vtable, &seed, 8);
  br_aesctr_drbg_generate(&ctx, buffer, 4096);
}

uint64_t RandGen::rand64() {
  // Generate random output
  if (idx > 4088) {
    br_aesctr_drbg_generate(&ctx, buffer, 4096);
    idx = 0;
  }
  uint64_t output = *(uint64_t *)(buffer + idx);
  idx += 8;
  return output;
}

uint32_t RandGen::rand32() {
  if (idx > 4092) {
    br_aesctr_drbg_generate(&ctx, buffer, 4096);
    idx = 0;
  }
  uint32_t output = *(uint32_t *)(buffer + idx);
  idx += 4;
  return output;
}
uint8_t RandGen::rand1() {
  if (idx > 4095) {
    br_aesctr_drbg_generate(&ctx, buffer, 4096);
    idx = 0;
  }
  uint8_t output = *(uint8_t *)(buffer + idx);
  idx += 1;
  return output & 1;
}

#endif