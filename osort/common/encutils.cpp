#include <utility>
#ifdef ENCLAVE_MODE
// #include "../Enclave.h"
#include "Enclave_t.h"
#endif

#ifndef NOOPENSSL
// #include <openssl/aes.h>
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

void handleErrors(void) {
  printf("AES Error\n");
  abort();
}

void aes_256_gcm_encrypt(uint64_t plaintextSize, uint8_t *plaintext,
                         const uint8_t key[AES_BLOCK_SIZE],
                         uint8_t iv[AES_BLOCK_SIZE],
                         uint8_t tag[AES_BLOCK_SIZE], uint8_t *ciphertext) {
  EVP_CIPHER_CTX *ctx;

  int len;

  int ciphertext_len;

  /* Create and initialise the context */
  if (!(ctx = EVP_CIPHER_CTX_new())) handleErrors();

  /* Initialise the encryption operation. */
  if (1 != EVP_EncryptInit_ex(ctx, EVP_aes_256_gcm(), NULL, NULL, NULL))
    handleErrors();

  // /*
  //  * Set IV length if default 12 bytes (96 bits) is not appropriate
  //  */
  // if (1 != EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_SET_IVLEN, iv_len, NULL))
  //   handleErrors();

  /* Initialise key and IV */
  if (1 != EVP_EncryptInit_ex(ctx, NULL, NULL, key, iv)) handleErrors();

  // /*
  //  * Provide any AAD data. This can be called zero or more times as
  //  * required
  //  */
  // if (1 != EVP_EncryptUpdate(ctx, NULL, &len, aad, aad_len)) handleErrors();

  /*
   * Provide the message to be encrypted, and obtain the encrypted output.
   * EVP_EncryptUpdate can be called multiple times if necessary
   */
  if (1 != EVP_EncryptUpdate(ctx, ciphertext, &len, plaintext, plaintextSize))
    handleErrors();
  ciphertext_len = len;

  /*
   * Finalise the encryption. Normally ciphertext bytes may be written at
   * this stage, but this does not occur in GCM mode
   */
  if (1 != EVP_EncryptFinal_ex(ctx, ciphertext + len, &len)) handleErrors();
  ciphertext_len += len;

  /* Get the tag */
  if (1 != EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_GET_TAG, 16, tag))
    handleErrors();

  /* Clean up */
  EVP_CIPHER_CTX_free(ctx);
}

bool aes_256_gcm_decrypt(uint64_t ciphertextSize, uint8_t *ciphertext,
                         const uint8_t key[AES_BLOCK_SIZE],
                         uint8_t iv[AES_BLOCK_SIZE],
                         uint8_t tag[AES_BLOCK_SIZE], uint8_t *plaintext) {
  EVP_CIPHER_CTX *ctx;
  int len;
  int plaintext_len;
  int ret;

  /* Create and initialise the context */
  if (!(ctx = EVP_CIPHER_CTX_new())) handleErrors();

  /* Initialise the decryption operation. */
  if (!EVP_DecryptInit_ex(ctx, EVP_aes_256_gcm(), NULL, NULL, NULL))
    handleErrors();

  // /* Set IV length. Not necessary if this is 12 bytes (96 bits) */
  // if (!EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_SET_IVLEN, iv_len, NULL))
  //   handleErrors();

  /* Initialise key and IV */
  if (!EVP_DecryptInit_ex(ctx, NULL, NULL, key, iv)) handleErrors();

  // /*
  //  * Provide any AAD data. This can be called zero or more times as
  //  * required
  //  */
  // if (!EVP_DecryptUpdate(ctx, NULL, &len, aad, aad_len)) handleErrors();

  /*
   * Provide the message to be decrypted, and obtain the plaintext output.
   * EVP_DecryptUpdate can be called multiple times if necessary
   */
  if (!EVP_DecryptUpdate(ctx, plaintext, &len, ciphertext, ciphertextSize))
    handleErrors();
  plaintext_len = len;

  /* Set expected tag value. Works in OpenSSL 1.0.1d and later */
  if (!EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_SET_TAG, 16, tag)) handleErrors();

  /*
   * Finalise the decryption. A positive return value indicates success,
   * anything else is a failure - the plaintext is not trustworthy.
   */
  ret = EVP_DecryptFinal_ex(ctx, plaintext + len, &len);

  /* Clean up */
  EVP_CIPHER_CTX_free(ctx);

  // if (ret > 0) {
  //   /* Success */
  //   plaintext_len += len;
  //   return plaintext_len;
  // } else {
  //   /* Verify failed */
  //   return -1;
  // }
  return ret > 0;
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

#include <openssl/evp.h>

#include <cstdint>
#include <cstring>

uint64_t secure_hash_with_salt(const uint8_t *data, size_t data_size,
                               const uint8_t (&salt)[16]) {
  EVP_MD_CTX *mdctx = EVP_MD_CTX_new();
  uint64_t result = 0;
  unsigned char hash[EVP_MAX_MD_SIZE];
  unsigned int lengthOfHash = 0;

  if (mdctx == NULL) {
    // Handle error
    return 0;
  }

  if (1 != EVP_DigestInit_ex(mdctx, EVP_sha256(), NULL)) {
    // Handle error
    EVP_MD_CTX_free(mdctx);
    return 0;
  }

  // Hash the salt
  if (1 != EVP_DigestUpdate(mdctx, salt, sizeof(salt))) {
    // Handle error
    EVP_MD_CTX_free(mdctx);
    return 0;
  }

  // Hash the data
  if (1 != EVP_DigestUpdate(mdctx, data, data_size)) {
    // Handle error
    EVP_MD_CTX_free(mdctx);
    return 0;
  }

  // Finalize the hash
  if (1 != EVP_DigestFinal_ex(mdctx, hash, &lengthOfHash)) {
    // Handle error
    EVP_MD_CTX_free(mdctx);
    return 0;
  }

  // Use the first 8 bytes of the hash as the result
  memcpy(&result, hash, sizeof(result));

  EVP_MD_CTX_free(mdctx);

  return result;
}

void read_rand(uint8_t *output, size_t size) {
  if (1 != RAND_bytes(output, size)) {
    printf("RAND_bytes failed\n");
  }
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

RandGen::RandGen() {
#if TCS_NUM <= 1
  new (this) RandGen(rd());
#endif
}

RandGen::RandGen(uint64_t seed) {
#if TCS_NUM <= 1
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
#endif
}

uint64_t RandGen::rand64() {
#if TCS_NUM > 1
  uint64_t output;
  sgx_read_rand((uint8_t *)&output, sizeof(output));
  return output;
#else
  // Generate random output
  if (idx > 4088) {
    br_aesctr_drbg_generate(&ctx, buffer, 4096);
    idx = 0;
  }
  uint64_t output = *(uint64_t *)(buffer + idx);
  idx += 8;
  return output;
#endif
}

uint32_t RandGen::rand32() {
#if TCS_NUM > 1
  uint32_t output;
  sgx_read_rand((uint8_t *)&output, sizeof(output));
  return output;
#else
  if (idx > 4092) {
    br_aesctr_drbg_generate(&ctx, buffer, 4096);
    idx = 0;
  }
  uint32_t output = *(uint32_t *)(buffer + idx);
  idx += 4;
  return output;
#endif
}
uint8_t RandGen::rand1() {
#if TCS_NUM > 1
  uint8_t output;
  sgx_read_rand((uint8_t *)&output, sizeof(output));
  return output & 1;
#else
  if (idx > 4095) {
    br_aesctr_drbg_generate(&ctx, buffer, 4096);
    idx = 0;
  }
  uint8_t output = *(uint8_t *)(buffer + idx);
  idx += 1;
  return output & 1;
#endif
}

void read_rand(uint8_t *output, size_t size) { sgx_read_rand(output, size); }

// Inside SGX enclave
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