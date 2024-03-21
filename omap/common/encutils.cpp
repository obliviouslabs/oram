#include <utility>
#ifdef ENCLAVE_MODE
// #include "../Enclave.h"
#include "Enclave_t.h"
#endif

#include <x86intrin.h>

#include "bearssl_aead.h"
#include "bearssl_hash.h"
#include "common/encutils.hpp"
#define memset_s(s, smax, c, n) memset(s, c, n);

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

#ifndef ENCLAVE_MODE_ENCLAVE
void handleErrors(void) {
  throw std::runtime_error("BearSSL error");
}
#endif

void __attribute__((noinline)) sgxsd_br_clear_stack() {
  uint8_t stack[4096];
  memset_s(&stack, sizeof(stack), 0, sizeof(stack));
  _mm256_zeroall();
}


uint64_t secure_hash_with_salt(const uint8_t* data, size_t data_size, const uint8_t (&salt)[16]) {
  uint64_t res;
    // Initialize the hash context
    br_sha256_context ctx;
    br_sha256_init(&ctx);

    // Hash the salt
    br_sha256_update(&ctx, salt, sizeof(salt));

    // Hash the data
    br_sha256_update(&ctx, &data, data_size);

    // Finalize the hash and get the result
    unsigned char hash[br_sha256_SIZE]; // br_sha256_SIZE is normally 32 for SHA-256
    br_sha256_out(&ctx, hash);

    // Copy the result to the output buffer, up to resSize bytes
    memcpy(&res, hash, 8);

    return res;
}

bool sgxsd_aes_gcm_run(bool encrypt,
                               const uint8_t p_key[SGXSD_AES_GCM_KEY_SIZE],
                               const void *p_src, uint32_t src_len, void *p_dst,
                               const uint8_t p_iv[SGXSD_AES_GCM_KEY_SIZE],
                               const void *p_aad, uint32_t aad_len,
                               uint8_t p_mac[SGXSD_AES_GCM_KEY_SIZE]) {
  if (p_key == NULL || ((p_src == NULL || p_dst == NULL) && src_len != 0) ||
      p_iv == NULL || (p_aad == NULL && aad_len != 0) || p_mac == NULL) {
    return 0;
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
    return 1;
  } else {
    if (p_dst != NULL) {
      memset_s(p_dst, src_len, 0, src_len);
    }
    return 0;
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
  return sgxsd_aes_gcm_run(false, key, ciphertext, ciphertextSize,
                                plaintext, iv, nullptr, 0, tag);
}

bool sgxsd_aes_ctr_run(bool encrypt,
                               const uint8_t p_key[SGXSD_AES_GCM_KEY_SIZE],
                               const void *p_src, uint32_t src_len, void *p_dst,
                               const uint8_t p_iv[SGXSD_AES_GCM_KEY_SIZE]) {
  if (p_key == NULL || ((p_src == NULL || p_dst == NULL) && src_len != 0) ||
      p_iv == NULL) {
    return 0;
  }
  br_aes_x86ni_ctr_keys aes_ctx;
  br_aes_x86ni_ctr_init(&aes_ctx, p_key, SGXSD_AES_GCM_KEY_SIZE);
  if (src_len != 0) {
    memmove(p_dst, p_src, src_len);
    br_aes_x86ni_ctr_run(&aes_ctx, p_iv, 0, p_dst, src_len);
  }

  // sgxsd_br_clear_stack();
  // memset_s(&aes_ctx, sizeof(aes_ctx), 0, sizeof(aes_ctx));
  return 1;
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
  return sgxsd_aes_ctr_run(false, key, ciphertext, ciphertextSize,
                                plaintext, iv);
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
        return; // Failure
    }

    size_t read = fread(output, 1, size, fp);
    fclose(fp);

    if (read != size) {
        perror("Failed to read enough bytes");
        // Handle the error, not enough data was read
        return; // Failure
    }

    return; // Success
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
