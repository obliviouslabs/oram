#pragma once
#include <inttypes.h>

#include <utility>

#include "common/encutils.hpp"
#include "common/tracing/tracer.hpp"
#include "common/utils.hpp"

namespace Concepts {
template <typename T>
concept Encryptable = requires(typename T::Encrypted_t et) {
  { et };
};
}  // namespace Concepts

template <typename T>
T prp(T val);

template <typename T>
// requires (IS_POD<T>())
struct Encrypted {
  // static_assert(IS_POD<T>());

  static constexpr uint64_t SIZE = sizeof(T);
  uint8_t data[SIZE];  // We don't need to adjust this because CTR modes don't
                       // need padding.
  uint8_t iv[IV_SIZE];

  // Encrypted& operator=(const Encrypted&) = delete;
  // Encrypted(const Encrypted&) = delete;

  INLINE void Encrypt(const T& in) {
    // PROFILE_F();
    GetRandIV(iv);
    aes_256_ctr_encrypt(
        SIZE, const_cast<uint8_t*>(reinterpret_cast<const uint8_t*>(&in)), iv,
        data);
  }

  INLINE void Decrypt(T& out) /*const*/ {
    // PROFILE_F();
    bool r =
        aes_256_ctr_decrypt(SIZE, data, iv, reinterpret_cast<uint8_t*>(&out));
    Assert(r);
    IGNORE_UNUSED(r);
  }

  INLINE void Encrypt(const T& in, const uint8_t key[AES_BLOCK_SIZE]) {
    // PROFILE_F();
    GetRandIV(iv);
    aex_ctx_t ctx;
    br_aes_x86ni_ctr_init(&ctx, key, AES_BLOCK_SIZE);
    aes_256_ctr_encrypt(
        SIZE, const_cast<uint8_t*>(reinterpret_cast<const uint8_t*>(&in)), iv,
        data, &ctx);
  }

  INLINE void Decrypt(T& out, const uint8_t key[AES_BLOCK_SIZE]) {
    aex_ctx_t ctx;
    br_aes_x86ni_ctr_init(&ctx, key, AES_BLOCK_SIZE);
    bool r = aes_256_ctr_decrypt(SIZE, data, iv,
                                 reinterpret_cast<uint8_t*>(&out), &ctx);
    Assert(r);
    IGNORE_UNUSED(r);
  }

#ifndef ENCLAVE_MODE
  friend std::ostream& operator<<(std::ostream& o, const Encrypted& x) {
    T v;
    const_cast<Encrypted&>(x).Decrypt(v);
    o << "E{" << v << "}";
    return o;
  }
#endif

  bool operator==(Encrypted& o) {
    T d1;
    Decrypt(d1);
    T d2;
    o.Decrypt(d2);
    return d1 == d2;
  }
};
static_assert(IS_POD<Encrypted<int>>());

template <typename T>
// requires (IS_POD<T>())
struct FreshEncrypted {
  // static_assert(IS_POD<T>());

  static constexpr uint64_t SIZE = sizeof(T);
  uint8_t data[SIZE];  // We don't need to adjust this because CTR modes don't
                       // need padding.
  uint8_t tag[MAC_SIZE];

  // Encrypted& operator=(const Encrypted&) = delete;
  // Encrypted(const Encrypted&) = delete;

  INLINE void Encrypt(const T& in, const uint8_t iv[IV_SIZE]) {
    // PROFILE_F();

    aes_256_gcm_encrypt(
        SIZE, const_cast<uint8_t*>(reinterpret_cast<const uint8_t*>(&in)), iv,
        tag, data);
  }

  INLINE void Decrypt(T& out, const uint8_t iv[IV_SIZE]) /*const*/ {
    bool r = aes_256_gcm_decrypt(SIZE, data, iv, tag,
                                 reinterpret_cast<uint8_t*>(&out));
    if (!r) {
      throw std::runtime_error("Authentication failed during decrypt.");
    }
  }

  INLINE void Encrypt(const T& in, const uint8_t iv[IV_SIZE],
                      const uint8_t key[AES_BLOCK_SIZE]) {
    // PROFILE_F();
    aex_ctx_t ctx;
    br_aes_x86ni_ctr_init(&ctx, key, AES_BLOCK_SIZE);
    aes_256_gcm_encrypt(
        SIZE, const_cast<uint8_t*>(reinterpret_cast<const uint8_t*>(&in)), iv,
        tag, data, &ctx);
  }

  INLINE void Decrypt(T& out, const uint8_t iv[IV_SIZE],
                      const uint8_t key[AES_BLOCK_SIZE]) {
    aex_ctx_t ctx;
    br_aes_x86ni_ctr_init(&ctx, key, AES_BLOCK_SIZE);
    bool r = aes_256_gcm_decrypt(SIZE, data, iv, tag,
                                 reinterpret_cast<uint8_t*>(&out), &ctx);
    if (!r) {
      throw std::runtime_error("Authentication failed during decrypt.");
    }
  }

#ifndef ENCLAVE_MODE
  friend std::ostream& operator<<(std::ostream& o, const FreshEncrypted& x) {
    o << "E{" << typeid(T).name() << "}";
    return o;
  }
#endif

  // not available as we don't have the iv.
  //
  bool operator==(FreshEncrypted& o) = delete;
};
static_assert(IS_POD<FreshEncrypted<int>>());

// Structure with same interfaces as Encrypted, but with
// everything as plaintext for debugging.
//
template <typename T>
struct NonEncrypted {
  static_assert(IS_POD<T>());

  static constexpr uint64_t SIZE = sizeof(T);
  T data;

  INLINE void Encrypt(const T& in) { data = in; }

  INLINE void Decrypt(T& out) /*const*/ { out = data; }

  INLINE void Encrypt(const T& in, const uint8_t key[AES_BLOCK_SIZE]) {
    data = in;
  }

  INLINE void Decrypt(T& out, const uint8_t key[AES_BLOCK_SIZE]) { out = data; }

#ifndef ENCLAVE_MODE
  friend std::ostream& operator<<(std::ostream& o, const NonEncrypted& x) {
    T v;
    const_cast<NonEncrypted&>(x).Decrypt(v);
    o << "E{" << v << "}";
    return o;
  }
#endif

  bool operator==(const NonEncrypted& o) const { return data == o.data; }
};
static_assert(IS_POD<NonEncrypted<int>>());

template <typename PublicData, typename PrivateData>
  requires(IS_POD<PublicData>()) && (IS_POD<PrivateData>()) &&
          ::Concepts::Encryptable<PrivateData>
struct MixedEncryptable {
  using PublicData_t = PublicData;
  using PrivateData_t = PrivateData;
  PublicData pub;
  PrivateData priv;

#ifndef ENCLAVE_MODE
  friend std::ostream& operator<<(std::ostream& o, const MixedEncryptable& x) {
    o << "ME{" << x.pub << ", " << x.priv << "}";
    return o;
  }
#endif

  bool operator==(const MixedEncryptable& o) const {
    return (pub == o.pub) * (priv == o.priv);
  }

  // Classes that extend this need to declare Encrypted_t.
  //
};