#pragma once
#include <cstdlib>
#include <cstring>
#include <string>
#ifndef ENCLAVE_MODE
#include <boost/multiprecision/cpp_int.hpp>
#endif
#include "common/encrypted.hpp"
#include "common/mov_intrinsics.hpp"

// #include "sgx_tcrypto.h"
#define ECP256_KEY_SIZE 32

// always use little endian when send through network
bool isLittleEndian() {
  const uint32_t one = 1;  // Use uint32_t instead of an integer literal
  return *((const uint8_t*)(&one)) == 1;
}

struct ETH_Addr {
  // 20 bytes

  uint32_t part[5];

  // define less operator
  bool operator<(const ETH_Addr& other) const {
    bool res = false;
    bool eq = true;
    for (int i = 0; i < 5; ++i) {
      res |= (eq & (part[i] < other.part[i]));
      eq &= part[i] == other.part[i];
    }
    return res;
  }

  bool operator==(const ETH_Addr& other) const {
    bool eq = true;
    for (int i = 0; i < 5; ++i) {
      eq &= (part[i] == other.part[i]);
    }
    return eq;
  }

  bool operator!=(const ETH_Addr& other) const { return !(*this == other); }

  template <typename T>
  static T parseHex(const char* hexStr) {
    T res = 0;
    static constexpr int numChar = sizeof(T) * 2;
    for (int i = 0; i < numChar; ++i) {
      T hexDigit = hexStr[i] - '0';
      obliMove(hexStr[i] > '9', hexDigit, (T)(hexStr[i] - ('A' - 10)));
      obliMove(hexStr[i] > 'F', hexDigit, (T)(hexStr[i] - ('a' - 10)));
      res = (res << 4) | hexDigit;
    }
    return res;
  }

  // we require set to be oblivious to avoid leaking the key
  void set(const std::string& hexStr) {
    // hexStr is required to be padded to 42 bytes long
    // 0x + 40 bytes
    if (hexStr.size() != 42 || hexStr[0] != '0' || hexStr[1] != 'x') {
      printf("invalid hex string: %s\n", hexStr.c_str());
      throw std::runtime_error("invalid hex string");
    }
    for (int i = 0; i < 5; ++i) {
      part[i] = parseHex<uint32_t>(hexStr.c_str() + 2 + i * 8);
    }
  }

  void ntoh() {
    if (isLittleEndian()) {
      return;
    }
    for (int i = 0; i < 5; ++i) {
      part[i] = __builtin_bswap32(part[i]);
    }
  }

  void hton() { ntoh(); }

  void print() const {
    printf("0x");
    for (int i = 0; i < 5; ++i) {
      printf("%08x", part[i]);
    }
  }

#ifndef TRUSTED_ENV
  // define friend << operator
  friend std::ostream& operator<<(std::ostream& os, const ETH_Addr& addr) {
    std::ostringstream oss;
    oss << "0x";
    for (int i = 0; i < 5; ++i) {
      // pad to 8 bytes
      oss << std::hex << std::setfill('0') << std::setw(8) << addr.part[i];
    }
    os << oss.str();
    return os;
  }

  friend std::istream& operator>>(std::istream& is, ETH_Addr& addr) {
    std::string str;
    is >> str;
    addr.set(str);
    return is;
  }
#endif
};

struct ERC20_Balance {
  uint64_t part[4] = {0};

  // ntoh
  void ntoh() {
    if (isLittleEndian()) {
      return;
    }
    for (int i = 0; i < 4; ++i) {
      part[i] = __builtin_bswap64(part[i]);
    }
  }

  // hton
  void hton() { ntoh(); }

#ifndef TRUSTED_ENV
  void set(const std::string& decStr) {
    boost::multiprecision::uint256_t bigInt(decStr);

    // Mask for extracting 64 bits
    boost::multiprecision::uint256_t mask =
        (boost::multiprecision::uint256_t(1) << 64) - 1;

    // Extracting 4 part of 64 bits each
    for (int i = 0; i < 4; ++i) {
      part[i] = static_cast<uint64_t>((bigInt >> (64 * i)) & mask);
    }
  }

  friend std::ostream& operator<<(std::ostream& os, const ERC20_Balance& bal) {
    boost::multiprecision::uint256_t bigInt = 0;
    for (int i = 0; i < 4; ++i) {
      bigInt += static_cast<boost::multiprecision::uint256_t>(bal.part[i])
                << (64 * i);
    }
    os << bigInt;
    return os;
  }

#endif
};

enum class QueryType { READ_BALANCE };

enum class CoinType { USDT };

typedef uint64_t Nounce;

struct Query {
  QueryType type;
  CoinType coinType;
  ETH_Addr addr;
  Nounce nounce;
};

struct Response {
  ERC20_Balance balance;
  uint64_t tillBlock;
  Nounce nounce;  // to prevent replay attack
  bool success;
};

struct ec256_public_t {
  uint8_t gx[ECP256_KEY_SIZE];
  uint8_t gy[ECP256_KEY_SIZE];
};

struct EncryptedQuery {
  FreshEncrypted<Query> encQuery;
  ec256_public_t senderPubKey;
  uint8_t iv[IV_SIZE];
};

struct EncryptedResponse {
  FreshEncrypted<Response> encResponse;
  uint8_t iv[IV_SIZE];
};

typedef ETH_Addr key_type;
typedef ERC20_Balance val_type;
typedef std::pair<key_type, val_type> pair_type;

constexpr uint64_t kv_read_batch_bytes = 1UL << 20;  // read 1 MB record a time

#ifndef TRUSTED_ENV
#include <openssl/bio.h>
#include <openssl/buffer.h>
#include <openssl/evp.h>

#include <cstdio>
#include <string>

std::string base64_encode(const unsigned char* buffer, size_t length) {
  BIO *bio, *b64;
  BUF_MEM* bufferPtr;

  b64 = BIO_new(BIO_f_base64());
  bio = BIO_new(BIO_s_mem());
  bio = BIO_push(b64, bio);

  BIO_set_flags(bio, BIO_FLAGS_BASE64_NO_NL);  // Ignore newlines - write
                                               // everything in one line
  BIO_write(bio, buffer, length);
  BIO_flush(bio);
  BIO_get_mem_ptr(bio, &bufferPtr);
  BIO_set_close(bio, BIO_NOCLOSE);

  // Create a new string from the buffer
  std::string encodedData(bufferPtr->data, bufferPtr->length);

  // Free memory
  BIO_free_all(bio);
  BUF_MEM_free(bufferPtr);

  return encodedData;
}

std::string base64_decode(const std::string& encodedData) {
  // Calculate the length of the decoded data
  size_t decodeLen = (encodedData.size() / 4) * 3;
  std::string buffer;
  buffer.resize(
      decodeLen);  // Resize the string to accommodate the decoded data

  BIO *b64, *bio;

  // Set up a BIO chain for base64 decoding
  b64 = BIO_new(BIO_f_base64());
  bio = BIO_new_mem_buf(encodedData.data(), encodedData.length());
  bio = BIO_push(b64, bio);

  // Do not use newlines to flush buffer
  BIO_set_flags(bio, BIO_FLAGS_BASE64_NO_NL);

  // Perform the decoding
  int length = BIO_read(bio, &buffer[0], encodedData.length());
  buffer.resize(length);  // Resize the string based on the actual length read

  // Free memory
  BIO_free_all(bio);

  return buffer;
}
#endif