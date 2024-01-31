#pragma once
#include <cstdlib>
#include <cstring>
#include <string>
#ifndef ENCLAVE_MODE
#include <boost/multiprecision/cpp_int.hpp>
#endif
#include "common/mov_intrinsics.hpp"
#include "external_memory/virtualvector.hpp"

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
    Assert(strlen(hexStr) == numChar);
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
      throw std::runtime_error("invalid hex string");
    }
    for (int i = 0; i < 5; ++i) {
      part[i] = parseHex<uint32_t>(hexStr.c_str() + 2 + i * 8);
    }
  }

#ifdef UNTRUSTED_ENV
  // define friend << operator
  friend std::ostream& operator<<(std::ostream& os, const ETH_Addr& addr) {
    os << "0x";
    for (int i = 0; i < 5; ++i) {
      // pad to 8 bytes
      os << std::hex << std::setfill('0') << std::setw(8) << addr.part[i];
    }
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

struct USDT_Balance {
  uint64_t part[4];

#ifdef UNTRUSTED_ENV
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

  friend std::ostream& operator<<(std::ostream& os, const USDT_Balance& bal) {
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

typedef ETH_Addr key_type;
typedef USDT_Balance val_type;
typedef std::pair<key_type, val_type> pair_type;

constexpr uint64_t kv_read_batch_bytes = 1UL << 20;  // read 1 MB record a time