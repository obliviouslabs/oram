#pragma once
#include <cstdlib>
#include <string>
#include <cstring>
#ifndef ENCLAVE_MODE
#include <boost/multiprecision/cpp_int.hpp>
#endif
#include "external_memory/virtualvector.hpp"
#include "common/mov_intrinsics.hpp"

struct ETH_Addr {
    // 20 bytes

    struct {
        uint64_t addr0_8;
        uint64_t addr9_16;
        uint32_t addr17_20;
    } addr;
    

    // define less operator
    bool operator<(const ETH_Addr& other) const {
        bool is0_8_less = addr.addr0_8 < other.addr.addr0_8;
        bool is0_8_equal = addr.addr0_8 == other.addr.addr0_8;
        bool is9_16_less = addr.addr9_16 < other.addr.addr9_16;
        bool is9_16_equal = addr.addr9_16 == other.addr.addr9_16;
        bool is17_20_less = addr.addr17_20 < other.addr.addr17_20;

        return is0_8_less | (is0_8_equal & (is9_16_less | (is9_16_equal & is17_20_less)));
    }

    bool operator==(const ETH_Addr& other) const {
        return (addr.addr0_8 == other.addr.addr0_8) & (addr.addr9_16 == other.addr.addr9_16) &
               (addr.addr17_20 == other.addr.addr17_20);
    }

    bool operator!=(const ETH_Addr& other) const {
        return (addr.addr0_8 != other.addr.addr0_8) | (addr.addr9_16 != other.addr.addr9_16) |
               (addr.addr17_20 != other.addr.addr17_20);
    }

    template <typename T>
    static T parseHex(const char* hexStr) {
        // hexStr is required to be padded to 42 bytes long
        // 0x + 40 bytes
        T res = 0;
        static constexpr int numChar = sizeof(T) * 2;
        Assert(strlen(hexStr) == numChar);
        for (int i = 0; i < numChar; ++i) {
            T hexDigit = hexStr[i] - '0';
            obliMove(hexDigit > '9', hexDigit, hexDigit - 'A');
            obliMove(hexDigit > 'F', hexDigit, hexDigit - 'a');
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
        addr.addr0_8 = parseHex<uint64_t>(hexStr.c_str() + 2);
        addr.addr9_16 = parseHex<uint64_t>(hexStr.c_str() + 18);
        addr.addr17_20 = parseHex<uint32_t>(hexStr.c_str() + 34);
    }

    #ifdef UNTRUSTED_ENV
    // define friend << operator
    friend std::ostream& operator<<(std::ostream& os, const ETH_Addr& addr) {
        os << "0x";
        os << std::hex << std::setfill('0') << std::setw(16) << addr.addr.addr0_8;
        os << std::hex << std::setfill('0') << std::setw(16) << addr.addr.addr9_16;
        os << std::hex << std::setfill('0') << std::setw(8) << addr.addr.addr17_20;
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
    uint64_t parts[4];

    #ifdef UNTRUSTED_ENV
    void set(const std::string& decStr) {
        boost::multiprecision::uint256_t bigInt(decStr);

        // Mask for extracting 64 bits
        boost::multiprecision::uint256_t mask = (boost::multiprecision::uint256_t(1) << 64) - 1;

        // Extracting 4 parts of 64 bits each
        for (int i = 0; i < 4; ++i) {
            parts[i] = static_cast<uint64_t>((bigInt >> (64 * i)) & mask);
        }
    }
    

    friend std::ostream& operator<<(std::ostream& os, const USDT_Balance& bal) {
        boost::multiprecision::uint256_t bigInt = 0;
        for (int i = 0; i < 4; ++i) {
            bigInt += static_cast<boost::multiprecision::uint256_t>(bal.parts[i]) << (64 * i);
        }
        os << bigInt;
        return os;
    }

    #endif
};

typedef ETH_Addr key_type;
typedef USDT_Balance val_type;
typedef std::pair<key_type, val_type> pair_type;

constexpr uint64_t kv_read_batch_bytes = 1UL << 20; // read 1 MB record a time