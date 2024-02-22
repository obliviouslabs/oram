#define CPPHTTPLIB_OPENSSL_SUPPORT
#include <openssl/ec.h>
#include <openssl/ecdh.h>
#include <openssl/pem.h>
#include <openssl/rand.h>
#include <stdio.h>
#include <stdlib.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#include "common.hpp"
#include "httplib.h"

// by Aashish Dugar
/**
 * See https://wiki.openssl.org/index.php/Elliptic_Curve_Diffie_Hellman
 * for source. The function below are based on the ecdh_low algorithm
 * described on that page and utilizes the OpenSSL low-level APIs for
 * Elliptic Curve Diffie Hellman key exchange algorithm.
 */

#include <assert.h>
#include <openssl/ec.h>
#include <openssl/ecdh.h>
#include <openssl/evp.h>
#include <stdio.h>

void testEncrypt() {
  struct Data {
    char data[32];
  };
  Data data;
  for (int i = 0; i < 32; ++i) {
    data.data[i] = 'a';
  }
  FreshEncrypted<Data> encData;
  uint8_t iv[IV_SIZE];
  for (int i = 0; i < IV_SIZE; ++i) {
    iv[i] = '0';
  }
  uint8_t key[32];
  for (int i = 0; i < 32; ++i) {
    key[i] = '1';
  }
  key[16] = '2';
  encData.Encrypt(data, key, iv);
  printf("Test Encrypted Query: ");
  for (int i = 0; i < sizeof(encData); ++i) {
    printf("%02x", ((uint8_t *)&encData)[i]);
  }
  printf("\n");
}

EC_KEY *create_key(void) {
  EC_KEY *key;
  if (NULL == (key = EC_KEY_new_by_curve_name(NID_X9_62_prime256v1))) {
    printf("Failed to create key curve\n");
    return NULL;
  }

  if (1 != EC_KEY_generate_key(key)) {
    printf("Failed to generate key\n");
    return NULL;
  }
  return key;
}

unsigned char *get_secret(EC_KEY *key, const EC_POINT *peer_pub_key,
                          size_t *secret_len) {
  int field_size;
  unsigned char *secret;

  field_size = EC_GROUP_get_degree(EC_KEY_get0_group(key));
  *secret_len = (field_size + 7) / 8;

  if (NULL == (secret = (unsigned char *)OPENSSL_malloc(*secret_len))) {
    printf("Failed to allocate memory for secret");
    return NULL;
  }

  *secret_len = ECDH_compute_key(secret, *secret_len, peer_pub_key, key, NULL);

  if (*secret_len <= 0) {
    OPENSSL_free(secret);
    return NULL;
  }
  return secret;
}

int generate_secure_iv(unsigned char *iv, int iv_length) {
  if (!iv || iv_length <= 0) {
    printf("Invalid IV buffer or length.\n");
    return 0;
  }

  // Generate random bytes for the IV
  if (RAND_bytes(iv, iv_length) != 1) {
    // RAND_bytes returns 1 on success, 0 otherwise.
    printf("Failed to generate secure IV.\n");
    return 0;
  }

  return 1;  // Success
}

std::string makeBalanceQueryBody(CoinType coinType, std::string &addr,
                                 const Nounce &nounce,
                                 const ec256_public_t &client_pub_key,
                                 const uint8_t *shared_secret) {
  Query query;
  EncryptedQuery encQuery;
  if (addr.length() < 2 || addr.substr(0, 2) != "0x") {
    printf("client invalid ethereum address %s\n", addr.c_str());
    throw std::invalid_argument("Invalid ethereum address (client)");
  }
  if (addr.length() < 42) {
    addr.insert(2, 42 - addr.length(), '0');
  }
  query.addr.set(addr);
  query.addr.hton();
  query.coinType = coinType;
  query.type = QueryType::READ_BALANCE;
  query.nounce = nounce;
  int status = generate_secure_iv(encQuery.iv, IV_SIZE);
  if (status != 1) {
    printf("Failed to generate secure IV\n");
    return "";
  }
  for (int i = 0; i < IV_SIZE; ++i) {
    printf("%02x", encQuery.iv[i]);
  }
  printf("\n");
  encQuery.encQuery.Encrypt(query, shared_secret, encQuery.iv);
  encQuery.senderPubKey = client_pub_key;
  return base64_encode((uint8_t *)&encQuery, sizeof(EncryptedQuery));
}

std::string decodeResponse(const std::string &response_base64,
                           const uint64_t &nounce,
                           const uint8_t *shared_secret) {
  std::string response = base64_decode(response_base64);
  if (response.size() != sizeof(EncryptedResponse)) {
    printf("Invalid encrypted response size\n");
    return "";
  }

  EncryptedResponse encResponse;
  memcpy(&encResponse, response.data(), sizeof(EncryptedResponse));
  Response responseObj;
  encResponse.encResponse.Decrypt(responseObj, shared_secret, encResponse.iv);
  if (responseObj.success) {
    if (responseObj.nounce != nounce) {
      return "Error: nounce mismatch, possible replay attack";
    }
    responseObj.balance.ntoh();
    std::ostringstream oss;
    oss << responseObj.balance << " till block " << responseObj.tillBlock;
    return oss.str();
  } else {
    return "Error: balance not found";
  }
}

EC_POINT *pub_key_to_ec_point(const EC_GROUP *group,
                              const struct ec256_public_t *pub_key) {
  if (!group || !pub_key) return NULL;

  // Create a new BN (Big Number) for the x and y coordinates
  BIGNUM *x = BN_new();
  BIGNUM *y = BN_new();
  if (!x || !y) {
    fprintf(stderr, "Failed to allocate BIGNUMs\n");
    if (x) BN_free(x);
    return NULL;
  }

  // Convert the x and y components from the SGX public key to BIGNUMs
  BN_bin2bn(pub_key->gx, sizeof(pub_key->gx), x);
  BN_bin2bn(pub_key->gy, sizeof(pub_key->gy), y);

  // Create a new EC point and set its coordinates
  EC_POINT *point = EC_POINT_new(group);
  if (!point ||
      !EC_POINT_set_affine_coordinates_GFp(group, point, x, y, NULL)) {
    fprintf(stderr, "Failed to create or set the EC_POINT\n");
    if (point) EC_POINT_free(point);
    point = NULL;
  }

  // Free the BIGNUM resources
  BN_free(x);
  BN_free(y);

  return point;
}

int ec_point_to_pub_key(const EC_GROUP *group, const EC_POINT *point,
                        struct ec256_public_t *pub_key) {
  if (!group || !point || !pub_key) return 0;

  BIGNUM *x = BN_new();
  BIGNUM *y = BN_new();
  if (!x || !y) {
    fprintf(stderr, "Failed to allocate BIGNUMs\n");
    if (x) BN_free(x);
    return 0;
  }

  // Use BN_CTX in the conversion process to manage temporary BIGNUMs
  BN_CTX *ctx = BN_CTX_new();
  if (!ctx) {
    fprintf(stderr, "Failed to create BN_CTX\n");
    BN_free(x);
    BN_free(y);
    return 0;
  }

  if (!EC_POINT_get_affine_coordinates_GFp(group, point, x, y, ctx)) {
    fprintf(stderr, "Failed to get affine coordinates\n");
    BN_free(x);
    BN_free(y);
    BN_CTX_free(ctx);
    return 0;
  }

  // Ensure the BIGNUMs fit into the ec256_public_t structure
  if (BN_num_bytes(x) > sizeof(pub_key->gx) ||
      BN_num_bytes(y) > sizeof(pub_key->gy)) {
    fprintf(stderr, "BIGNUMs do not fit into ec256_public_t\n");
    BN_free(x);
    BN_free(y);
    BN_CTX_free(ctx);
    return 0;
  }

  // Convert the BIGNUMs to byte arrays
  BN_bn2bin(x, pub_key->gx + sizeof(pub_key->gx) -
                   BN_num_bytes(x));  // Adjust for big endian
  BN_bn2bin(y, pub_key->gy + sizeof(pub_key->gy) -
                   BN_num_bytes(y));  // Adjust for big endian

  // Cleanup
  BN_free(x);
  BN_free(y);
  BN_CTX_free(ctx);

  return 1;  // Success
}

void interactive(httplib::Client &cli, const ec256_public_t &client_pub_key,
                 const unsigned char *shared_secret) {
  Nounce nounce = 0;
  while (true) {
    std::string addr;
    std::cout << "Enter an ethereum address: " << std::endl;
    std::cin >> addr;

    std::string body = makeBalanceQueryBody(CoinType::USDT, addr, nounce,
                                            client_pub_key, shared_secret);
    auto res = cli.Post("/secure", body, "application/octet-stream");
    if (res && res->status == 200) {
      std::string resValue = decodeResponse(res->body, nounce, shared_secret);
      std::cout << "Response: " << resValue << std::endl;
    } else {
      std::cout << "Error: " << res.error() << std::endl;
    }
    ++nounce;
  }
}

void stressTest(httplib::Client &cli, const ec256_public_t &client_pub_key,
                const unsigned char *shared_secret) {
  // measure time
  Nounce nounce = 0;
  std::chrono::high_resolution_clock::time_point t1 =
      std::chrono::high_resolution_clock::now();
  for (uint64_t round = 0; round < 1000; ++round) {
    std::string addr = "0x9086a22abf1a5f072242669fccfe4536c24c4084";
    std::string body = makeBalanceQueryBody(CoinType::USDT, addr, nounce,
                                            client_pub_key, shared_secret);
    auto res = cli.Post("/secure", body, "application/octet-stream");
    if (res && res->status == 200) {
      std::string resValue = decodeResponse(res->body, nounce, shared_secret);
      std::cout << "Response: " << resValue << std::endl;
    } else {
      std::cout << "Error: " << res.error() << std::endl;
    }
    ++nounce;
  }
  std::chrono::high_resolution_clock::time_point t2 =
      std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
  std::cout << "Time taken: " << duration << " microseconds for " << 1000
            << " queries" << std::endl;
}

int main(int argc, char **argv) {
  // Example usage
  // std::string host = "example.com";
  // testEncrypt();
  httplib::Client cli("localhost", 8080);
  auto res = cli.Get("/public_key");
  if (res && res->status == 200) {
    std::cout << "Response: " << res->body << std::endl;
  } else {
    std::cout << "Error: " << res.error() << std::endl;
    return -1;
  }

  ec256_public_t server_pub_key;
  std::string server_pub_key_str = base64_decode(res->body);
  if (server_pub_key_str.size() != sizeof(ec256_public_t)) {
    printf("Invalid server public key size %lu bytes, expected %lu bytes\n",
           server_pub_key_str.size(), sizeof(ec256_public_t));
    return -1;
  }
  printf("Server public key received\n");
  memcpy(&server_pub_key, server_pub_key_str.data(), sizeof(ec256_public_t));
  EC_GROUP *group = EC_GROUP_new_by_curve_name(NID_X9_62_prime256v1);
  EC_POINT *server_pub_key_ec_point =
      pub_key_to_ec_point(group, &server_pub_key);
  // Here, after receiving the server's public key, perform ECDH, AES
  // encryption, and Base64 encoding as necessary for the actual application
  // logic.
  EC_KEY *client_key = create_key();
  if (!client_key) {
    printf("Failed to generate client ECC key\n");
    return -1;
  }
  printf("Client ECC key generated\n");

  const EC_POINT *client_pub_key_ec_point = EC_KEY_get0_public_key(client_key);
  if (!client_pub_key_ec_point) {
    printf("Failed to serialize client public key ec point\n");
    return -1;
  }
  ec256_public_t client_pub_key;
  if (!ec_point_to_pub_key(group, client_pub_key_ec_point, &client_pub_key)) {
    printf("Failed to serialize client public key\n");
    return -1;
  }

  printf("Client public key serialized\n");
  size_t shared_secret_len;
  unsigned char *shared_secret =
      get_secret(client_key, server_pub_key_ec_point, &shared_secret_len);
  // unsigned char *shared_secret =
  //     compute_shared_secret(client_key, (unsigned char *)&server_pub_key,
  //                           sizeof(ec256_public_t), &shared_secret_len);
  if (!shared_secret) {
    printf("Failed to compute shared secret\n");
    return -1;
  }

  printf("Shared secret computed\n");
  interactive(cli, client_pub_key, shared_secret);
  // stressTest(cli, client_pub_key, shared_secret);
  return 0;
}