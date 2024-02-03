#define TRUSTED_ENV 1
#include "../Enclave.h"
#include "Enclave_t.h"
////
#include <functional>
#include <unordered_map>

#include "../../common.hpp"
#include "crypt.hpp"
#include "oram/omap.hpp"

using namespace ORAM;
EM::Backend::MemServerBackend* EM::Backend::g_DefaultBackend = nullptr;
OMap<key_type, val_type> omap;
sgx_ec256_private_t private_key;

void ecall_gen_key_pair(uint8_t pubkey[64]) {
  sgx_ec256_public_t public_key;
  sgx_status_t status = generate_key_pair(&private_key, &public_key);
  if (status != SGX_SUCCESS) {
    abort();
  }
  memcpy(pubkey, &public_key, 64);
}

void ecall_omap_init(uint64_t N, uint64_t initSize) {
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
  size_t BackendSize = N * sizeof(pair_type) * 50;
  EM::Backend::g_DefaultBackend =
      new EM::Backend::MemServerBackend(BackendSize);
  omap.SetSize(N);
  std::unique_ptr<uint8_t> kv_buffer(new uint8_t[kv_read_batch_bytes]);
  uint64_t bufferBytes = 0;
  uint64_t bufferOffset = 0;
  try {
    EM::VirtualVector::VirtualReader<std::pair<key_type, val_type>> reader(
        initSize, [&](uint64_t i) {
          if (bufferOffset >= bufferBytes) {
            ocall_Fetch_Next_KV_Batch(&bufferBytes, kv_buffer.get(),
                                      kv_read_batch_bytes);
            bufferOffset = 0;
          }
          std::pair<key_type, val_type> res;
          memcpy(&res.first, kv_buffer.get() + bufferOffset, sizeof(key_type));
          bufferOffset += sizeof(key_type);
          memcpy(&res.second, kv_buffer.get() + bufferOffset, sizeof(val_type));
          bufferOffset += sizeof(val_type);
          return res;
        });
    omap.InitFromReader(reader);
  } catch (std::exception& e) {
    printf("exception: %s\n", e.what());
    abort();
  }
  return;
}

int ecall_omap_find(uint8_t* key, uint8_t* val, uint32_t keyLength,
                    uint32_t valLength) {
  const key_type& k = *reinterpret_cast<key_type*>(key);
  const val_type& v = *reinterpret_cast<val_type*>(val);
  bool res = omap.find(k, *reinterpret_cast<val_type*>(val));
  return (int)res;
}

int ecall_omap_insert(uint8_t* key, uint8_t* val, uint32_t keyLength,
                      uint32_t valLength) {
  const key_type& k = *reinterpret_cast<key_type*>(key);
  const val_type& v = *reinterpret_cast<val_type*>(val);
  bool res = omap.insert(k, v);
  return (int)res;
}

int ecall_omap_delete(uint8_t* key, uint32_t keyLength) {
  const key_type& k = *reinterpret_cast<key_type*>(key);
  bool res = true;  // omap.erase(k);
  return (int)res;
}

int ecall_omap_update(uint8_t* key, uint8_t* val, uint32_t keyLength,
                      uint32_t valLength) {
  const key_type& k = *reinterpret_cast<key_type*>(key);
  const val_type& v = *reinterpret_cast<val_type*>(val);
  bool res = omap.update(k, v);
  return (int)res;
}

void ecall_handle_encrypted_query(uint8_t* encryptedQuery,
                                  uint8_t* encryptedResponse,
                                  uint32_t encryptedQueryLength,
                                  uint32_t encryptedResponseLength) {
  if (encryptedQueryLength != sizeof(EncryptedQuery)) {
    printf("wrong query size\n");
    return;
  }
  if (encryptedResponseLength != sizeof(EncryptedResponse)) {
    printf("wrong response size\n");
    return;
  }
  sgx_ec256_dh_shared_t shared_key;
  EncryptedQuery* encQuery = reinterpret_cast<EncryptedQuery*>(encryptedQuery);
  // sgx_ec256_public_t a_public_key ;
  sgx_status_t status = compute_shared_key(
      &private_key,
      reinterpret_cast<sgx_ec256_public_t*>(&encQuery->senderPubKey),
      &shared_key);
  uint8_t* shared_key_ptr = reinterpret_cast<uint8_t*>(&shared_key);

  // decrypt query using shared key
  Query query;
  encQuery->encQuery.Decrypt(query, shared_key_ptr, encQuery->iv);

  val_type v;
  Response response;
  EncryptedResponse encResponse;
  query.addr.ntoh();
  response.success = omap.find(query.addr, response.balance);

  sgx_read_rand(encResponse.iv, IV_SIZE);
  encResponse.encResponse.Encrypt(response, shared_key_ptr, encResponse.iv);
  memcpy(encryptedResponse, &encResponse, sizeof(EncryptedResponse));
}