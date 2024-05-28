#define TRUSTED_ENV 1
#include "../Enclave.h"
#include "Enclave_t.h"
////
#include <functional>
#include <unordered_map>

#include "../../common.hpp"
#include "crypt.hpp"
#include "odsl/omap.hpp"
#include "sgx_error.h"
#include "sgx_report.h"
#include "sgx_spinlock.h"
#include "sgx_trts.h"
#include "sgx_tseal.h"
#include "sgx_utils.h"

uint32_t enclave_create_report(const sgx_target_info_t* p_qe3_target,
                               const sgx_report_data_t* p_data,
                               sgx_report_t* p_report) {
  return sgx_create_report(p_qe3_target, p_data, p_report);
}

using namespace ODSL;
EM::Backend::MemServerBackend* EM::Backend::g_DefaultBackend = nullptr;
OMap<key_type, val_type> omap;
// lock for the global OMAP
// TODO: use fine grained lock within omap
sgx_spinlock_t omap_lock = SGX_SPINLOCK_INITIALIZER;
sgx_ec256_private_t private_key;
uint64_t globalLastBlock = 0;

class OMapCritical {
 public:
  OMapCritical() { sgx_spin_lock(&omap_lock); }
  ~OMapCritical() { sgx_spin_unlock(&omap_lock); }
};

// return sealed private key size
uint32_t ecall_gen_key_pair(uint8_t pubkey[64], uint8_t sealedPrivKey[1024]) {
  static_assert(sizeof(ec256_public_t) == 64);
  sgx_ec256_public_t public_key;
  sgx_status_t status = generate_key_pair(&private_key, &public_key);
  ec256_public_t public_key_big_endian = convert_to_ec256_public_t(public_key);
  if (status != SGX_SUCCESS) {
    abort();
  }
  memcpy(pubkey, &public_key_big_endian, 64);
  size_t sealed_data_size;
  uint8_t* sealed_data;
  status = seal_private_key(&private_key, &sealed_data, &sealed_data_size);
  if (status != SGX_SUCCESS) {
    printf("seal private key failed\n");
    abort();
  }
  if (sealed_data_size > 1024) {
    printf("sealed data size too large\n");
    abort();
  }
  memcpy(sealedPrivKey, sealed_data, sealed_data_size);
  free(sealed_data);
  return sealed_data_size;
}

int ecall_set_private_key(uint8_t* sealedPrivKey, uint32_t sealedPrivKeySize) {
  sgx_status_t status =
      unseal_private_key(sealedPrivKey, sealedPrivKeySize, &private_key);
  if (status != SGX_SUCCESS) {
    printf("unseal private key failed, generating new key\n");
    return 0;
  }
  return 1;
}

void ecall_omap_init(uint64_t N, uint64_t initSize) {
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
  size_t BackendSize = N * sizeof(pair_type) * 50;
  EM::Backend::g_DefaultBackend =
      new EM::Backend::MemServerBackend(BackendSize);
  OMapCritical section;
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
  OMapCritical section;
  bool res = omap.Find(k, *reinterpret_cast<val_type*>(val));
  return (int)res;
}

int ecall_omap_insert(uint8_t* key, uint8_t* val, uint32_t keyLength,
                      uint32_t valLength) {
  const key_type& k = *reinterpret_cast<key_type*>(key);
  const val_type& v = *reinterpret_cast<val_type*>(val);
  OMapCritical section;
  bool res = omap.Insert(k, v);
  return (int)res;
}

int ecall_omap_delete(uint8_t* key, uint32_t keyLength) {
  const key_type& k = *reinterpret_cast<key_type*>(key);
  bool res = omap.Erase(k);
  return (int)res;
}

int ecall_omap_update(uint8_t* key, uint8_t* val, uint32_t keyLength,
                      uint32_t valLength) {
  const key_type& k = *reinterpret_cast<key_type*>(key);
  const val_type& v = *reinterpret_cast<val_type*>(val);
  OMapCritical section;
  bool res = omap.Insert(k, v);
  return (int)res;
}

void ecall_set_last_block(uint64_t lastBlock) { globalLastBlock = lastBlock; }

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
  // to little endianess
  sgx_ec256_public_t a_public_key = convert_to_sgx_ec256_public_t(
      *reinterpret_cast<ec256_public_t*>(&encQuery->senderPubKey));

  sgx_status_t status =
      compute_shared_key(&private_key, &a_public_key, &shared_key);
  if (status != SGX_SUCCESS) {
    printf("compute_shared_key failed\n");
    return;
  }
  uint8_t* shared_key_ptr = reinterpret_cast<uint8_t*>(&shared_key);

  // reverse the shared secret to match openssl endianess
  for (int i = 0; i < 16; ++i) {
    std::swap(shared_key.s[i], shared_key.s[31 - i]);
  }
  // decrypt query using shared key
  Query query;
  encQuery->encQuery.Decrypt(query, encQuery->iv, shared_key_ptr);
  val_type v;
  Response response;
  EncryptedResponse encResponse;
  response.nounce = query.nounce;
  {
    OMapCritical section;
    response.tillBlock = globalLastBlock;  // TODO add pending status
    uint64_t start, end;
    ocall_measure_time(&start);
    bool found = omap.Find(query.addr, response.balance);
    ocall_measure_time(&end);
    response.queryTime = end - start;
    obliMove(!found, response.balance, ERC20_Balance());
    hton((uint8_t*)&response.tillBlock, sizeof(response.tillBlock));
    hton((uint8_t*)&response.queryTime, sizeof(response.queryTime));
    response.success = true;
  }

  sgx_read_rand(encResponse.iv, IV_SIZE);
  encResponse.encResponse.Encrypt(response, encResponse.iv, shared_key_ptr);
  memcpy(encryptedResponse, &encResponse, sizeof(EncryptedResponse));
}