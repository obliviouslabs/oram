#include <stdio.h>

#include <string>
#include <thread>
#include <vector>
#define UNTRUSTED_ENV 1
#include "../../common.hpp"
#include "../App.h"
#include "Enclave_u.h"
#ifdef DISK_IO
#include "external_memory/server/enclaveFileServer_untrusted.hpp"
#else
#include "external_memory/server/enclaveMemServer_untrusted.hpp"
#endif
#include <boost/beast/core/detail/base64.hpp>

// #include "async_https_server.hpp"
#define CPPHTTPLIB_OPENSSL_SUPPORT
#include "../../httplib.h"
#include "db_updater.hpp"
#include "kv_db.hpp"
#include "sgx_tcrypto.h"
// since most balances are small, it is more space efficient to directly store
// them as string
using DB_ = KV_DB<std::string, std::string>;

DB_* db;
typename DB_::Iterator dbIt;
std::string publicKeyBase64;
std::counting_semaphore sem(1);  // may set to tcs_max_num

uint64_t ocall_Fetch_Next_KV_Batch(uint8_t* data, uint64_t batchBytes) {
  uint64_t bytes = 0;
  for (; dbIt.isValid() &&
         bytes < batchBytes - sizeof(key_type) - sizeof(val_type);
       dbIt.next()) {
    key_type key;
    val_type val;
    std::string keyStr = dbIt.key();
    if (keyStr.substr(0, 2) != "0x" || keyStr.size() != 42) {
      continue;
    }
    // pad to 42 bytes
    key.set(keyStr);
    val.set(dbIt.value());

    memcpy(data + bytes, &key, sizeof(key_type));
    bytes += sizeof(key_type);
    memcpy(data + bytes, &val, sizeof(val_type));
    bytes += sizeof(val_type);
  }
  return bytes;
}

void updateDBAndORAM(const std::string& logPath,
                     std::unordered_map<std::string, uint256_t>& unstableDiff) {
  std::vector<std::pair<std::string, std::string>> insertList;
  std::vector<std::pair<std::string, std::string>> updateList;
  std::vector<std::string> deleteList;
  uint64_t lastBlock = updateDBFromLog(db, logPath, unstableDiff, insertList,
                                       updateList, deleteList);
  printf("insert %lu keys, update %lu keys, delete %lu keys\n",
         insertList.size(), updateList.size(), deleteList.size());
  SemaphoreLock lock(sem);
  for (auto& kv : updateList) {
    key_type key;
    val_type val;
    key.set(kv.first);
    val.set(kv.second);
    int findFlag = false;
    int ret =
        ecall_omap_update(global_eid, &findFlag, (uint8_t*)&key, (uint8_t*)&val,
                          sizeof(key_type), sizeof(val_type));
    if (ret != SGX_SUCCESS) abort();
    if (findFlag == 0) {
      std::cerr << "key not found when update" << std::endl;
    }
  }

  for (auto& kv : insertList) {
    key_type key;
    val_type val;
    key.set(kv.first);
    val.set(kv.second);
    int findFlag = false;
    int ret =
        ecall_omap_insert(global_eid, &findFlag, (uint8_t*)&key, (uint8_t*)&val,
                          sizeof(key_type), sizeof(val_type));
    if (ret != SGX_SUCCESS) abort();
    if (findFlag == 1) {
      std::cerr << "key already exists when insert" << std::endl;
    }
  }

  for (auto& key : deleteList) {
    key_type k;
    k.set(key);
    int findFlag = false;
    int ret = ecall_omap_delete(global_eid, &findFlag, (uint8_t*)&k,
                                sizeof(key_type));
    if (ret != SGX_SUCCESS) abort();
    if (findFlag == 0) {
      std::cerr << "key not found when delete" << std::endl;
    }
  }
  int ret = ecall_set_last_block(global_eid, lastBlock);
  if (ret != SGX_SUCCESS) abort();
}

void periodicUpdateFromLog() {
  std::string logPath = "tx.log";
  std::unordered_map<std::string, uint256_t> unstableDiff;
  while (true) {
    updateDBAndORAM(logPath, unstableDiff);
    std::this_thread::sleep_for(std::chrono::seconds(5));
  }
}

void handleEncryptedQuery(uint8_t* encryptedQueryPtr,
                          uint8_t* encryptedResponsePtr) {
  uint32_t encryptedQueryLength = sizeof(EncryptedQuery);
  uint32_t encryptedResponseLength = sizeof(EncryptedResponse);
  SemaphoreLock lock(sem);
  sgx_status_t ret = ecall_handle_encrypted_query(
      global_eid, encryptedQueryPtr, encryptedResponsePtr, encryptedQueryLength,
      encryptedResponseLength);
  if (ret != SGX_SUCCESS) abort();
  // send encryptedResponse to client
}

void handleEncryptedQuery(EncryptedQuery& encryptedQuery,
                          EncryptedResponse& encryptedResponse) {
  uint8_t* encryptedQueryPtr = (uint8_t*)&encryptedQuery;
  uint8_t* encryptedResponsePtr = (uint8_t*)&encryptedResponse;
  handleEncryptedQuery(encryptedQueryPtr, encryptedResponsePtr);
}

void handleEncryptedQuery(const std::string& encryptedQueryBase64,
                          std::string& encryptedResponseBase64) {
  std::string encryptedQuery = base64_decode(encryptedQueryBase64);
  if (encryptedQuery.size() != sizeof(EncryptedQuery)) {
    throw std::runtime_error("Invalid encrypted query size");
  }
  uint8_t* encryptedQueryPtr = (uint8_t*)encryptedQuery.data();
  EncryptedResponse encryptedResponse;
  uint8_t* encryptedResponsePtr = (uint8_t*)&encryptedResponse;
  handleEncryptedQuery(encryptedQueryPtr, encryptedResponsePtr);
  encryptedResponseBase64 =
      base64_encode(encryptedResponsePtr, sizeof(EncryptedResponse));
}

std::string getServerPublicKeyBase64() { return publicKeyBase64; }

void InitKeys(DB_* db) {
  sgx_status_t ret = SGX_ERROR_UNEXPECTED;
  std::string sealedPrivateKeyBase64;
  int success = 0;
  if (db->get("public_key", publicKeyBase64) &&
      db->get("sealed_private_key", sealedPrivateKeyBase64)) {
    printf("Recovered public key and sealed private key from db\n");
    printf("public key: %s\n", publicKeyBase64.c_str());
    printf("sealed private key: %s\n", sealedPrivateKeyBase64.c_str());
    std::string sealedPrivateKey = base64_decode(sealedPrivateKeyBase64);
    ret = ecall_set_private_key(global_eid, &success,
                                (uint8_t*)sealedPrivateKey.data(),
                                sealedPrivateKey.size());
    if (ret != SGX_SUCCESS) {
      printf("ecall_set_private_key ecall failed\n");
      abort();
    }
  }
  if (!success) {
    ec256_public_t public_key;
    uint32_t sealedDataSize = 0;
    uint8_t sealedData[1024];  // ecall will fail it detected overflow
    ret = ecall_gen_key_pair(global_eid, &sealedDataSize, (uint8_t*)&public_key,
                             sealedData);

    if (ret != SGX_SUCCESS) {
      printf("ecall_gen_key_pair failed\n");
      abort();
    }
    publicKeyBase64 =
        base64_encode((uint8_t*)&public_key, sizeof(ec256_public_t));
    sealedPrivateKeyBase64 = base64_encode(sealedData, sealedDataSize);
    db->put("public_key", publicKeyBase64);
    db->put("sealed_private_key", sealedPrivateKeyBase64);
  }
}

void InitDB(const char* dbPath) {
  try {
    db = new DB_(dbPath);
    dbIt = db->getIterator();
    dbIt.seekToFirst();
    if (!dbIt.isValid()) {
      printf("db is empty\n");
      return;
    }
  } catch (const std::runtime_error& e) {
    dbIt.reset();
    if (db) {
      delete db;
      db = NULL;
    }
    std::cerr << e.what() << std::endl;
    return;
  }
}

void OutputDB(std::string filePath) {
  std::ofstream ofs(filePath);
  DBMetaData metaData = readMetaData(db);
  ofs << metaData.lastBlock << " " << metaData.lastTxIdx << " "
      << metaData.recordCount << std::endl;
  if (!ofs.is_open()) {
    std::cerr << "Failed to open file " << filePath << std::endl;
    return;
  }
  auto it = db->getIterator();
  for (it.seekToFirst(); it.isValid(); it.next()) {
    if (it.key().substr(0, 2) != "0x") {
      continue;
    }
    ofs << it.key() << " " << it.value() << std::endl;
  }
  ofs.close();
}

void DeleteDB(void) {
  if (db) {
    delete db;
    db = NULL;
  }
}

std::string GetPublicKeyBase64(void) { return publicKeyBase64; }

void ActualMain(const char* dbPath) {
  sgx_status_t ret = SGX_ERROR_UNEXPECTED;
  InitDB(dbPath);
  // OutputDB("usdt_latest.txt");
  // DeleteDB();
  // return;
  InitKeys(db);

  DBMetaData metaData;
  try {
    metaData = readMetaData(db);
  } catch (const std::runtime_error& e) {
    std::cerr << e.what() << std::endl;
    return;
  }
  printf("lastBlock = %lu, recordCount = %lu\n", metaData.lastBlock,
         metaData.recordCount);
  size_t mapSize = (size_t)(metaData.recordCount * 1.2);
  ret = ecall_omap_init(global_eid, mapSize, metaData.recordCount);
  if (ret != SGX_SUCCESS) abort();
  ret = ecall_set_last_block(global_eid, metaData.lastBlock);
  if (ret != SGX_SUCCESS) abort();
  dbIt.reset();
  // auto it = db->getIterator();
  // for (it.seekToFirst(); it.isValid(); it.next()) {
  //   key_type key;
  //   val_type val;
  //   std::string keyStr = it.key();
  //   if (keyStr.substr(0, 2) != "0x" || keyStr.size() != 42) {
  //     std::cout << "key not addr: " << keyStr << std::endl;
  //     continue;
  //   }
  //   key.set(keyStr);
  //   int findFlag = false;
  //   int ret =
  //       ecall_omap_find(global_eid, &findFlag, (uint8_t*)&key,
  //       (uint8_t*)&val,
  //                       sizeof(key_type), sizeof(val_type));
  //   if (ret != SGX_SUCCESS) abort();
  //   if (findFlag == 0) {
  //     std::cerr << "key not found when find " << keyStr << std::endl;
  //   }
  // }

  std::thread updateThread(periodicUpdateFromLog);
  httplib::Server svr;
  svr.Get("/public_key",
          [](const httplib::Request& /*req*/, httplib::Response& res) {
            res.set_content(getServerPublicKeyBase64(), "text/plain");
          });
  svr.Post("/secure", [](const httplib::Request& req, httplib::Response& res) {
    printf("Received encrypted request\n");
    std::string encryptedResponseBase64;
    try {
      handleEncryptedQuery(req.body, encryptedResponseBase64);
      res.set_content(encryptedResponseBase64, "text/plain");
    } catch (const std::exception& e) {
      res.status = 500;
      res.set_content("Error processing the encrypted request", "text/plain");
    }
  });
  printf("Server listening on port 8080\n");
  bool success = svr.listen("localhost", 8080);
  printf("Server listening ends.\n");
  updateThread.join();

  DeleteDB();

  // printf("Did not abort\n");
}
