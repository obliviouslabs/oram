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
#include "kvDB.hpp"
#include "sgx_tcrypto.h"

// since most balances are small, it is more space efficient to directly store
// them as string
using DB_ = KV_DB<std::string, std::string>;

DB_* db;
typename DB_::Iterator dbIt;
ec256_public_t public_key;

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

void updateDBAndORAM(const std::string& logPath) {
  std::vector<std::pair<std::string, std::string>> insertList;
  std::vector<std::pair<std::string, std::string>> updateList;
  std::vector<std::string> deleteList;
  updateDBFromLog(db, logPath, insertList, updateList, deleteList);
  printf("insert %lu keys, update %lu keys, delete %lu keys\n",
         insertList.size(), updateList.size(), deleteList.size());
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
}

void periodicUpdateFromLog() {
  std::string logPath = "tx.log";
  while (true) {
    updateDBAndORAM(logPath);
    std::this_thread::sleep_for(std::chrono::seconds(5));
  }
}

void handleEncryptedQuery(uint8_t* encryptedQueryPtr,
                          uint8_t* encryptedResponsePtr) {
  uint32_t encryptedQueryLength = sizeof(EncryptedQuery);
  uint32_t encryptedResponseLength = sizeof(EncryptedResponse);
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

std::string getServerPublicKeyBase64() {
  std::string publicKeyStr =
      base64_encode((uint8_t*)&public_key, sizeof(sgx_ec256_public_t));
  return publicKeyStr;
}

// void handleHttpRequest(const httplib::Request& request,
//                        httplib::Response& response) {
//   if (request.target() == "/public_key") {
//     // Client is requesting the server's public key
//     response.result(http::status::ok);
//     response.set(http::field::content_type, "text/plain");
//     response.body() = getServerPublicKeyBase64();
//   } else if (request.target() == "/secure") {
//     // Client is sending an encrypted request
//     try {
//       // Encrypt the response
//       std::string encryptedResponseData;
//       handleEncryptedQuery(request.body(), encryptedResponseData);

//       response.result(http::status::ok);
//       response.set(http::field::content_type, "text/plain");
//       response.body() = encryptedResponseData;
//     } catch (const std::exception& e) {
//       // Handle decryption errors or other exceptions
//       response.result(http::status::internal_server_error);
//       response.set(http::field::content_type, "text/plain");
//       response.body() = "Error processing the encrypted request";
//     }
//   } else {
//     // Unrecognized request path
//     response.result(http::status::not_found);
//     response.set(http::field::content_type, "text/plain");
//     response.body() = "Not Found";
//   }

//   // Prepare the payload for transmission
//   response.prepare_payload();
// }

void ActualMain(void) {
  sgx_status_t ret = SGX_ERROR_UNEXPECTED;

  ret = ecall_gen_key_pair(global_eid, (uint8_t*)&public_key);
  if (ret != SGX_SUCCESS) abort();
  size_t mapSize = 1e5;
  try {
    db = new DB_("./db_rcc");
    dbIt = db->getIterator();
    dbIt.seekToFirst();
    if (!dbIt.isValid()) {
      printf("db is empty\n");
      return;
    }
    // Example: Insert a key-value pair
    // db->put("0xb0255cc0ad302545bd1dad0b39af5b1492f7a4b3",
    // "4348029885273000000000000");

    // Example: Iterate over database entries
  } catch (const std::runtime_error& e) {
    dbIt.reset();
    if (db) {
      delete db;
      db = NULL;
    }
    std::cerr << e.what() << std::endl;
    return;
  }

  DBMetaData metaData;
  try {
    metaData = readMetaData(db);
  } catch (const std::runtime_error& e) {
    std::cerr << e.what() << std::endl;
    return;
  }
  printf("lastBlock = %lu, recordCount = %lu\n", metaData.lastBlock,
         metaData.recordCount);

  ret = ecall_omap_init(global_eid, mapSize, metaData.recordCount);
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
  svr.listen("localhost", 1234);
  updateThread.join();

  if (db) {
    delete db;
    db = NULL;
  }

  // printf("Did not abort\n");
}
