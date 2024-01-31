#include <stdio.h>

#include <thread>
#include <vector>
#include <string>
#define UNTRUSTED_ENV 1
#include "../../common.hpp"
#include "../App.h"
#include "Enclave_u.h"
#ifdef DISK_IO
#include "external_memory/server/enclaveFileServer_untrusted.hpp"
#else
#include "external_memory/server/enclaveMemServer_untrusted.hpp"
#endif
#include "kvDB.hpp"

// since most balances are small, it is more space efficient to directly store them as string
using DB_ = KV_DB<std::string, std::string>;

DB_* db;
typename DB_::Iterator dbIt;

uint64_t ocall_Fetch_Next_KV_Batch(uint8_t* data, uint64_t batchBytes) {
  uint64_t bytes = 0;
  for (;dbIt.isValid() && bytes < batchBytes - sizeof(key_type) - sizeof(val_type); dbIt.next()){
    key_type key;
    val_type val;
    std::string keyStr = dbIt.key();
    if (keyStr.substr(0, 2) != "0x") {
      continue;
    }
    // pad to 42 bytes
    if (keyStr.size() != 42) {
      int padLen = 42 - keyStr.size();
      keyStr.insert(2, padLen, '0');
    }
    key.set(keyStr);
    val.set(dbIt.value());

    memcpy(data + bytes, &key, sizeof(key_type));
    bytes += sizeof(key_type);
    memcpy(data + bytes, &val, sizeof(val_type));
    bytes += sizeof(val_type);
  }
  return bytes;
}

void ActualMain(void) {
  sgx_status_t ret = SGX_ERROR_UNEXPECTED;
  size_t mapSize = 10000;
  try {
        db = new DB_("./db");
        dbIt = db->getIterator();
        dbIt.seekToFirst();
        if (!dbIt.isValid()) {
            printf("db is empty\n");
            return;
        }
        // Example: Insert a key-value pair
        // db->put("0xb0255cc0ad302545bd1dad0b39af5b1492f7a4b3", "4348029885273000000000000");
        

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

    uint64_t lastBlock = 0, recordCount = 0;
    try {
        lastBlock = std::stoull(db->get("lastBlock"));
        recordCount = std::stoull(db->get("recordCount"));
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        return;
    }
    printf("lastBlock = %lu, recordCount = %lu\n", lastBlock, recordCount);

  ret = ecall_omap_init(global_eid, mapSize, recordCount);
  if (ret != SGX_SUCCESS) abort();
  dbIt.reset();
  if (db) {
    delete db;
    db = NULL;
  }

  int findFlag = false;
  key_type key;
  val_type val;
  key.set(std::string("0xff2785eef607697b5fed88333a851d9add933625"));

  ret = ecall_omap_find(global_eid, &findFlag, (uint8_t*)&key, (uint8_t*)&val,
                        sizeof(key_type), sizeof(val_type));
  if (ret != SGX_SUCCESS) abort();
  if (findFlag == 0) {
    std::cout << "not found" << std::endl;
    abort();
  }
  std::cout << "found: " << val << std::endl;
  val.set(std::string("123456"));
  ret = ecall_omap_insert(global_eid, &findFlag, (uint8_t*)&key, (uint8_t*)&val,
                        sizeof(key_type), sizeof(val_type));
  if (ret != SGX_SUCCESS) abort();
  if (findFlag == 0) {
    std::cout << "not found" << std::endl;
    abort();
  }
  ret = ecall_omap_find(global_eid, &findFlag, (uint8_t*)&key, (uint8_t*)&val,
                        sizeof(key_type), sizeof(val_type));
  if (ret != SGX_SUCCESS) abort();
  if (findFlag == 0) {
    std::cout << "not found" << std::endl;
    abort();
  }
  std::cout << "found: " << val << std::endl;

  // printf("Did not abort\n");
}
