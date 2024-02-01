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
#include "db_updater.hpp"
#include "kvDB.hpp"

// since most balances are small, it is more space efficient to directly store
// them as string
using DB_ = KV_DB<std::string, std::string>;

DB_* db;
typename DB_::Iterator dbIt;

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

void ActualMain(void) {
  sgx_status_t ret = SGX_ERROR_UNEXPECTED;
  size_t mapSize = 6e6;
  try {
    db = new DB_("./db_usdt");
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

  std::string logPath = "tx.log";
  while (true) {
    updateDBAndORAM(logPath);
    std::this_thread::sleep_for(std::chrono::seconds(5));
  }

  if (db) {
    delete db;
    db = NULL;
  }

  // printf("Did not abort\n");
}
