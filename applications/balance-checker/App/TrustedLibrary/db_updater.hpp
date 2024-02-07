#pragma once
#include <boost/multiprecision/cpp_int.hpp>
#include <fstream>
using boost::multiprecision::uint256_t;
#include "lock_utils.hpp"
#include "kvDB.hpp"
using DB_ = KV_DB<std::string, std::string>;
struct DBMetaData {
  uint64_t lastBlock;
  uint64_t lastTxIdx;
  uint64_t recordCount;
};

DBMetaData readMetaData(DB_* db) {
  std::string lastBlockStr, lastTxIdxStr, recordCountStr;
  bool lastBlockFlag = db->get("lastBlock", lastBlockStr);
  if (!lastBlockFlag) {
    throw std::runtime_error("lastBlock not specified");
  }
  bool recordCountFlag = db->get("recordCount", recordCountStr);
  if (!recordCountFlag) {
    throw std::runtime_error("recordCount not specified");
  }

  bool lastTxIdxFlag = db->get("lastTxIdx", lastTxIdxStr);
  if (!lastTxIdxFlag) {
    throw std::runtime_error("lastTxIdx not specified");
  }
  DBMetaData metaData;
  metaData.lastBlock = std::stoull(lastBlockStr);
  metaData.recordCount = std::stoull(recordCountStr);
  metaData.lastTxIdx = std::stoull(lastTxIdxStr);
  return metaData;
}

// returns the last block scanned
uint64_t updateDBFromLog(
    DB_* db, const std::string& logPath,
    std::vector<std::pair<std::string, std::string>>& insertList,
    std::vector<std::pair<std::string, std::string>>& updateList,
    std::vector<std::string>& deleteList) {
  FileLocker logLocker(logPath);
  int logFd = logLocker.getFd();
  std::ifstream logFile(logPath);

  DBMetaData metaData = readMetaData(db);
  std::unordered_map<std::string, uint256_t> balanceMap;
  // log file in format: blockNumber, txHash, txIndex, from, to, value,
  // timestamp e.g. 16455165
  // 0x8a32784c5a32d1fdf2aff9c1a977c0389a536bc6550fe65711dc72c977bcd4c1	157
  // 0x8d12A197cB00D4747a1fe03395095ce2A5CC6819
  // 0xC44D85575607a609C1d7F49819754722ca0BbC97	38863933877536673917104
  // 2023-01-21T12:26:35

  std::string line;
  bool changed = false;
  while (std::getline(logFile, line)) {
    std::istringstream iss(line);
    std::string txHash, from, to, timestamp;
    uint64_t blockNumber, txIndex;
    uint256_t value;
    if (!(iss >> blockNumber >> txHash >> txIndex >> from >> to >> value >>
          timestamp)) {
      throw std::runtime_error("Error reading log file");
    }

    if (blockNumber < metaData.lastBlock ||
        (blockNumber == metaData.lastBlock && txIndex <= metaData.lastTxIdx)) {
      continue;
    }

    Assert(from.substr(0, 2) != "0x", "from address not hex");
    Assert(to.substr(0, 2) != "0x", "to address not hex");
    metaData.lastBlock = blockNumber;
    metaData.lastTxIdx = txIndex;
    changed = true;
    if (from.size() != 42) {
      int padLen = 42 - from.size();
      from.insert(2, padLen, '0');
    }
    if (to.size() != 42) {
      int padLen = 42 - to.size();
      to.insert(2, padLen, '0');
    }
    balanceMap[from] -= value;
    balanceMap[to] += value;
  }
  if (!changed) {
    // prune log
    clearFile(logFd);
    return metaData.lastBlock;
  }
  // ensure atomic write
  auto writeBatch = db->newWriteBatch();
  for (auto& [addr, balance] : balanceMap) {
    // printf("addr = %s, diff balance = %s\n", addr.c_str(),
    //        balance.str().c_str());
    if (balance == 0) {
      continue;
    }
    std::string originalBalanceStr;
    bool found = db->get(addr, originalBalanceStr);
    // printf("found = %d\n", found);
    if (found) {
      // printf("originalBalanceStr = %s\n", originalBalanceStr.c_str());
      uint256_t originalBalance(originalBalanceStr);
      balance += originalBalance;
    }
    // TODO: delete address with zero balance when oram supports delete
    // if (balance == 0) {
    //   writeBatch.Delete(addr);
    //   --metaData.recordCount;
    //   deleteList.push_back(addr);
    //   continue;
    // } else
    if (found) {
      updateList.push_back({addr, balance.str()});
    } else {
      ++metaData.recordCount;
      insertList.push_back({addr, balance.str()});
    }
    writeBatch.Put(addr, balance.str());
  }
  writeBatch.Put("lastBlock", std::to_string(metaData.lastBlock));
  writeBatch.Put("lastTxIdx", std::to_string(metaData.lastTxIdx));
  writeBatch.Put("recordCount", std::to_string(metaData.recordCount));
  db->writeBatch(writeBatch);
  clearFile(logFd);
  return metaData.lastBlock;
}
