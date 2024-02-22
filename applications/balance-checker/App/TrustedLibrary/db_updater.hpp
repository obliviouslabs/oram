#pragma once
#include <boost/multiprecision/cpp_int.hpp>
#include <fstream>
using boost::multiprecision::uint256_t;
#include "kv_db.hpp"
#include "lock_utils.hpp"
using DB_ = KV_DB<std::string, std::string>;
struct DBMetaData {
  uint64_t lastBlock;
  uint64_t lastTxIdx;
  uint64_t lastStableBlock;
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
  metaData.lastStableBlock = metaData.lastBlock - 10;
  return metaData;
}

// returns the last block scanned
uint64_t updateDBFromLog(
    DB_* db, const std::string& logPath,
    std::unordered_map<std::string, uint256_t>& unstableDiff,
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
  struct Tx {
    bool isDuplicate = false;
    uint64_t blockNumber;
    uint64_t txIndex;
    std::string from;
    std::string to;
    uint256_t value;
  };
  std::string line;
  bool changed = false;
  std::vector<Tx> txs;
  int lineCount = 0;
  while (std::getline(logFile, line)) {
    ++lineCount;
    std::string txHash;
    std::string timestamp;
    std::istringstream iss(line);
    txs.emplace_back();
    auto& lastTx = txs.back();
    if (!(iss >> lastTx.blockNumber >> txHash >> lastTx.txIndex >>
          lastTx.from >> lastTx.to >> lastTx.value >> timestamp)) {
      throw std::runtime_error("Error reading log file");
    }
    if (lastTx.blockNumber <= metaData.lastStableBlock) {
      txs.pop_back();
      continue;
    }

    // if (blockNumber < metaData.lastBlock ||
    //     (blockNumber == metaData.lastBlock && txIndex <= metaData.lastTxIdx))
    //     {
    //   continue;
    // }
    if (lastTx.from.length() < 2 || lastTx.from.substr(0, 2) != "0x") {
      printf("Invalid ethereum from address %s\n", lastTx.from.c_str());
      printf("line = %s\n", line.c_str());
      throw std::invalid_argument("Invalid ethereum address");
    }
    if (lastTx.to.length() < 2 || lastTx.to.substr(0, 2) != "0x") {
      printf("Invalid ethereum to address %s\n", lastTx.to.c_str());
      printf("line = %s\n", line.c_str());
      throw std::invalid_argument("Invalid ethereum address");
    }

    // metaData.lastBlock = blockNumber;
    // metaData.lastTxIdx = txIndex;
    // changed = true;
    if (lastTx.from.size() != 42) {
      int padLen = 42 - lastTx.from.size();
      lastTx.from.insert(2, padLen, '0');
    }
    if (lastTx.to.size() != 42) {
      int padLen = 42 - lastTx.to.size();
      lastTx.to.insert(2, padLen, '0');
    }
    // balanceMap[from] -= value;
    // balanceMap[to] += value;
  }
  printf("read %d lines from log, keeps %lu txs\n", lineCount, txs.size());
  if (!txs.empty() && txs.back().blockNumber >= metaData.lastBlock) {
    metaData.lastBlock = txs.back().blockNumber;
    metaData.lastTxIdx = txs.back().txIndex;
    metaData.lastStableBlock = metaData.lastBlock - 10;
  } else {
    clearFile(logFd);
    if (!txs.empty())
      printf("last tx block = %lu, prev last block = %lu\n",
             txs.back().blockNumber, metaData.lastBlock);
    return metaData.lastBlock;
  }
  uint64_t minBlock = metaData.lastBlock;
  uint64_t minIdx = metaData.lastTxIdx;
  int duplicateCount = 0;
  for (auto txIt = txs.rbegin(); txIt != txs.rend(); ++txIt) {
    if (txIt->blockNumber > minBlock ||
        (txIt->blockNumber == minBlock && txIt->txIndex >= minIdx)) {
      txIt->isDuplicate = true;
      ++duplicateCount;
      continue;
    }
    minBlock = txIt->blockNumber;
    minIdx = txIt->txIndex;
  }
  printf("duplicateCount = %d\n", duplicateCount);
  auto txIt = txs.begin();
  for (; txIt != txs.end(); ++txIt) {
    if (txIt->isDuplicate) {
      continue;
    }
    if (txIt->blockNumber > metaData.lastStableBlock) {
      break;
    }
    balanceMap[txIt->from] -= txIt->value;
    balanceMap[txIt->to] += txIt->value;
  }

  // update db with stable changes
  // ensure atomic write
  auto writeBatch = db->newWriteBatch();
  for (const auto& [addr, balance] : balanceMap) {
    // printf("addr = %s, diff balance = %s\n", addr.c_str(),
    //        balance.str().c_str());
    if (balance == 0) {
      continue;
    }
    std::string originalBalanceStr;
    bool found = db->get(addr, originalBalanceStr);
    uint256_t newBalance = balance;
    // printf("found = %d\n", found);
    if (found) {
      // printf("originalBalanceStr = %s\n", originalBalanceStr.c_str());
      uint256_t originalBalance(originalBalanceStr);
      newBalance += originalBalance;
    }
    if (newBalance == 0) {
      writeBatch.Delete(addr);
      --metaData.recordCount;
      // deleteList.push_back(addr);
      continue;
    } else if (found) {
      // updateList.push_back({addr, balance.str()});
    } else {
      ++metaData.recordCount;
      // insertList.push_back({addr, balance.str()});
    }
    writeBatch.Put(addr, newBalance.str());
  }
  writeBatch.Put("lastBlock", std::to_string(metaData.lastBlock));
  writeBatch.Put("lastTxIdx", std::to_string(metaData.lastTxIdx));
  writeBatch.Put("recordCount", std::to_string(metaData.recordCount));

  /* balanceMap == newStable - oldStable */

  // revert the previous unstable changes
  for (const auto& [addr, balance] : unstableDiff) {
    balanceMap[addr] -= balance;
  }
  /* balanceMap == newStable - oldUnstable */

  std::unordered_map<std::string, uint256_t> newUnstableDiff;
  // include the new unstable changes, and set the new unstable diff
  for (; txIt != txs.end(); ++txIt) {
    if (txIt->isDuplicate) {
      continue;
    }
    auto& tx = *txIt;
    balanceMap[tx.from] -= tx.value;
    balanceMap[tx.to] += tx.value;
    newUnstableDiff[tx.from] -= tx.value;
    newUnstableDiff[tx.to] += tx.value;
  }

  /* balanceMap == newUnstable - oldUnstable */

  for (auto& [addr, balance] : balanceMap) {
    if (balance == 0) {
      continue;
    }

    std::string dbOriginalBalanceStr;
    bool dbFound = db->get(addr, dbOriginalBalanceStr);

    // printf("found = %d\n", found);
    uint256_t omapOriginalBalance = unstableDiff[addr];
    if (dbFound) {
      // printf("originalBalanceStr = %s\n", originalBalanceStr.c_str());
      uint256_t dbOriginalBalance(dbOriginalBalanceStr);
      omapOriginalBalance += dbOriginalBalance;
    }
    bool found = omapOriginalBalance != 0;
    balance += omapOriginalBalance;
    if (balance == 0) {
      Assert(found);
      deleteList.push_back(addr);
    } else if (found) {
      updateList.push_back({addr, balance.str()});
    } else {
      insertList.push_back({addr, balance.str()});
    }
  }
  db->writeBatch(writeBatch);
  clearFile(logFd);
  std::swap(unstableDiff, newUnstableDiff);
  return metaData.lastBlock;
}
