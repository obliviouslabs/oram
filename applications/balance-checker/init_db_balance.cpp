#include <fstream>

#include "App/TrustedLibrary/kv_db.hpp"

using DB_ = KV_DB<std::string, std::string>;

void initFromBalanceFile(const char* DBPath, const char* balancePath) {
  DB_* db = new DB_(DBPath);
  std::ifstream balanceFile(balancePath);
  std::string line;
  uint64_t lastBlock = 0, lastTxIdx = 0, recordCount = 0;
  std::getline(balanceFile, line);
  std::cout << line << std::endl;
  std::istringstream iss(line);
  if (!(iss >> lastBlock >> lastTxIdx >> recordCount)) {
    printf("Error reading balance file meta data\n");
    abort();
  }
  printf("lastBlock = %lu, recordCount = %lu\n", lastBlock, recordCount);
  while (std::getline(balanceFile, line)) {
    std::istringstream iss(line);
    std::string addr;
    std::string balance;
    if (!(iss >> addr >> balance)) {
      printf("Error reading balance file\n");
      abort();
    }
    if (addr.size() != 42) {
      int padLen = 42 - addr.size();
      addr.insert(2, padLen, '0');
    }
    db->put(addr, balance);
  }
  db->put("lastBlock", std::to_string(lastBlock + 10));
  db->put("lastTxIdx", std::to_string(lastTxIdx));
  db->put("recordCount", std::to_string(recordCount));
  delete db;
}

int main(int argc, char** argv) {
  if (argc != 3) {
    printf("Usage: %s <DBPath> <balancePath>\n", argv[0]);
    return 1;
  }
  initFromBalanceFile(argv[1], argv[2]);
  return 0;
}