#include "../Enclave.h"
#include "Enclave_t.h"
////
#include <functional>
#include <unordered_map>

#include "../../common.hpp"
#include "oram/omap.hpp"

using namespace ORAM;
EM::Backend::MemServerBackend* EM::Backend::g_DefaultBackend = nullptr;
OMap<key_type, val_type> omap;

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