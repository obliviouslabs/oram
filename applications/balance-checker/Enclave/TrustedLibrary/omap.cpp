#include "../Enclave.h"
#include "Enclave_t.h"
////
#include <functional>
#include <unordered_map>

#include "../../common.hpp"
#include "oram/omap.hpp"

using namespace ORAM;
EM::Backend::MemServerBackend* EM::Backend::g_DefaultBackend = nullptr;
OMap<uint64_t, int64_t> omap;

void ecall_omap_init(uint64_t N) {
  printf("init omap with size %lu\n", N);
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
  size_t BackendSize = N * 1000;
  EM::Backend::g_DefaultBackend =
      new EM::Backend::MemServerBackend(BackendSize);
  omap.SetSize(N);
  try {
    EM::VirtualVector::VirtualReader<std::pair<uint64_t, int64_t>> reader(
        N,
        [](uint64_t i) { return std::pair<uint64_t, int64_t>(i * 10, i * 3); });
    omap.InitFromReader(reader);

  } catch (std::exception& e) {
    printf("exception: %s\n", e.what());
  }
  return;
}

int ecall_omap_find(uint8_t* key, uint8_t* val, uint32_t keyLength,
                    uint32_t valLength) {
  const key_type& k = *reinterpret_cast<key_type*>(key);
  bool res = omap.find(k, *reinterpret_cast<val_type*>(val));
  return (int)res;
}