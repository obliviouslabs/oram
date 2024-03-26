#include "interface/recoram_interface.hpp"

#include "odsl/omap.hpp"
#include "odsl/par_omap.hpp"
#include "odsl/recursive_oram.hpp"

using namespace ODSL;

ORAMBindingSingleton::ORAMBindingSingleton() {
  uint64_t BackendSize = 1e10;
  EM::Backend::g_DefaultBackend =
      new EM::Backend::MemServerBackend(BackendSize);
  oram = nullptr;
}

void ORAMBindingSingleton::InitORAM(uint32_t size) {
  // ASSERT(oram == nullptr);
  oram = (void*)(new RecursiveORAM<uint64_t, uint32_t>(size));
  ((RecursiveORAM<uint64_t, uint32_t>*)oram)->InitDefault(0);
}

void ORAMBindingSingleton::Write(uint32_t addr, uint64_t val) {
  ((RecursiveORAM<uint64_t, uint32_t>*)oram)->Write(addr, val);
}

uint64_t ORAMBindingSingleton::Read(uint32_t addr) {
  uint64_t ret;
  ((RecursiveORAM<uint64_t, uint32_t>*)oram)->Read(addr, ret);
  return ret;
}
