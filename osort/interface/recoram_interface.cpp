#include "oram/omap.hpp"
#include "oram/par_omap.hpp"
#include "oram/recursive_oram.hpp"
#include "interface/recoram_interface.hpp"

using namespace ODSL;

ORAMBindingSingleton::ORAMBindingSingleton() {    
  uint64_t BackendSize = 1e10;
  EM::Backend::g_DefaultBackend = new EM::Backend::MemServerBackend(BackendSize);
  oram = nullptr;
}

void ORAMBindingSingleton::InitORAM(uint64_t size) {
  // ASSERT(oram == nullptr);
  oram = (void*) (new RecursiveORAM<uint64_t, uint32_t>(size));
}

void ORAMBindingSingleton::Write(uint32_t addr, uint64_t val) {
  ((RecursiveORAM<uint64_t, uint32_t>*)oram)->Write(addr, val);
}

uint64_t ORAMBindingSingleton::Read(uint32_t addr) {
  uint64_t ret;
  ((RecursiveORAM<uint64_t, uint32_t>*)oram)->Read(addr, ret);
  return ret;
}

