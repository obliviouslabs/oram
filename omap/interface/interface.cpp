#include <stdint.h>
#include <stdio.h>
void HelloWorld(uint32_t num) { printf("Hello, world %u!\n", num); }

#include "interface/common_interface.hpp"
#include "external_memory/server/serverBackend.hpp"
void ResetBackend(uint64_t size) {
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
  EM::Backend::g_DefaultBackend = new EM::Backend::MemServerBackend(size);
}

void DeleteBackend() {
  if (EM::Backend::g_DefaultBackend) {
    delete EM::Backend::g_DefaultBackend;
  }
}