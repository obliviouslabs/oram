#pragma once
#include <stdint.h>
#include <stdio.h>

#include "omap_generic_interface.hpp"
void HelloWorld(uint32_t num) { printf("Hello, world %u!\n", num); }

#include "external_memory/server/serverBackend.hpp"
#include "interface/common_interface.hpp"
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
