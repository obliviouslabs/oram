#pragma once
#include "common/cpp_extended.hpp"

#ifndef SERVER__CACHE_SIZE
#define SERVER__CACHE_SIZE 8192
#endif

#define ORAM_USE_INRAM_SERVER true

#ifndef ENCLAVE_SIZE
#define ENCLAVE_SIZE 128
#endif
#ifndef DEFAULT_HEAP_SIZE
#define DEFAULT_HEAP_SIZE ((uint64_t)ENCLAVE_SIZE * 0xD0000UL)
#endif

#define FLAMEGRAPHS_BASE_FOLDER "./quality/flamegraphs/"