#pragma once
#include <stdint.h>

void ResetBackend(uint64_t size);
void DeleteBackend();

void HelloWorld(uint32_t num);

template <typename K>
struct TestKey {
  K key;
  uint32_t val;
};