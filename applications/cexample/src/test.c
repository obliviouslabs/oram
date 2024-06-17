#include <stdio.h>
#include <string.h>

#include "interface/omap_c.h"

void testOMap() {
  OblivMap omap = NewOblivMap(20, 32);
  InitEmpty(omap, 100);
  StartInit(omap, 100);
  uint8_t key[20] = {0};
  uint8_t val[32] = {0};
  memcpy(key, "key1", 4);
  memcpy(val, "val1", 4);
  Insert(omap, key, val);
  FinishInit(omap);
  memcpy(key, "key2", 4);
  memcpy(val, "val2", 4);
  Insert(omap, key, val);
  memcpy(key, "key1", 4);
  uint8_t val1[32] = {0};
  Find(omap, key, val1);
  printf("val1: %s\n", val1);
  uint8_t val2[32] = {0};
  memcpy(key, "key2", 4);
  Find(omap, key, val2);
  printf("val2: %s\n", val2);
  Destroy(omap);
}

int main() {
  testOMap();
  return 0;
}