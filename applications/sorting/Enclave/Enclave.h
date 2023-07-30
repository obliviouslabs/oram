#pragma once

#include <stdlib.h>
#include <assert.h>
#ifndef NOOPENSSL
    #define NOOPENSSL
#endif
#if defined(__cplusplus)
extern "C" {
#endif

void printf(const char *fmt, ...);

#if defined(__cplusplus)
}
#endif
