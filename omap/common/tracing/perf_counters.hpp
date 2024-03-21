#ifndef F
#error Define F before including this file
#endif
#define _PERF_COUNTERS_CALLGUARD_ 1

#ifndef NDEBUG
#include "perf_counters.hxx"
#endif

#include "perf_counters_release.hxx"

#undef _PERF_COUNTERS_CALLGUARD_
#undef F