#pragma once
#include <iostream>

#ifndef ENCLAVE_MODE
// #define ENABLE_PERF_COUNTERS 1
#endif

struct PerfCounters {
#ifdef ENABLE_PERF_COUNTERS
#define F(_, name, ...) uint64_t name = 0;
#include "common/tracing/perf_counters.hpp"

  void Reset() {
#define F(_, name, ...) this->name = 0;
#include "common/tracing/perf_counters.hpp"
  }

  void Log(std::ostream& ofs) {
#define F(iscomputed, name, description, expression)         \
  if constexpr (iscomputed) {                                \
    ofs << description << ": " << (expression) << std::endl; \
  } else {                                                   \
    ofs << description << ": " << (name) << std::endl;       \
  }
#include "common/tracing/perf_counters.hpp"
  }
#endif
};

#ifdef ENABLE_PERF_COUNTERS
#define PERFCTR_RESET() g_PerfCounters.Reset()
#define PERFCTR_LOG() g_PerfCounters.Log(std::cout)
#define PERFCTR_LOG_TO(to) g_PerfCounters.Log(to)
#define PERFCTR_INCREMENT(name) g_PerfCounters.name++;
#define PERFCTR_INCREMENT_BY(name, by) g_PerfCounters.name += (by);
#else
#define PERFCTR_RESET()
#define PERFCTR_LOG()
#define PERFCTR_LOG_TO(...)
#define PERFCTR_INCREMENT(...)
#define PERFCTR_INCREMENT_BY(...)
#endif

extern PerfCounters g_PerfCounters;