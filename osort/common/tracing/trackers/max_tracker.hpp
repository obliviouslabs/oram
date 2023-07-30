#pragma once
#include "common/tracing/tracer.hpp"

struct MaxTracker {
  using eint_t = std::underlying_type_t<EventId>;

  uint64_t maxes[TOTAL_TRACKERS];
  uint64_t counts[TOTAL_TRACKERS];
  uint64_t totals[TOTAL_TRACKERS];
  
  MaxTracker() {
    std::cout << "Building MaxTracker" << std::endl;
    for (int i=0; i<TOTAL_TRACKERS; i++) {
      maxes[i] = 0;
      counts[i] = 0;
      totals[i] = 0;
    }
  }
  
  INLINE void Update(const EventId& eventId, uint64_t val) {
    if (g_disableTracing) return;
    const eint_t eid = static_cast<eint_t>(eventId);
    maxes[eid] = std::max(maxes[eid], val);
    totals[eid] += val;
    counts[eid] ++;
  }

  ~MaxTracker() {
    Reset();
  }

  void Log() {
    #define F(name) { \
      eint_t i = static_cast<eint_t>(EventId::name); \
      if (counts[i] != 0) { \
        std::cout << (#name) \
                  << " " << maxes[i] \
                  << " " << totals[i] / (double)counts[i]\
                  << std::endl; \
      }}
    #include "../tracer_events.hxx"
  }

  void Reset() {
    Log();
    for (int i=0; i<TOTAL_TRACKERS; i++) {
      maxes[i] = 0;
      counts[i] = 0;
      totals[i] = 0;
    }
  }
};

extern MaxTracker g_MaxTracker;
