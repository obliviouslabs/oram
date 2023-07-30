#pragma once
#include "common/tracing/tracer.hpp"

struct TimeTracker {
  using BlockTracer_t = BlockTracer<TimeTracker, std::chrono::_V2::system_clock::time_point>;
  using eint_t = std::underlying_type_t<EventId>;
  static TimeTracker g_Tracker;

  double totals[TOTAL_TRACKERS];
  uint64_t count[TOTAL_TRACKERS];
  uint64_t blockRecCount;
  BlockTracer_t* stack[MAX_TRACKER_REC];
  
  explicit TimeTracker() : stack{0} {
    std::cout << "Building tracker with " << TOTAL_TRACKERS << " total tracers." << std::endl;
    for (int i=0; i<TOTAL_TRACKERS; i++) {
      totals[i] = 0;
      count[i] = 0;
    }
    blockRecCount = 0;
  }

  INLINE void BeginBlock(BlockTracer_t& that) {
    Assert(blockRecCount < MAX_TRACKER_REC);
    that.metadata = std::chrono::system_clock::now();
    stack[blockRecCount] = &that;
    blockRecCount += 1;
  }

  INLINE void EndBlock() {
    blockRecCount -= 1;
    stack[blockRecCount]->Finish();
  }
  
  INLINE void Measure(const EventId eventId, const std::chrono::_V2::system_clock::time_point& start) {
    auto end = std::chrono::system_clock::now();
    if (g_disableTracing) return;
    const eint_t eid = static_cast<eint_t>(eventId);
    std::chrono::duration<double> diff = end - start;
    count[eid]++;
    totals[eid] += diff.count();
  }

  ~TimeTracker() {
    Reset();
  }

  void Log() {
    std::cout << "TimeTracker::Log:" << std::endl;
    #define F(name) { \
      eint_t i = static_cast<eint_t>(EventId::name); \
      if (count[i] != 0) { \
        std::cout << (#name) \
          << " " << totals[i] / (std::max((uint64_t)1,count[i])) \
          << " " << totals[i] \
          << " " << count[i] \
          << std::endl; \
      }}
    #include "../tracer_events.hxx"
  }

  void Reset(bool log=true) {
    if (log) {
      Log();
    }
    for (int i=0; i<TOTAL_TRACKERS; i++) {
      totals[i] = 0;
      count[i] = 0;
    }
  }
};

