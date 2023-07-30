#pragma once
#include <chrono>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <string>
#include <type_traits>
#include <csignal>
#include <cstring>

#include "common/defs.hpp"

#ifndef ENCLAVE_MODE
  #include <boost/stacktrace.hpp>

  extern bool g_disableTracing;


  enum class EventId : int {
    #define F(name)  name,
    #include "./tracer_events.hxx"
    LAST_EVENT
  };

  consteval EventId CE_GetFuncEventId(const char* funcName = SOURCE_LOCATION::current().function_name()) {
    #define F(name) if (CE_StringEquals(funcName, # name)) return EventId::name;
    #include "./tracer_events.hxx"
    return EventId::LAST_EVENT;
  }

  // BlockTracer:
  // 
  template<typename T, typename Metadata>
  struct BlockTracer {
    Metadata metadata;
    EventId eventId;

    BlockTracer() = default;
    BlockTracer( const BlockTracer& ) = delete; // non construction-copyable
    BlockTracer& operator=( const BlockTracer& ) = delete; // non copyable

    INLINE BlockTracer(const EventId eventId) : eventId(eventId) {
      Assert(eventId != EventId::LAST_EVENT);
      T::g_Tracker.BeginBlock(*this);
    }

    INLINE void Finish() {
      Assert(eventId != EventId::LAST_EVENT);
      T::g_Tracker.Measure(eventId, metadata);
      eventId = EventId::LAST_EVENT;
    }

    INLINE ~BlockTracer() {
      if (eventId != EventId::LAST_EVENT) {
        T::g_Tracker.EndBlock();
      }
    }
  };

  constexpr uint64_t TOTAL_TRACKERS = static_cast<std::underlying_type_t<EventId> >(EventId::LAST_EVENT);
  constexpr uint64_t MAX_TRACKER_REC = 128;


  /// Code related to initialization
  //
  extern void print_stacktrace(int signum);
  struct OnExitHandlerInstaller {
    int theStart;

    void Init() {
      theStart = 10;
      std::cerr << "Installing signal handler" << std::endl;
      ::signal(SIGSEGV, &print_stacktrace);
      ::signal(SIGABRT, &print_stacktrace);
      ::signal(SIGFPE, &print_stacktrace);
    }

    OnExitHandlerInstaller() {
      Init();
    }
  };
  extern OnExitHandlerInstaller g_OnExit;
  //
  //
  

  #include "trackers/time_tracker.hpp"
  #include "trackers/max_tracker.hpp"
  #include "trackers/profiler.hpp"

#endif

#ifdef ENCLAVE_MODE
#ifdef ENABLE_PROFILING
#undef ENABLE_PROFILING
#endif
#endif

#ifdef ENABLE_PROFILING
#define PROFILER_SET(val) g_disableProfiling = !val;
#define PROFILER_RESET(...) Profiler::g_Tracker.Reset(__VA_ARGS__);
#define PROFILE_F() typename Profiler::BlockProfiler_t _CONCAT_LINE(profiler)(CE_GetFuncEventId(__FUNCTION__));
#define PROFILE_B(name) typename Profiler::BlockProfiler_t _CONCAT_LINE(_profiler_##name)(EventId::name);
#define END_PROFILE_B() Profiler::g_Tracker.EndBlock();
#define SWITCH_PROFILE_B(name) END_PROFILE_B(); PROFILE_B(name);
#else
#define PROFILER_SET(val) 
#define PROFILER_RESET(...)
#define PROFILE_F()
#define PROFILE_B(name) 
#define END_PROFILE_B() 
#define SWITCH_PROFILE_B(name) 
#endif

#if defined ENABLE_PROFILING || defined ENABLE_VALUE_TRACKING
#define PROFILE_V(varName, val) Profiler::g_Tracker.TrackValue(EventId::varName, val)
#else
#define PROFILE_V(varName, val)
#endif

#if defined ENABLE_TRACING
#define TRACER_SET(val) g_disableTracing = !val
#define TRACER_RESET() TimeTracker::g_Tracker.Reset(); g_MaxTracker.Reset();
#define TRACE_MAX(varName, val) g_MaxTracker.Update(EventId::varName, val);
#define TRACE_F() typename TimeTracker::BlockTracer_t _CONCAT_LINE(tracer)(CE_GetFuncEventId(__FUNCTION__));
#define TRACE_B(name) typename TimeTracker::BlockTracer_t _CONCAT_LINE(_tracer_##name)(EventId::name);
#define END_TRACE_B() TimeTracker::g_Tracker.EndBlock();
#define SWITCH_TRACE_B(name) END_TRACE_B(); TRACE_B(name);
#else
#define TRACER_SET(val)
#define TRACER_RESET()
#define TRACE_F()
#define TRACE_B(name) 
#define END_TRACE_B() 
#define SWITCH_TRACE_B(name) 
#define TRACE_MAX(varName, val)
#endif

#include "common/tracing/perf.hpp"