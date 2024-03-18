#include "common/tracing/tracer.hpp"

#ifndef ENCLAVE_MODE

void print_stacktrace(int signum) {
  ::signal(signum, SIG_DFL);
  std::cerr << "Stack trace:\n" << boost::stacktrace::stacktrace() << '\n';
  Profiler::g_Tracker.Reset();
  TimeTracker::g_Tracker.Reset();
  g_MaxTracker.Reset();
  ::raise(SIGABRT);
}

TimeTracker TimeTracker::g_Tracker;
Profiler Profiler::g_Tracker;

MaxTracker g_MaxTracker;
bool g_disableTracing = false;
bool g_disableProfiling = true;
OnExitHandlerInstaller g_OnExit;

#endif

PerfCounters g_PerfCounters;
