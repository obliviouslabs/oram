#ifndef _PERF_COUNTERS_CALLGUARD_
#error Include "perf_counter.hpp" instead
#endif

// F(iscomputed, name, description, expression)
F(false, ORAM_STACK_OVERFLOW, "Number of times ORAM stash was leaked due to overflow", 0)
F(false, OHMAP_DEAMORT_OVERFLOW, "Number of times OHMAP stash size was leaked due to OHMAP deamortization stash overflow", 0)
