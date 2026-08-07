#ifndef _PERF_COUNTERS_CALLGUARD_
#error Include "perf_counter.hpp" instead
#endif

// F(iscomputed, name, description, expression)
F(false, CIRCUITORAM_OVERFLOW, "Number of times ORAM stash was leaked due to overflow", 0)
F(false, CIRCUITORAM_STASH_HASH_BUILDS, "Number of one-time Circuit ORAM stash hash builds", 0)
F(false, CIRCUITORAM_STASH_HASH_CANDIDATE_OVERFLOWS, "Number of one-time Circuit ORAM stash hash candidate overflows", 0)
F(false, CIRCUITORAM_STASH_HASH_FALLBACKS, "Number of Circuit ORAM batches that fell back after all stash hash candidates overflowed", 0)
F(false, CIRCUITORAM_STASH_HASH_MAIN_ENTRIES_PROBED, "Number of one-time stash hash main entries probed", 0)
F(false, CIRCUITORAM_STASH_HASH_OVERFLOW_ENTRIES_PROBED, "Number of one-time stash hash overflow entries probed", 0)
F(false, CIRCUITORAM_STASH_HASH_RESTORES, "Number of one-time Circuit ORAM stash hash restores", 0)
F(false, OHMAP_DEAMORT_OVERFLOW, "Number of times OHMAP stash size was leaked due to OHMAP deamortization stash overflow", 0)
