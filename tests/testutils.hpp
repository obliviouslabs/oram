#include <gtest/gtest.h>
#include "common/tracing/tracer.hpp"
#ifdef NDEBUG
  #define ISDEBUG false
#else
  #define ISDEBUG true
#endif

#define DEBUG_ONLY_TEST() if (!ISDEBUG) { GTEST_SKIP(); } 
#define RELEASE_ONLY_TEST() if (ISDEBUG) { GTEST_SKIP(); } 

class FooEnvironment : public ::testing::Environment {
 public:
  ~FooEnvironment() override {}

  // Override this to define how to set up the environment.
  void SetUp() override {
    g_OnExit.Init();
  }

  // Override this to define how to tear down the environment.
  void TearDown() override {}
};

::testing::Environment* const foo_env = ::testing::AddGlobalTestEnvironment(new FooEnvironment);
