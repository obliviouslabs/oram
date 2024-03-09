searchfiles()

include_directories (${PROJECT_SOURCE_DIR}/osort)
include_directories (${OPENSSL_INCLUDE_DIR})

add_executable(
  test_par_omap
  par_omap.cpp
)
target_compile_definitions(test_par_omap PRIVATE BOOST_STACKTRACE_USE_ADDR2LINE)
target_link_libraries(
  test_par_omap
  dl
  gtest_main
  common
  OpenMP::OpenMP_CXX
)

add_executable(
  test_omap
  omap.cpp
)
target_compile_definitions(test_omap PRIVATE BOOST_STACKTRACE_USE_ADDR2LINE)
target_link_libraries(
  test_omap
  dl
  gtest_main
  common
)

add_executable(
  test_cuckoo
  cuckoo.cpp
)
target_compile_definitions(test_cuckoo PRIVATE BOOST_STACKTRACE_USE_ADDR2LINE)
target_link_libraries(
  test_cuckoo
  dl
  gtest_main
  common
  OpenMP::OpenMP_CXX
)

add_executable(
  test_basic
  basic.cpp
)
target_compile_definitions(test_basic PRIVATE BOOST_STACKTRACE_USE_ADDR2LINE)
target_link_libraries(
  test_basic
  dl
  gtest_main
  common
)

add_executable(
  test_oram
  oram.cpp
)
target_compile_definitions(test_oram PRIVATE BOOST_STACKTRACE_USE_ADDR2LINE)
target_link_libraries(
  test_oram
  dl
  gtest_main
  common
  OpenMP::OpenMP_CXX
)

add_executable(
  test_algo
  sort.cpp
)
target_compile_definitions(test_algo PRIVATE BOOST_STACKTRACE_USE_ADDR2LINE)
target_link_libraries(
  test_algo
  dl
  gtest_main
  common
)


add_executable(
  test_basic_perf
  perf_sort.cpp
)
target_compile_definitions(test_basic_perf PRIVATE BOOST_STACKTRACE_USE_ADDR2LINE)
target_link_libraries(
  test_basic_perf
  dl
  gtest_main
  common
)

include(GoogleTest)
gtest_discover_tests(test_algo)
gtest_discover_tests(test_basic_perf)
