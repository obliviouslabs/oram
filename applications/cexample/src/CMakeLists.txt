cmake_minimum_required(VERSION 3.10)
project(COMAP LANGUAGES C)

set(CMAKE_C_STANDARD 99)

# Add the executable
add_executable(main test.c)

# Link against the C++ library
target_link_libraries(main PRIVATE oraminterface dl
gtest_main
common
OpenMP::OpenMP_CXX)