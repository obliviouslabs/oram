cmake_minimum_required(VERSION 3.10)
project(gobinding)

set(CMAKE_CXX_STANDARD 20)

# Find SWIG
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})
include(cmake/bearssl.cmake)
include(cmake/boost.cmake)

# Set SWIG module and include directories
set_source_files_properties(odsl.i PROPERTIES CPLUSPLUS ON)
include_directories("../../../omap")
set(CMAKE_SWIG_FLAGS "-c++" "-go" "-cgo" "-intgosize" "64")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG ")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -ggdb -fno-inline -no-pie")
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -march=native -mtune=native")



set(CPP_PARENT_PATH "../../../omap")


set(SOURCES odsl.i ${CPP_PARENT_PATH}/common/encutils.cpp ${CPP_PARENT_PATH}/common/tracing/tracer.cpp ${CPP_PARENT_PATH}/external_memory/server/serverBackend.cpp)

# Add SWIG module
swig_add_library(mylib TYPE SHARED LANGUAGE go SOURCES ${SOURCES})

# Set output name
set_target_properties(mylib PROPERTIES OUTPUT_NAME "libodsl")

# Set output directory
set_target_properties(mylib PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/go)

target_include_directories(mylib PUBLIC ${BEARSSL_INCLUDE_DIRS} ${BEARSSL_INCLUDE_DIRS}/bearssl)

# Link SWIG module and additional sources
target_link_libraries(mylib INTERFACE dl ${BEARSSL_LIBRARY})

# Generate Ninja build files
set_property(TARGET mylib PROPERTY SWIG_GENERATED TRUE)
set_property(SOURCE odsl.i PROPERTY GENERATED TRUE)
