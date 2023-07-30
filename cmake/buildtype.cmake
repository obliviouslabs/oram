set(default_build_type "Debug")
# set(default_build_type "Release")
# if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
#   set(default_build_type "Debug")
# endif()

if(ENCLAVE_BUILD AND ENCLAVE_BUILD EQUAL "1")
  # We only do enclave builds in Release mode
  #
  # We build in CXX20 mode, but in actual builds we will use stl from sgx, which 
  # is c++14 stl minus all functions that cannot exists inside an enclave.
  #
  set(CMAKE_CXX_STANDARD 20)
  set(default_build_type "Release")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DENCLAVE_MODE")
else()
  set(CMAKE_CXX_STANDARD 20)
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()


if("${CMAKE_BUILD_TYPE}" EQUAL "Release")
  # This flag breaks addr2line:
  #
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

message(STATUS "Build type is ${CMAKE_BUILD_TYPE} (default is '${default_build_type}').")