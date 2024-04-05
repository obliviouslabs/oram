CMAKE_MINIMUM_REQUIRED(VERSION 3.8.0 FATAL_ERROR)

find_package(SWIG REQUIRED)
include(UseSWIG)
find_program(GO_EXECUTABLE go REQUIRED)

if(POLICY CMP0078)
  cmake_policy(SET CMP0078 NEW)
endif()
project(
  Example
  LANGUAGES CXX
)


#########################################################
# 第一步产生Example.so，这里只需要.cpp .h文件
#########################################################
add_library(Example SHARED "")

target_sources(Example
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/simpleclass.h>
    $<INSTALL_INTERFACE:include/simpleclass.h>
  PRIVATE
    src/simpleclass.cpp
)

target_include_directories(Example
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(Example PUBLIC cxx_std_11)

add_library(${PROJECT_NAME}::Example ALIAS Example)





############################################################
## 第二步产生wrap.cxx文件，这里只需要*.i及其include的头文件 
###########################################################
set(CMAKE_SWIG_FLAGS -cgo -intgosize 64)

set_property(SOURCE simplelib.i PROPERTY CPLUSPLUS ON) #没有这一条会生成wrap.c
set_property(SOURCE simplelib.i PROPERTY COMPILE_OPTIONS
  -package simplelib)

swig_add_library(exampleForGo
  LANGUAGE go
  TYPE OBJECT
  OUTPUT_DIR ${PROJECT_BINARY_DIR}/go/${PROJECT_NAME}
  SOURCES simplelib.i)


target_include_directories(exampleForGo PRIVATE ${PROJECT_SOURCE_DIR}/include)
#simplelib.i中有include “*.h”, 这里为这些头文件提供路径

set_target_properties(exampleForGo PROPERTIES
  SWIG_INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/include #SWIG也添加include路径
  SWIG_USE_TARGET_INCLUDE_DIRECTORIES ON
  POSITION_INDEPENDENT_CODE ON)
target_link_libraries(exampleForGo PRIVATE ${PROJECT_NAME}::Example)

#############################################################################
# 第三步复制一些文件到wrapper所在目录，供go build使用
#############################################################################
file(COPY "${CMAKE_SOURCE_DIR}/include/simpleclass.h" DESTINATION "${PROJECT_BINARY_DIR}/go/${PROJECT_NAME}")
file(COPY "${CMAKE_SOURCE_DIR}/lib.go" DESTINATION "${PROJECT_BINARY_DIR}/go/${PROJECT_NAME}")

#############################################################################
# 复制动态链接库libExample.so到系统目录 /usr/local/lib
#############################################################################

install(TARGETS Example DESTINATION ${CMAKE_INSTALL_LIBDIR})