searchfiles()

include_directories (${OPENSSL_INCLUDE_DIR})
add_library (common common/encutils.cpp common/tracing/tracer.cpp external_memory/server/serverBackend.cpp)
target_compile_definitions(common PRIVATE BOOST_STACKTRACE_USE_ADDR2LINE)
target_include_directories(common PUBLIC .)
target_link_libraries (common dl OpenSSL::SSL OpenSSL::Crypto)

add_library (oraminterface interface/recoram_interface.cpp common)
target_compile_definitions(oraminterface PRIVATE BOOST_STACKTRACE_USE_ADDR2LINE)
target_include_directories(oraminterface PUBLIC .)
target_link_libraries (oraminterface dl OpenSSL::SSL OpenSSL::Crypto OpenMP::OpenMP_CXX)
