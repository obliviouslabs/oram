searchfiles()

include_directories (${OPENSSL_INCLUDE_DIR})
add_library (common common/encutils.cpp common/tracing/tracer.cpp external_memory/server/serverBackend.cpp)
target_compile_definitions(common PRIVATE BOOST_STACKTRACE_USE_ADDR2LINE)
target_include_directories(common PUBLIC .)
target_include_directories(common PUBLIC ${BEARSSL_INCLUDE_DIRS})
target_link_libraries (common dl OpenSSL::SSL OpenSSL::Crypto ${BEARSSL_LIBRARY})

