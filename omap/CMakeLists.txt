searchfiles()

add_library (common common/encutils.cpp common/tracing/tracer.cpp external_memory/server/serverBackend.cpp)
target_compile_definitions(common PRIVATE BOOST_STACKTRACE_USE_ADDR2LINE)
target_include_directories(common PUBLIC .)
# HACK: in ubuntu bearssl is included as /bearsslrootdir whereas in other distros it is included as /bearsslrootdir/bearssl. This makes it work on both even though we will only use one of the includes on each distro:
target_include_directories(common PUBLIC ${BEARSSL_INCLUDE_DIRS} ${BEARSSL_INCLUDE_DIRS}/bearssl)
target_link_libraries (common dl ${BEARSSL_LIBRARY})

