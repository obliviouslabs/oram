get_property(INCLUDE_FILES GLOBAL PROPERTY INCLUDE_FILES_G)
get_property(SOURCE_FILES_NM GLOBAL PROPERTY SOURCE_FILES_NM_G)

add_custom_target(
        cppcheck
        COMMAND cppcheck
        --enable=warning,performance,portability,style,information,missingInclude
        --std=c++20
        --library=googletest.cfg
        # --template="[{severity}][{id}] {message} {callstack} \(On {file}:{line}\)"
        -I../omap/ -I../tests/
        --verbose
        --quiet
        ${SOURCE_FILES_NM}
)

add_custom_target(
        cppcheck-builder
        COMMAND cppcheck 
        --xml --xml-version=2 
        --enable=all 
        -I../omap/ -I../tests/
        --std=c++20 
        --library=googletest.cfg
        -DBOOST_STACKTRACE_USE_ADDR2LINE 
        ${SOURCE_FILES_NM} --suppress=missingIncludeSystem
        --check-config 
        2> quality/cppcheck.xml
        COMMAND 
        cppcheck-codequality
        --input-file=quality/cppcheck.xml
        --output-file=quality/cppcheck.json
)