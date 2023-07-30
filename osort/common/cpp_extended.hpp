#pragma once

#ifndef ENCLAVE_MODE
#include <chrono>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <csignal>
#include <ctime>
#define NENCLAVE(...) __VA_ARGS__
#else
#include <cstdio>
#define NENCLAVE(...)
// #define consteval constexpr
// For enclave mode, we assume the following ocalls are implemented:
// ocall_print_string(char*)
// ocall_FAIL()
//
#endif

#include <iostream>
#include <utility>
#include <inttypes.h>

// struct std::ostringstream;
// struct std::ostream;
#ifndef ENCLAVE_MODE
template<bool logTime, bool logFile>
inline void Log_Recursive(const char* file, int line, std::ostringstream& msg)
{
    
        auto now_c = std::chrono::system_clock::now();
        if constexpr (logTime) {
            std::cerr << "[" << now_c.time_since_epoch().count() / 1000000000 << "." << ((now_c.time_since_epoch().count() % 1000000000) / 1000) << "]";
        }
        if constexpr (logFile) {
            std::cerr << "[" << file << ":" << line << "]: ";
        }
        std::cerr << msg.str();

        // ocall_print_string(msg.str().c_str());

}

template<typename T>
inline void Log_One(std::ostringstream& msg, T value)
{
    msg << value;
}

// c++20 is not compatible with c++14 libraries for the << for unsigned int's, so we overload with to_string:
//

template<>
inline void Log_One<uint64_t>(std::ostringstream& msg, uint64_t value)
{
    msg << std::to_string(value);
}

template<>
inline void Log_One<uint32_t>(std::ostringstream& msg, uint32_t value)
{
    msg << std::to_string(value);
}

template<>
inline void Log_One<int32_t>(std::ostringstream& msg, int32_t value)
{
    msg << std::to_string(value);
}

template<>
inline void Log_One<int64_t>(std::ostringstream& msg, int64_t value)
{
    msg << std::to_string(value);
}

template<>
inline void Log_One<bool>(std::ostringstream& msg, bool value)
{
    msg << std::to_string(value);
}

// "Recursive" variadic function
template<bool logTime, bool logFile, typename T, typename... Args>
void Log_Recursive(const char* file, int line, std::ostringstream& msg,
                   T value, const Args&... args)
{
    Log_One(msg, value);
    Log_Recursive<logTime, logFile>(file, line, msg, args...);
}

template<bool logTime=true, bool logFile=true, typename... Args>
void LogWrapper(const char* file, int line, const Args&... args)
{
    std::ostringstream msg;
    Log_Recursive<logTime, logFile>(file, line, msg, args...);
}

#endif


template<typename T>
inline constexpr bool IS_POD() {
    return std::is_trivial<T>() && std::is_standard_layout<T>();
}

#define static_false(msg) []<bool flag = false>() {static_assert(flag, msg);}();


template <typename F>
struct _defer_class {
    _defer_class(F&& f) : _f(std::forward<F>(f)) {}
    ~_defer_class() { _f(); }
    typename std::remove_reference<F>::type _f;
};

template <typename F>
inline _defer_class<F> _create_defer_class(F&& f) {
    return _defer_class<F>(std::forward<F>(f));
}

#define _CONCAT_INNER(x, n) x##n
#define _CONCAT(x, n) _CONCAT_INNER(x, n)
#define _CONCAT_LINE(name) _CONCAT(name, __LINE__)

#define DEFER(e) \
    auto _CONCAT_LINE(_defer_var_) = _create_defer_class([&](){ e; })


#ifndef __clang__
    #define CLANG_OR_GCC(a,b) b
    #define CE_StringEquals(a,b) strcmp(a,b)==0
#else
    #define CLANG_OR_GCC(a,b) a
    #pragma clang diagnostic ignored "-Wc99-designator"
    #pragma clang diagnostic ignored "-Wunused-command-line-argument"
    consteval bool strings_equal(char const * a, char const * b) {
        return *a == *b && (*a == '\0' || strings_equal(a + 1, b + 1));
    }
    #define CE_StringEquals(a,b) strings_equal(a,b)
#endif

#ifndef ENCLAVE_MODE
    #ifndef __clang__            
        #include <source_location>
        #define SOURCE_LOCATION std::source_location        
    #else        
        #include <experimental/source_location>
        #define SOURCE_LOCATION std::experimental::fundamentals_v2::source_location    
    #endif
#endif


#define IGNORE_UNUSED(x)  (void)x;
#define INLINE inline __attribute__((always_inline))
#define RESTRICT __restrict__


#ifndef NDEBUG
#define DEBUG_ONLY(...) __VA_ARGS__
#else
#define DEBUG_ONLY(...)
#endif


// A preprocessor argument counter
#define COUNT(...) COUNT_I(__VA_ARGS__, 9, 8, 7, 6, 5, 4, 3, 2, 1,)
#define COUNT_I(_9,_8,_7,_6,_5,_4,_3,_2,_1,X,...) X
// Preprocessor paster
#define GLUE(A,B) GLUE_I(A,B)
#define GLUE_I(A,B) A##B
// chained caller
#define NAMED_VALUES(...) GLUE(NAMED_VALUES_,COUNT(__VA_ARGS__))(__VA_ARGS__)
// chain
#define NAMED_VALUES_1(a) ", " #a "=",a
#define NAMED_VALUES_2(a,...) ", " #a "=",a,NAMED_VALUES_1(__VA_ARGS__)
#define NAMED_VALUES_3(a,...) ", " #a "=",a,NAMED_VALUES_2(__VA_ARGS__)
#define NAMED_VALUES_4(a,...) ", " #a "=",a,NAMED_VALUES_3(__VA_ARGS__)
#define NAMED_VALUES_5(a,...) ", " #a "=",a,NAMED_VALUES_4(__VA_ARGS__)
#define NAMED_VALUES_6(a,...) ", " #a "=",a,NAMED_VALUES_5(__VA_ARGS__)
#define NAMED_VALUES_7(a,...) ", " #a "=",a,NAMED_VALUES_6(__VA_ARGS__)
#define NAMED_VALUES_8(a,...) ", " #a "=",a,NAMED_VALUES_7(__VA_ARGS__)
#define NAMED_VALUES_9(a,...) ", " #a "=",a,NAMED_VALUES_8(__VA_ARGS__)

#ifndef ENCLAVE_MODE
#define X_FAIL() std::raise(SIGABRT);
#else
#define X_FAIL() ocall_FAIL();
#endif

#ifndef NDEBUG
#ifndef ENCLAVE_MODE
#define Assert(expr, ...) if (expr) [[likely]] {} else { X_LOG("Assertion violated: {", #expr, "}" __VA_OPT__(, NAMED_VALUES(__VA_ARGS__))); X_FAIL(); }
#define X_LOG(...) LogWrapper(__FILE__, __LINE__, __VA_ARGS__, "\n")
#define X_LOG_SIMPLE(...) LogWrapper<false,false>(__FILE__, __LINE__, __VA_ARGS__)
#else
#define Assert(expr, ...) if (expr) [[likely]] {} else { printf("Assertion violated: { %s }\n", #expr); X_FAIL(); }
#define X_LOG(...) 
#define X_LOG_SIMPLE(...)
#endif

#else

#define Assert(...)
#define X_LOG(...) 
#define X_LOG_SIMPLE(...)

#endif

#ifndef NDEBUG
/* When debugging is enabled, these form aliases to useful functions */
#define dbg_printf(...) printf(__VA_ARGS__)

#else
/* When debugging is disabled, no code gets generated for these */
#define dbg_printf(...)
#endif

#define STMT( stuff ) do { stuff } while (false);
#define TRACE_FUNCTION(...)  STMT( \
X_LOG_SIMPLE("[>> ", __func__," ", NAMED_VALUES(__VA_ARGS__), "]\n"); \
); \
auto fname = __func__; \
DEFER( STMT( X_LOG("[<< ", fname, "]"); ); );

#define COMMA ,