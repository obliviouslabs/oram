# Oblivious Map
An implementation of external memory efficient, cpu instruction and memory access trace parallel oblivious maps.

# Prerequisites
Install cmake, ninja and intel sgx sdk, or use the cppbuilder docker image.

## How to build the builder docker image
```bash
docker build -t cppbuilder:latest ./tools/docker/cppbuilder
```

## How to enter the docker environment to run unit tests
```bash
docker run -it --rm -v $PWD:/builder -u $(id -u) cppbuilder
```

## How to enter the docker environment to run algorithms in enclave
```bash
docker run -v /tmp/omapbackend:/ssdmount --privileged -it --rm -v $PWD:/builder -p 8080:8080 cppbuilder
```

## How to run the unit tests
```bash
rm -rf build
cmake -B build -G Ninja
ninja -C build
ninja -C build test
```

## How to run the unit tests in release mode

```bash
rm -rf build # Needed after the CC/CXX export or after changing the CMAKE_BUILD_TYPE
export CC=/usr/bin/gcc
export CXX=/usr/bin/g++
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
ninja -C build
```

## Build the omap example enclave (hardware mode)
```bash
source /startsgxenv.sh
cd applications/omap
make
```

## Build the omap example enclave (simulation mode)
```bash
source /startsgxenv.sh
cd applications/omap
make SGX_MODE=SIM
```

## Run a sample script to test runtime of omap algorithms
```bash
cd applications/omap
./algo_runner.sh
```

## Folder structure high level details

omap - C++ oblivious map library code
tests - C++ tests modules
applications - Enclaves example of omap
tools - tools used to generate graphs or test sets
tools/docker - dockerfiles used for reproducible builds

### omap folder structure

odsl - core library code of oblivious data structures
algorithm - algorithmic building blocks for oblivious data structures
common - common c++ utilies, cpu abstractions, cryptography abstractions and tracing code
external_memory - external memory vector abstraction
external_memory/server - server abstraction for different external memory scenarios (sgx, file system, ram)


### Profiling code

1) Compile with ENABLE_PROFILING

2) For the functions that need profiling, add PROFILE_F(); at as the first line of the function code. Additionally add the function name to trace_events.hxx

3) Use PROFILER_SET(false); to disable profiling, use PROFILER_RESET() to write the profile to the log file (see profiling related functions in profiling_test to confirm).

4) Use any of the tools in "Links to view flamegraph files" above to look at the profiling, adjust uncached IO time based on the results of the benchmarks enclave (benchmark_sgx).

