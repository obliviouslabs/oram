rm -rf build # Needed after the CC/CXX export or after changing the CMAKE_BUILD_TYPE
export CC=/usr/bin/gcc
export CXX=/usr/bin/g++
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
ninja -C build
cp build/odsl.go odsl.go
sudo cp go/libodsl.so /usr/local/lib/