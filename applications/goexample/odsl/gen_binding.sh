rm -rf build # Needed after the CC/CXX export or after changing the CMAKE_BUILD_TYPE
export CC=/usr/bin/gcc
export CXX=/usr/bin/g++
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
ninja -C build
cp build/odsl.go odsl.go
LD_PATH=/usr/local/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LD_PATH
sudo cp go/libodsl.so $LD_PATH
cp go/libodsl.so $LD_PATH