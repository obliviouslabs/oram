#!/bin/bash

ROCKSDB_VERSION=8.10.0

#Run as a root user
if [ "$EUID" -ne 0 ]
  then echo "Please run as root (with sudo command)"
  exit
fi

if ! command -v wget &> /dev/null
then
    echo "wget command could not be found. installing..."
    apt install -y wget
fi

if ! command -v make &> /dev/null
then
    echo "make command could not be found. installing..."
    apt install -y cmake
fi

if ! command -v g++ &> /dev/null
then
    echo "g++ command could not be found. installing..."
    apt install -y g++
fi

echo "installing required dependancies..."
sudo apt install -y libgflags-dev libsnappy-dev zlib1g-dev libbz2-dev liblz4-dev libzstd-dev

echo "downloading rocksdb..."
pushd /tmp || return
wget https://github.com/facebook/rocksdb/archive/refs/tags/v${ROCKSDB_VERSION}.tar.gz
tar -xvf v${ROCKSDB_VERSION}.tar.gz && cd rocksdb-${ROCKSDB_VERSION} || return

# Ignore GCC warnings
export CXXFLAGS='-Wno-error=deprecated-copy -Wno-error=pessimizing-move -Wno-error=class-memaccess'

# Build as a shared library
make shared_lib

# The following command installs the shared library in /usr/lib/ and the header files in /usr/include/rocksdb/:
make install-shared INSTALL_PATH=/usr
popd || return

# cleanup
rm -rf /tmp/rocksdb-${ROCKSDB_VERSION} /tmp/v${ROCKSDB_VERSION}.tar.gz

echo "installation successful"

# To uninstall
#make uninstall INSTALL_PATH=/usr