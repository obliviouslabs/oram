FROM ubuntu:22.04
ENV TZ=Europe/Lisbon
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt update -y && apt upgrade -y
RUN apt install -y wget gnupg

RUN echo 'deb [arch=amd64] https://download.01.org/intel-sgx/sgx_repo/ubuntu jammy main' | tee /etc/apt/sources.list.d/intel-sgx.list
RUN wget -qO - https://download.01.org/intel-sgx/sgx_repo/ubuntu/intel-sgx-deb.key | apt-key add -

RUN apt update -y && apt upgrade -y \
&& apt install -y build-essential lld make git wget unzip cmake ninja-build gdb \
        ocaml ocamlbuild automake autoconf libtool wget python-is-python3 libssl-dev git cmake perl \
        python3 python3-pip lcov gcovr cppcheck \
&& apt install -y clang clang-13 libc++-13-dev libc++abi-13-dev clang-format clang-format-13 \
        libx86-dev libclang-common-13-dev libclang-common-14-dev \
&& apt install -y libgtest-dev libboost-dev libssl-dev \
&& apt install -y doxygen \
&& apt install -y python3-matplotlib python3-numpy bc \
&& apt install -y libssl-dev libcurl4-openssl-dev protobuf-compiler libprotobuf-dev \
        debhelper cmake reprepro unzip pkgconf libboost-dev libboost-system-dev libboost-thread-dev \
        protobuf-c-compiler libprotobuf-c-dev lsb-release libsystemd0 \
build-essential lld make git wget unzip cmake ninja-build gdb \
        ocaml ocamlbuild automake autoconf libtool wget python-is-python3 libssl-dev git cmake perl \
        python3 python3-pip lcov gcovr cppcheck \
clang clang-13 libc++-13-dev libc++abi-13-dev clang-format clang-format-13 \
libx86-dev libclang-common-13-dev libclang-common-14-dev libgtest-dev libboost-dev libssl-dev  doxygen  python3-matplotlib python3-numpy  libssl-dev libcurl4-openssl-dev protobuf-compiler libprotobuf-dev \
debhelper cmake reprepro unzip pkgconf libboost-dev libboost-system-dev libboost-thread-dev \
protobuf-c-compiler libprotobuf-c-dev lsb-release libsystemd0 

RUN pip install cppcheck-codequality
RUN pip install requests
RUN pip install web3
RUN pip install tqdm

RUN apt-get install -y libsgx-urts libsgx-launch libsgx-enclave-common
RUN apt-get install -y libsgx-epid libsgx-quote-ex libsgx-dcap-ql libsgx-dcap-ql-dev


RUN wget 'https://github.com/boyter/scc/releases/download/v3.0.0/scc-3.0.0-x86_64-unknown-linux.zip'
RUN unzip 'scc-3.0.0-x86_64-unknown-linux.zip' -d /
RUN chmod +x /scc

RUN mkdir /dockerfiles
RUN cd /dockerfiles && wget https://download.01.org/intel-sgx/latest/linux-latest/distro/ubuntu22.04-server/sgx_linux_x64_sdk_2.19.100.3.bin \
&& chmod +x ./* \
&& echo yes | ./sgx_linux_x64_sdk_2.19.100.3.bin
RUN echo "source /dockerfiles/sgxsdk/environment" > /startsgxenv.sh
RUN chmod +x /startsgxenv.sh

COPY rocksDBInstall.sh /rocksDBInstall.sh
RUN chmod +x /rocksDBInstall.sh && \
    /rocksDBInstall.sh

SHELL ["/bin/bash", "-c"] 
RUN ln -s /usr/lib/libsgx_urts.so /usr/lib/libsgx_urts.so.2

RUN apt update -y && apt upgrade -y && apt install -y libbearssl-dev

WORKDIR /builder/
CMD ["bash"]
