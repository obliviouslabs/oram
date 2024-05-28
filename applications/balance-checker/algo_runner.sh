#!/bin/bash
source /startsgxenv.sh
export LD_LIBRARY_PATH=/usr/local/lib/:$LD_LIBRARY_PATH

cp ../../omap/common/encutils.cpp ./Enclave/TrustedLibrary/encutils.cpp


SGX_MODE=HW # HW or SIM
# g++ init_db_balance.cpp -o init_db_balance -L/usr/local/lib -lrocksdb
# ./init_db_balance ./db ./rcc_balance.txt 
# ./init_db_balance ./db_usdt ./usdt_balance.txt
# Algorithms:
MIN_ENCLAVE_SIZE=16384 # enclave size in MB
MAX_ENCLAVE_SIZE=16384
TCS_NUM=2
CORE_ID=0 # the cpu core id to run the program
DISK_IO=0 # 0: no disk IO, 1: disk IO
DB_PATH=./db_usdt

if [ $SGX_MODE = HW ]
then
    SGX_PRERELEASE=1
else
    SGX_PRERELEASE=0
fi

if [ $DISK_IO = 1 ]
then
    DISK_TAG=_DISK
fi

if [ $MAX_ENCLAVE_SIZE != 128 ]
then
    ENCLAVE_SIZE_TAG=_${MIN_ENCLAVE_SIZE}_${MAX_ENCLAVE_SIZE}
fi

FILENAME=OMap${DISK_TAG}${ENCLAVE_SIZE_TAG}.out
if [ -z "$1" ]; then
rm -f $FILENAME
echo "output to "${FILENAME}
fi
for (( encsize=$MIN_ENCLAVE_SIZE; encsize<=$MAX_ENCLAVE_SIZE; encsize*=2 ))
do
heapsizeB=$(( encsize * 1000000 ))
hex_encsize=$(printf '%x\n' $heapsizeB)

sed -i "/.*<Heap.*/c\  <HeapMaxSize>0x"${hex_encsize}"</HeapMaxSize>" ./Enclave/Enclave.config.xml

make clean
make SGX_MODE=$SGX_MODE SGX_PRERELEASE=$SGX_PRERELEASE DISK_IO=$DISK_IO ENCLAVE_SIZE=$encsize
if [[ $1 = 1 ]]; then
    taskset -c ${CORE_ID} ./omap.elf $DB_PATH
    sleep 1
else
    taskset -c ${CORE_ID} stdbuf -oL nohup ./omap.elf $DB_PATH &>> $FILENAME < /dev/null
fi
done

sed -i '/.*<Heap.*/c\  <HeapMaxSize>0x7A00000</HeapMaxSize>' ./Enclave/Enclave.config.xml