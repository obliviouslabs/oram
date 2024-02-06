#!/bin/bash
source /startsgxenv.sh
export LD_LIBRARY_PATH=/usr/local/lib/:$LD_LIBRARY_PATH
SGX_MODE=HW # HW or SIM
# g++ init_db_balance.cpp -o init_db_balance -L/usr/local/lib -lrocksdb
# ./init_db_balance ./db ./rcc_balance.txt 
# ./init_db_balance ./db_usdt ./usdt_balance.txt
# Algorithms:
MIN_ENCLAVE_SIZE=4096 # enclave size in MB
MAX_ENCLAVE_SIZE=4096
IO_ROUNDs=(1) # number of rounds encryption/decryption is performed, used to get breakdown
CORE_ID=0 # the cpu core id to run the program
DISK_IO=0 # 0: no disk IO, 1: disk IO

for IO_ROUND in ${IO_ROUNDs[@]}; do
if [ $IO_ROUND = 0 ]
then
    IO_TAG=_MOCK_IO
fi

if [ $IO_ROUND -gt 1 ]
then
    IO_TAG=_${IO_ROUND}IO
fi

if [ $DISK_IO = 1 ]
then
    DISK_TAG=_DISK
fi

if [ $MAX_ENCLAVE_SIZE != 128 ]
then
    ENCLAVE_SIZE_TAG=_${MIN_ENCLAVE_SIZE}_${MAX_ENCLAVE_SIZE}
fi


FILENAME=OMap${IO_TAG}${DISK_TAG}${ENCLAVE_SIZE_TAG}.out
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
make SGX_MODE=$SGX_MODE SGX_PRERELEASE=1 IO_ROUND=$IO_ROUND DISK_IO=$DISK_IO ENCLAVE_SIZE=$encsize
if [[ $1 = 1 ]]; then
    taskset -c ${CORE_ID} ./omap.elf
    sleep 1
else
    taskset -c ${CORE_ID} stdbuf -oL nohup ./omap.elf &>> $FILENAME < /dev/null
fi
done
done

sed -i '/.*<Heap.*/c\  <HeapMaxSize>0x7A00000</HeapMaxSize>' ./Enclave/Enclave.config.xml