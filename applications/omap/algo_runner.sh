#!/bin/bash
source /startsgxenv.sh

cp ../../omap/common/encutils.cpp ./Enclave/TrustedLibrary/encutils.cpp

SGX_MODE=SIM # HW or SIM

# Algorithms:
MIN_ENCLAVE_SIZE=8192 # enclave size in MB
MAX_ENCLAVE_SIZE=8192
CORE_ID=5 # the cpu core id to run the program
DISK_IO=0 # 0: no disk IO, 1: disk IO
TCS_NUM=2

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


FILENAME=OMap${IO_TAG}${DISK_TAG}${ENCLAVE_SIZE_TAG}.out
if [ -z "$1" ]; then
rm -f $FILENAME
echo "output to "${FILENAME}
fi
for (( encsize=$MIN_ENCLAVE_SIZE; encsize<=$MAX_ENCLAVE_SIZE; encsize*=2 ))
do
heapsizeB=$(( encsize * 1048576 ))
hex_encsize=$(printf '%x\n' $heapsizeB)

sed -i "/.*<Heap.*/c\  <HeapMaxSize>0x"${hex_encsize}"</HeapMaxSize>" ./Enclave/Enclave.config.xml

sed -i "/.*<TCSNum.*/c\  <TCSNum>"${TCS_NUM}"</TCSNum>" ./Enclave/Enclave.config.xml

make clean
make SGX_MODE=$SGX_MODE SGX_PRERELEASE=$SGX_PRERELEASE DISK_IO=$DISK_IO ENCLAVE_SIZE=$encsize TCS_NUM=$TCS_NUM
if [[ $1 = 1 ]]; then
    ./omap.elf
    sleep 1
else
    stdbuf -oL nohup ./omap.elf &>> $FILENAME < /dev/null
fi
done

sed -i '/.*<Heap.*/c\  <HeapMaxSize>0x7A00000</HeapMaxSize>' ./Enclave/Enclave.config.xml
sed -i "/.*<TCSNum.*/c\  <TCSNum>1</TCSNum>" ./Enclave/Enclave.config.xml