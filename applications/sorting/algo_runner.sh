#!/bin/bash
source /startsgxenv.sh

SGX_MODE=HW # HW or SIM

# Algorithms:
ALGOs=(CABUCKETSORT)
MIN_ELEMENT_SIZE=128 # element size in bytes
MAX_ELEMENT_SIZE=128
MIN_SIZE=1000000    # input size in number of elements
MAX_SIZE=10000000
MIN_ENCLAVE_SIZE=128 # enclave size in MB
MAX_ENCLAVE_SIZE=128
IO_ROUNDs=(1) # number of rounds encryption/decryption is performed, used to get breakdown
CORE_ID=5 # the cpu core id to run the program
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

for ALGO in ${ALGOs[@]}; do
FILENAME=${ALGO}_${MIN_ELEMENT_SIZE}_${MAX_ELEMENT_SIZE}_${MIN_SIZE}_${MAX_SIZE}${IO_TAG}${DISK_TAG}${ENCLAVE_SIZE_TAG}.out
if [ -z "$1" ]; then
rm -f $FILENAME
echo "output to "${FILENAME}
fi
for (( encsize=$MIN_ENCLAVE_SIZE; encsize<=$MAX_ENCLAVE_SIZE; encsize*=2 ))
do
heapsizeB=$(( encsize * 1000000 ))
hex_encsize=$(printf '%x\n' $heapsizeB)

sed -i "/.*<Heap.*/c\  <HeapMaxSize>0x"${hex_encsize}"</HeapMaxSize>" ./Enclave/Enclave.config.xml
for (( s=$MIN_ELEMENT_SIZE; s<=$MAX_ELEMENT_SIZE; s=s*3/2 ))
do
    make clean
    make SGX_MODE=$SGX_MODE SGX_PRERELEASE=1 ELEMENT_SIZE=$s ALGO=$ALGO MIN_SIZE=$MIN_SIZE MAX_SIZE=$MAX_SIZE IO_ROUND=$IO_ROUND DISK_IO=$DISK_IO ENCLAVE_SIZE=$encsize
    if [[ $1 = 1 ]]; then
        taskset -c ${CORE_ID} ./sort.elf
        sleep 1
    else
        taskset -c ${CORE_ID} stdbuf -oL nohup ./sort.elf &>> $FILENAME < /dev/null
        sleep 1
        last_line=$(tail -n 1 $FILENAME)

        last_token=$(echo -e "$last_line" | cut -f 3)

        # Check if duration is greater than 4000
        if (( $(echo "$last_token > 4000" | bc -l) )); then
            break
        fi
    fi
done
done
done
done
sed -i '/.*<Heap.*/c\  <HeapMaxSize>0x7A00000</HeapMaxSize>' ./Enclave/Enclave.config.xml