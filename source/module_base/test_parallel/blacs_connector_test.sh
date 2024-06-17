#!/bin/bash -e

np=`cat /proc/cpuinfo | grep "cpu cores" | uniq| awk '{print $NF}'`
echo "nprocs in this machine is $np"

for i in {8..2};
do
    if [[ $i -gt $np ]];then
        continue
    fi
    echo "TEST in parallel, nprocs=$i"
    mpirun -np $i ./blacs_connector
    break    
done

