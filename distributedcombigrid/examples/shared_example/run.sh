#!/bin/bash
export LD_LIBRARY_PATH=../../../lib/sgpp/:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=$1

NGROUP=$(grep ngroup ctparam | awk -F"=" '{print $2}')
NPROCS=$(grep nprocs ctparam | awk -F"=" '{print $2}')

mpirun.mpich -n $(($NGROUP*$NPROCS+1)) ./shared_example $2
