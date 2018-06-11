#!/bin/bash
export LD_LIBRARY_PATH=/home/shreyas/Thesis/combi/lib/sgpp:$LD_LIBRARY_PATH

NGROUP=$(grep ngroup ctparam | awk -F"=" '{print $2}')
NPROCS=$(grep nprocs ctparam | awk -F"=" '{print $2}')

mpirun -n $(($NGROUP*$NPROCS+1)) ./combi_example
