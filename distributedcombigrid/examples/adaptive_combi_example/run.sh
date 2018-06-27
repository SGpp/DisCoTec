#!/bin/bash
export LD_LIBRARY_PATH=/home/simon/Uni/master/idp/comb/lib/sgpp:$LD_LIBRARY_PATH

NGROUP=$(grep ngroup ctparam | awk -F"=" '{print $2}')
NPROCS=$(grep nprocs ctparam | awk -F"=" '{print $2}')

mpirun -n $(($NGROUP*$NPROCS+1)) --oversubscribe ./combi_example
#mpirun -n $(($NGROUP*$NPROCS+1)) --oversubscribe urxvt -e gdb ./combi_example
