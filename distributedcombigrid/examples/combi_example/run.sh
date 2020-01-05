#!/bin/bash
export PATH_TO_SGPP=`pwd`/../../../
export LD_LIBRARY_PATH=$PATH_TO_SGPP/lib/sgpp:$LD_LIBRARY_PATH

NGROUP=$(grep ngroup ctparam | awk -F"=" '{print $2}')
NPROCS=$(grep nprocs ctparam | awk -F"=" '{print $2}')

mpirun.mpich -n $(($NGROUP*$NPROCS+1)) -l ./combi_example
