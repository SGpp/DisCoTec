#!/bin/bash

MPI_PATH=/opt/mpich/bin/
LIB_SGPP_DIR=../../
export LD_LIBRARY_PATH=$LIB_SGPP_DIR/lib/sgpp:$LD_LIBRARY_PATH

maxNgroup=2
maxNprocs=2
total=$(( 2*(maxNgroup * maxNprocs + 1) ))

numCombinations=10

PORT=9999
PORT_USAGE=$(lsof -P -i tcp:$PORT)

if [ "$PORT_USAGE" == "" ]; then
  $MPI_PATH/mpirun -n $total -l ./test_distributedcombigrid_boost --run_test=thirdLevel -- # additional params
else
  echo -e "Port $PORT is already in use: \n$PORT_USAGE"
fi


# Use this to simplify Debug:
#$MPI_PATH/mpirun -n $total -l xterm -hold -e gdb -ex run --args ./test_distributedcombigrid_boost --run_test=thirdLevel
