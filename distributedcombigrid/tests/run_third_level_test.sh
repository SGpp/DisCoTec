#!/bin/bash

MPI_PATH=/opt/mpich/bin/
SGPP_DIR=../../
LIB_GLPK_DIR=$SGPP_DIR/glpk/lib/

export LD_LIBRARY_PATH=$SGPP_DIR/lib/sgpp:$LIB_GLPK_DIR:$LD_LIBRARY_PATH

# set these to the maximum of all third_level tests, such that each test has
# enough procs available
maxNgroup=2
maxNprocs=4
total=$(( 2*(maxNgroup * maxNprocs + 1) ))

numCombinations=10

PORT=9999
PORT_USAGE=$(lsof -P -i tcp:$PORT)

if [ "$PORT_USAGE" == "" ]; then
  $MPI_PATH/mpirun -n $total -l ./test_distributedcombigrid_boost --run_test=thirdLevel -- # additional params

  # Use this to simplify Debug:
  #$MPI_PATH/mpirun -n $total -l xterm -hold -e gdb -ex run --args ./test_distributedcombigrid_boost --run_test=thirdLevel
else
  echo -e "Port $PORT is already in use: \n$PORT_USAGE"
fi


