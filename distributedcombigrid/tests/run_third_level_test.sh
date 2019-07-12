#!/bin/bash

MPI_PATH=/opt/mpich/bin/
LIB_SGPP_DIR=../../
LIB_SIMPLE_AMQP_DIR=$LIB_SGPP_DIR/distributedcombigrid/SimpleAmqpClient/build
export LD_LIBRARY_PATH=$LIB_SGPP_DIR/lib/sgpp:$LIB_SIMPLE_AMQP_DIR:$LD_LIBRARY_PATH

ngroup=2
nprocs=2
total=$(( 2*(ngroup * nprocs + 1) ))

# Use this to simplify Debug:
#$MPI_PATH/mpirun -n $total -l xterm -hold -e gdb -ex run --args ./test_distributedcombigrid_boost --run_test=thirdLevel -- $ngroup $nprocs # additional params

$MPI_PATH/mpirun -n $total -l ./test_distributedcombigrid_boost --run_test=thirdLevel -- $ngroup $nprocs # additional params
