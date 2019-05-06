#!/bin/bash

MPI_PATH=/opt/mpich/bin/
LIB_SGPP_DIR=../../
LIB_SIMPLE_AMQP_DIR=$LIB_SGPP_DIR/distributedcombigrid/SimpleAmqpClient/build
export LD_LIBRARY_PATH=$LIB_SGPP_DIR/lib/sgpp:$LIB_SIMPLE_AMQP_DIR:$LD_LIBRARY_PATH

#$MPI_PATH/mpirun -n 6 -l xterm -hold -e ./test_distributedcombigrid_boost --run_test=thirdLevel -- # additional params
$MPI_PATH/mpirun -n 6 -l ./test_distributedcombigrid_boost --run_test=thirdLevel -- # additional params
