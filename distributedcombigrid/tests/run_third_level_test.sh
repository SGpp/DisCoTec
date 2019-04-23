#!/bin/bash

LIB_SGPP_DIR=../../
LIB_SIMPLE_AMQP_DIR=$LIB_SGPP_DIR/distributedcombigrid/SimpleAmqpClient/build
export LD_LIBRARY_PATH=$LIB_SGPP_DIR/lib/sgpp:$LIB_SIMPLE_AMQP_DIR:$LD_LIBRARY_PATH

mpirun -n 3 -l xterm -hold -e gdb -ex run --args ./test_distributedcombigrid_boost --run_test=thirdLevel -- # additional params
#mpirun -n 9 -l ./test_distributedcombigrid_boost --run_test=reduce -- # additional params
