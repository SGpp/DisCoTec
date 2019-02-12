#!/bin/bash
MPI_DIR=/opt/mpich/bin/
LIB_SGPP_DIR=../../../
LIB_SIMPLE_AMQP_DIR=$LIB_SGPP_DIR/distributedcombigrid/SimpleAmqpClient/build

export LD_LIBRARY_PATH=$LIB_SGPP_DIR/lib/sgpp:$LIB_SIMPLE_AMQP_DIR:$LD_LIBRARY_PATH

# TODO 
#
NGROUP=$(grep ngroup ctparam | awk -F"=" '{print $2}')
NPROCS=$(grep nprocs ctparam | awk -F"=" '{print $2}')

$MPI_DIR/mpirun -n $(($NGROUP*$NPROCS+1)) ./combi_example $1
