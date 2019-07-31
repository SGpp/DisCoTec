#!/bin/bash
LIB_SGPP_DIR=/home/marci/UNI/HIWI/combi/
LIB_SIMPLE_AMQP_DIR=$LIB_SGPP_DIR/distributedcombigrid/SimpleAmqpClient/build
export LD_LIBRARY_PATH=$LIB_SGPP_DIR/lib/sgpp:$LIB_SIMPLE_AMQP_DIR:$LD_LIBRARY_PATH
MPI_PATH=/opt/mpich/bin/


paramfile="ctparam"
# allows to read the parameter file from the arguments.
# Useful for testing the third level combination on a single system
if [ $# -ge 1 ] ; then
   paramfile=$1
fi

ngroup=$(grep ngroup $paramfile | awk -F"=" '{print $2}')
nprocs=$(grep nprocs $paramfile | awk -F"=" '{print $2}')

mpiprocs=$((ngroup*nprocs+1))

$MPI_PATH/mpirun -n "$mpiprocs" ./combi_example $paramfile

# Use for debugging
#$MPI_PATH/mpirun -n "$mpiprocs" xterm -hold -e gdb -ex run --args ./combi_example $paramfile
