#!/bin/bash
SGPP_DIR=../../../
LIB_BOOST_DIR=/usr/lib/

export LD_LIBRARY_PATH=$SGPP_DIR/lib/sgpp\
                       :$LIB_BOOST_DIR\
                       :$LD_LIBRARY_PATH

paramfile="ctparam"
# allows to read the parameter file from the arguments.
# Useful for testing the third level combination on a single system
if [ $# -ge 1 ] ; then
   paramfile=$1
fi

ngroup=$(grep ngroup $paramfile | awk -F"=" '{print $2}')
nprocs=$(grep nprocs $paramfile | awk -F"=" '{print $2}')

mpiprocs=$((ngroup*nprocs+1))


# On Helium
#mpiexec.mpich·-n·"$mpiprocs"·./combi_example·$paramfile
# Use for debugging
#mpiexec.mpich -n "$mpiprocs" xterm -hold -e gdb -ex run --args ./combi_example $paramfile

# General
MPI_PATH=/opt/mpich/bin/
$MPI_PATH/mpirun -n "$mpiprocs" ./combi_example $paramfile
# Use for debugging
#$MPI_PATH/mpirun -n "$mpiprocs" xterm -hold -e gdb -ex run --args ./combi_example $paramfile

# On HLRS
#aprun -n "$mpiprocs" ./combi_example $paramfile


