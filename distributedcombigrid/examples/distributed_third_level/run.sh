#!/bin/bash

#PBS -N test_withstats
#PBS -l select=1:node_type=rome:mpiprocs=128
#PBS -l walltime=00:40:00


SGPP_DIR=/lustre/hpe/ws10/ws10.1/ws/ipvpolli-widely/DisCoTec/
LIB_BOOST_DIR=
LIB_GLPK=$SGPP_DIR/glpk/lib/

export LD_LIBRARY_PATH=$SGPP_DIR/lib/sgpp:$LIB_GLPK:$LIB_BOOST_DIR:$LD_LIBRARY_PATH

module switch mpt/2.23 openmpi/4.0.5
module load boost/1.70.0 # for boost

# Change to the direcotry that the job was submitted from
cd $PBS_O_WORKDIR

paramfile="ctparam_smaller_tl_system0"
# allows to read the parameter file from the arguments.
# Useful for testing the third level combination on a single system
if [ $# -ge 1 ] ; then
   paramfile=$1
fi

ngroup=$(grep ngroup $paramfile | awk -F"=" '{print $2}')
nprocs=$(grep nprocs $paramfile | awk -F"=" '{print $2}')

mpiprocs=$((ngroup*nprocs+1))


# General
#mpirun -n "$mpiprocs" omplace -v -ht spread -c 0-127:bs=128+st=128 ./xthi
mpirun -n "$mpiprocs" ./combi_example $paramfile
#mpirun -n "$mpiprocs" --bind-to core ./combi_example $paramfile
#mpirun -n "$mpiprocs" omplace -v -ht spread -c 0-127:bs=128+st=128 ./combi_example $paramfile
# Use for debugging
#mpirun -n "$mpiprocs" xterm -hold -e gdb -ex run --args ./combi_example $paramfile

