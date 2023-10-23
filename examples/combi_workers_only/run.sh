#!/bin/bash
#SBATCH --exclusive
#SBATCH --time=00:10:00

NGROUP=$(grep ngroup ctparam | awk -F"=" '{print $2}')
NPROCS=$(grep nprocs ctparam | awk -F"=" '{print $2}')

echo $(git log -1 --pretty=format:"%H")

cat ctparam

export OMP_NUM_THREADS=4;export OMP_PLACES=cores;export OMP_PROC_BIND=close

# heaptrack call; use to make sure no unexpected allocations are made
# can be obtained with `wget https://download.kde.org/stable/heaptrack/1.4.0/heaptrack-v1.4.0-x86_64.AppImage`
# mpiexec --rank-by core --map-by node:PE=${OMP_NUM_THREADS} -n $(($NGROUP*$NPROCS)) ./heaptrack-v1.4.0-x86_64.AppImage ./combi_example_worker_only ctparam

mpiexec --rank-by core --map-by node:PE=${OMP_NUM_THREADS} --report-bindings -n $(($NGROUP*$NPROCS)) ./combi_example_worker_only ctparam
