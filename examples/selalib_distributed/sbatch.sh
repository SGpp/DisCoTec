#!/bin/bash
# Job Name and Files (also --job-name)
#SBATCH -J bsl_highres_combi
#Output and error (also --output, --error):
#SBATCH -o ./%x.%j.out
#SBATCH -e ./%x.%j.err
#Initial working directory (also --chdir):
#SBATCH -D ./
# Wall clock limit:
#SBATCH --time=72:00:00
#SBATCH --no-requeue

#SBATCH --ntasks=9

NGROUP=$(grep ngroup ctparam | awk -F"=" '{print $2}')
NPROCS=$(grep nprocs ctparam | awk -F"=" '{print $2}')

export OMP_NUM_THREADS=1;export OMP_PLACES=cores;export OMP_PROC_BIND=close

. ~/epyc/spack-newpackage/share/spack/setup-env.sh
spack load hdf5 /3f3o2w6

mpiexec.openmpi --rank-by core --map-by node:PE=${OMP_NUM_THREADS} -n $(($NGROUP*$NPROCS)) ./selalib_distributed_workers_only
#mpiexec.openmpi --rank-by core --map-by node:PE=${OMP_NUM_THREADS} -n $(($NGROUP*$NPROCS+1)) ./selalib_distributed
