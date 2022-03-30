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

export LD_LIBRARY_PATH=$(pwd)/../../../lib/sgpp:$(pwd)/../../../glpk/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=1

. ~/spack/share/spack/setup-env.sh
spack load boost@1.74.0
spack load hdf5@1.10.5

mpiexec.openmpi -n $SLURM_NTASKS ./selalib_distributed
