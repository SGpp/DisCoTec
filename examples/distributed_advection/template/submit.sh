#!/bin/bash
# Job Name and Files (also --job-name)
#SBATCH -J advect
#Output and error (also --output, --error):
#SBATCH -o ./%x.%j.out
#SBATCH -e ./%x.%j.err
#Initial working directory (also --chdir):
#SBATCH -D ./
#Notification and type
#SBATCH --mail-type=END
# Wall clock limit:
#SBATCH --time=72:00:00
#SBATCH --no-requeue

#SBATCH --ntasks=9

mpiexec.openmpi -n $SLURM_NTASKS ./advection_example
