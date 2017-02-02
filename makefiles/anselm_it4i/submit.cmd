#!/bin/bash
#PBS -q qprod
#PBS -N GENE
#PBS -l select=1:ncpus=16
##PBS -A OPEN-0-0

module load PrgEnv-intel fftw3

export OMP_NUM_THREADS=1
export MACHINE=anselm_it4i

cd $PBS_O_WORKDIR

mpirun -np 16 -hostfile $PBS_NODEFILE ./gene_anselm_it4i


