#!/bin/bash
### batch jobs run in /bin/bash per default
### MAXWT 24
### PROCSPERNODE 16
### SUBMITCMD qsub
### 
###
### check 'man qsub' for all available options
###
#PBS -q qprod
#PBS -N JOBNAME
#PBS -l select=NODES:ncpus=PPN
##PBS -A OPEN-0-0

module load PrgEnv-intel fftw3

export OMP_NUM_THREADS=1
export MACHINE=anselm_it4i

cd $PBS_O_WORKDIR

mpirun -np NMPIPROCS -hostfile $PBS_NODEFILE ./gene_anselm_it4i
