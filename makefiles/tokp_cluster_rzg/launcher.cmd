#!/bin/bash
### batch jobs run in /bin/bash per default
### MAXWT 24
### PROCSPERNODE 16
### SUBMITCMD qsub
### 
###
### check 'man qsub' for all available options
###
### job has the same environmnent variables as 
### the submission shell (currently disabled)
###$ -V
### join stdout and stderr
#$ -j y
### start the job in the current working dir
#$ -cwd
### send notification (on (a)bort, (b)egin, (e)nd, (n)ever)
#$ -m n
###
### start a parallel environment with multiples of 8 processes
#$ -pe impi_hydra NMPIPROCS
###
### wallclock limit
#$ -l h_rt=WALLTIME
###
### the job's name
#$ -N JOBNAME

#load environment
module load impi intel mkl fftw petsc slepc
module load hdf5-mpi

### set OMP_NUM_THREADS
export OMP_NUM_THREADS=$THREADS_PER_MPI_TASK

### set MKL to serial mode 
export MKL_SERIAL=yes

### start program
### Note: when running with OMP_NUM_THREADS>1, you may need to replace
### $NSLOTS by the actual number of mpi processes and set 
### -perhost <<MPI-PROCS-PER-HOST>> flag

mpiexec -l -n NMPIPROCS ./gene_tokp_cluster_rzg
