#!/bin/bash
### batch jobs run in /bin/bash per default
### MAXWT 24
### PROCSPERNODE 8
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
#$ -pe impi4 NMPIPROCS
###
### wallclock limit
#$ -l h_rt=WALLTIME
###
### the job's name
#$ -N JOBNAME

#load environment
module load impi/4.0.3 intel mkl fftw petsc slepc
module load hdf5-mpi

### set OMP_NUM_THREADS, switch I_MPI_PIN_DOMAIN to omp
### to run with OpenMP and use the -perhost option in the
### mpiexec call
### (see also /afs/ipp/.cs/intel/impi/4.0.0/doc/Reference_Manual.pdf;
### especially chapter 3.2)
export OMP_NUM_THREADS=1
#export I_MPI_PIN_DOMAIN=omp


### set MKL to serial mode
export MKL_SERIAL=yes

### enable/disable Shared Receive Queue
### (currently we observe sporadical hang-ups if enabled)
export MV2_USE_SRQ=0

### start program
### Note: when running with OMP_NUM_THREADS>1, you may need to replace
### NMPIPROCS by the actual number of mpi processes ... 

mpiexec -l -n NMPIPROCS ./gene_gp

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np NMPIPROCS --ppn 8
