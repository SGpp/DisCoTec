#!/bin/bash
# DO NOT USE environment = COPY_ALL
#@ wall_clock_limit = 1:00:00
#@ job_type = parallel
###@ job_type = MPICH  for intel MPI
#@ job_name = GENE
#@ class = general
#@ network.MPI = sn_all,not_shared,us
#@ node = 1
#@ tasks_per_node = 8
#@ output = $(job_name).$(jobid).out
#@ error = $(job_name).$(jobid).err
#@ queue

. /etc/profile
. /etc/profile.d/modules.sh

export MKL_SERIAL=YES
export OMP_NUM_THREADS=1

#setup of environment
#module unload mpi.ibm
#module load mpi.intel
module unload fftw hdf5
module load fftw/mpi hdf5/mpi
module load petsc/3.3c slepc/3.3

mpiexec -n 8 ./gene_supermuc

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np 8 --ppn 8 --mps 4
