#!/bin/bash
### MAXWT 24
### SUBMITCMD llsubmit
### PROCSPERNODE 16
#@ wall_clock_limit = WALLTIME
#@ job_type = parallel
#@ job_name = JOBNAME
#@ class = general
#@ network.MPI = sn_all,not_shared,us
#@ node = NODES
#@ tasks_per_node = PPN
#@ output = $(job_name).$(jobid).out
#@ error = $(job_name).$(jobid).err
#@ queue

. /etc/profile
. /etc/profile.d/modules.sh

#setup of environment
#module unload mpi.ibm
#module load mpi.intel
module unload fftw hdf5
module load fftw/mpi hdf5/mpi
module load petsc/3.3c slepc/3.3

export MKL_SERIAL=YES
export OMP_NUM_THREADS=1

mpiexec -n NMPIPROCS ./gene_supermuc
