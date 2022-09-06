#!/bin/bash
# Job Name and Files (also --job-name)
#SBATCH -J tests
#Output and error (also --output, --error):
#SBATCH -o ./%x.%j.out
#SBATCH -e ./%x.%j.err
#Initial working directory (also --chdir):
#SBATCH -D ./
# Wall clock limit:
#SBATCH --time=00:20:00
#SBATCH --no-requeue
#SBATCH --partition=micro
#Setup of execution environment
#SBATCH --export=NONE 
#SBATCH --account=pn34mi

#Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48

##fixed frequency, no dynamic adjustment
###SBATCH --ear=off
#optional: keep job within one island
#SBATCH --switches=1

. ./../../setenv.sh

mpiexec -np 9 ./test_distributedcombigrid_boost --log_level=message --run_test=ftolerance
mpiexec -np 9 ./test_distributedcombigrid_boost --log_level=message --run_test=fullgrid
mpiexec -np 9 ./test_distributedcombigrid_boost --log_level=message --run_test=hierarchization
mpiexec -np 9 ./test_distributedcombigrid_boost --log_level=message --run_test=loadmodel
mpiexec -np 9 ./test_distributedcombigrid_boost --log_level=message --run_test=reduce
mpiexec -np 9 ./test_distributedcombigrid_boost --log_level=message --run_test=rescheduling
mpiexec -np 9 ./test_distributedcombigrid_boost --log_level=message --run_test=stats
mpiexec -np 9 ./test_distributedcombigrid_boost --log_level=message --run_test=task
mpiexec -np 9 ./test_distributedcombigrid_boost --log_level=message --run_test=distributedfullgrid
mpiexec -np 9 ./test_distributedcombigrid_boost --log_level=message --run_test=integration
mpiexec -np 9 ./test_distributedcombigrid_boost --log_level=message --run_test=distributedsparsegrid
mpiexec -np 9 ./test_distributedcombigrid_boost --log_level=message --run_test=thirdLevel
