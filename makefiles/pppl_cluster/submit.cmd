#!/bin/tcsh
#
#PBS -N STDTST
#
#PBS -m n
#
##PBS -M $USER@pppl.gov
#
## --- run the job on 'n' nodes with 'y' processors per node (ppn)
#PBS -l nodes=1:ppn=8
#
#--- memory required for job
##PBS -l mem=2000mb
#
# --- specify the period of time the job will run in the format
#PBS -l walltime=01:00:00
#
# --- do not rerun this job if it fails
#PBS -r n
#
# export all my environment variables to the job
##PBS -V
#
#PBS -j oe

# --- $NPROCS is required by mpirun.
set NPROCS=`wc -l < $PBS_NODEFILE`

cd $PBS_O_WORKDIR

setenv OMP_NUM_THREADS 1

module load pathscale openmpi
module load fftw/3.1.2 lapack-blas
module load acml/4.4.0-open64-64bit
module load scalapack/1.8.0
module load petsc_complex/3.2-p6
module load slepc_complex/3.2-p3

mpiexec -n $NPROCS ./gene_pppl_cluster


### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np 8 --ppn 8 --mps 4
