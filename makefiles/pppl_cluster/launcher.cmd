#!/bin/tcsh
### PROCSPERNODE 8
### MAXWT 24
### SUBMITCMD qsub
#
## --- run the job on 'n' nodes with 'y' processors per node (ppn)
#PBS -l nodes=NODES:ppn=PPN
#
# --- specify the period of time the job will run in the format
#PBS -l walltime=WALLTIME
#
#PBS -N JOBNAME
#
##PBS -M $USER@pppl.gov
#PBS -m n
#
##PBS -l mem=2000mb
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
