#############################################################
#$ -S /bin/bash
#$ -j y
#$ -cwd
#$ -m n
#$ -M $USER@ipp.mpg.de 
#$ -N GENE_bench
#$ -pe impi_hydra_cuda 32
#$ -l h_rt=00:30:00

export MACHINE=odin

module load impi
module load cuda
module load intel
module load perflib/2.1
export PERFLIB_OUTPUT_FORMAT=xml
export I_MPI_PIN_PROCESSOR_LIST=allcores

### set OMP_NUM_THREADS
export OMP_NUM_THREADS=1

### set MKL to serial mode
export MKL_SERIAL=yes

### set MPI algorithms for Intel MPI ###
#export I_MPI_ADJUST_REDUCE=4

### start program
mpiexec -np $NSLOTS ./gene_odin

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np $NSLOTS --mps 4
