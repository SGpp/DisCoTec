#############################################################
#$ -S /bin/bash
#$ -j y
#$ -cwd
#$ -m n
#$ -M $USER@ipp.mpg.de 
#$ -N GENE_bench
#$ -pe impi 2048
#$ -l h_rt=04:00:00
#$ -R yes

export MACHINE=aims_cluster

module load impi
#export PATH=/u/rzgbench/mv2-1.4-intel/bin:$PATH
#NP=$(cat $TMPDIR/machines |wc -l)
#cat $TMPDIR/machines


### set OMP_NUM_THREADS
export OMP_NUM_THREADS=1

### set MKL to serial mode
export MKL_SERIAL=yes

### set MPI algorithms for Intel MPI ###
export I_MPI_ADJUST_REDUCE=4

### start program
#mpiexec -n 8 ./gene_aims_cluster
mpiexec -np $NSLOTS ./gene_aims_cluster

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np $NSLOTS --mps 4
