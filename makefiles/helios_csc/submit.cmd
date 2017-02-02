#!/bin/sh
#SBATCH --ntasks=1024       ### Total number of MPI tasks (16 cores/node)
##SBATCH -N 64              ### Number of nodes; sometimes req'd         
#SBATCH --cpus-per-task=1   ### Number of threads per task (OpenMP)
#SBATCH --time=04:00:00     ### wall clock time
#SBATCH --job-name=GENE
##SBATCH -A <Myproject>     ### Uncomment for specific project budget 
                            ### * check with "id -nG" for your choices
			    ### * add "-0" to project name for free
			    ###   low priority queue

##uncomment next line and specify faulty node(s) XYZ to be excluded
##SBATCH --exclude=heliosXYZ

##list nodes for debugging
printenv SLURM_NODELIST

##uncomment next line for core files
#ulimit -c unlimited

##uncomment ulimit line to increase stack size 
##(default had been just 10MB in the past, should be 512MB now)
#ulimit -s <your-desired-stack-size>

##bullxmpi (if not already set in the login shell)
module load bullxmpi intel

##uncomment the following lines for intelmpi (if not already set in the login shell)
#module load intelmpi intel
#export I_MPI_EXTRA_FILESYSTEM=enable
#export I_MPI_EXTRA_FILESYSTEM_LIST=lustre
#mpdboot

##set openmp threads
export OMP_NUM_THREADS=1

##execute
mpiexec -n $SLURM_NTASKS ./gene_helios_csc

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np $SLURM_NTASKS --ppn 16 --mps 4
