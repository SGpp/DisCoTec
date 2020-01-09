#The job name:
#PBS -N GENE

# Currently 24 hours is the maximum wall-clock limit.
#PBS -l walltime=02:00:00

# Remember that Lindgren uses an architecture with 24 cores per node; please
# ask for full nodes!  In order to avoid complications with auto-
# parallelisation, you may want to limit yourself to allocating a number of
# cores conforming to n_cores==24*(2**N) for some integer N. 
# The number of cores to be allocated is set below:
#PBS -l mppwidth=24

#PBS -e error_file.e
#PBS -o output_file.o

# Change to the work directory
cd $PBS_O_WORKDIR

# Add module support and modules:
. /opt/modules/default/etc/modules.sh
module unload PrgEnv-pgi petsc-complex slepc-complex                    #
module load PrgEnv-intel petsc-complex/3.3.00 
module add subversion fftw xt-libsci

# Find the required number of cores automatically:
export NPROC=`qstat -f $PBS_JOBID | awk '/mppwidth/ {print $3}'`

# If you have told GENE to run parallel GENE tasks, uncomment the line
# below and add  -N $NTASK  to the aprun command!
#export NTASK=`qstat -f $PBS_JOBID | awk '/mppnppn/  {print $3}'`

# Run GENE:
aprun -n $NPROC ./gene_lindgren > gene.out 2>&1

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np $NPROC --ppn 24 --mps 4
