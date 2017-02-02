### MAXWT 24
### PROCSPERNODE 24
### SUBMITCMD qsub
#PBS -N JOBNAME

# Currently 24 hours is the maximum wall-clock limit.
#PBS -l walltime=WALLTIME

# Number of cores to be allocated is NMPIPROCS. Remember that Lindgren uses
# an architecture with 24 cores per node; please ask for full nodes to 
# avoid complications!
#PBS -l mppwidth=NMPIPROCS

#PBS -e error_file.e
#PBS -o output_file.o

# Change to the work directory
cd $PBS_O_WORKDIR

# Add modules:
. /opt/modules/default/etc/modules.sh
module unload PrgEnv-pgi petsc-complex slepc-complex                    #
module load PrgEnv-intel petsc-complex/3.3.00 
module add subversion fftw xt-libsci

# Run the executable named gene 
# and write the output into gene.out
aprun -n NMPIPROCS ./gene_lindgren > gene.out 2>&1
