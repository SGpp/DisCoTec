## MAXWT 24
## PROCSPERNODE 24
## SUBMITCMD qsub
#PBS -N JOBNAME
#PBS -q regular
#PBS -l mppwidth=NMPIPROCS
#PBS -l walltime=WALLTIME
#PBS -j eo
#PBS -V
#PBS -m abe

cd $PBS_O_WORKDIR

export NPROC=`qstat -f $PBS_JOBID | awk '/mppwidth/ {print $3}'`
export NTASK=`qstat -f $PBS_JOBID | awk '/mppnppn/  {print $3}'`
export OMP_NUM_THREADS=1

# Launch the parallel job to the allocated compute nodes
aprun -n $NPROC ./gene_edison