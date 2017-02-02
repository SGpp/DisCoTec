### MAXWT 24
### PROCSPERNODE 32
### SUBMITCMD qsub

#!/bin/bash
#PBS -N JOBNAME
##Total number of procs:
#PBS -l mppwidth=NMPIPROCS
##Number of procs per node (32 max):
#PBS -l mppnppn=32
#PBS -l walltime=WALLTIME
  
# Change to the direcotry that the job was submitted from
cd $PBS_O_WORKDIR

export NPROC=`qstat -f $PBS_JOBID | awk '/mppwidth/ {print $3}'`
export NTASK=`qstat -f $PBS_JOBID | awk '/mppnppn/  {print $3}'`
export OMP_NUM_THREADS=1

# Launch the parallel job to the allocated compute nodes
aprun -n $NPROC -N $NTASK ./gene_hermit_xe6