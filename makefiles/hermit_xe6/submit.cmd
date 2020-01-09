#!/bin/bash
#PBS -N GENE
##Total number of procs:
#PBS -l mppwidth=32
##Number of procs per node (32 max):
#PBS -l mppnppn=32
#PBS -l walltime=00:30:00             
  
# Change to the direcotry that the job was submitted from
cd $PBS_O_WORKDIR

export NPROC=`qstat -f $PBS_JOBID | awk '/mppwidth/ {print $3}'`
export NTASK=`qstat -f $PBS_JOBID | awk '/mppnppn/  {print $3}'`
export OMP_NUM_THREADS=1

# Launch the parallel job to the allocated compute nodes
aprun -n $NPROC -N $NTASK ./gene_hermit_xe6

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np $NPROC --ppn $NTASK --mps 4
