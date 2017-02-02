#!/bin/sh
### SUBMITCMD bsub<
### PROCSPERNODE 8
### MAXWT 24

## submit via:  bsub < submit.cmd
## check queue: bjobs -u all
## check job:   bpeek <job_id>
## kill job:    bkill <job_id>
#BSUB -a openmpi
#BSUB -q priority
#BSUB -n NMPIPROCS
#BSUB -o job.%J.out
#BSUB -e job.%J.out
#BSUB -J JOBNAME
#BSUB -c WALLTIME

export OMP_NUM_THREADS=1
mpirun -np NMPIPROCS --hostfile $LSB_DJOB_HOSTFILE ./gene_ITM-Gateway
