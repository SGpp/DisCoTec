#!/bin/sh
## submit via:  bsub < submit.cmd
## check queue: bjobs -u all
## check job:   bpeek <job_id>
## kill job:    bkill <job_id>
#BSUB -a openmpi
#BSUB -q priority
#BSUB -n 8
#BSUB -o job.%J.out
#BSUB -e job.%J.out
#BSUB -J GENE
#BSUB -c 00:30

export OMP_NUM_THREADS=1
mpirun -np `cat $LSB_DJOB_HOSTFILE | wc -l` --hostfile $LSB_DJOB_HOSTFILE ./gene_ITM-Gateway


### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np 8 --ppn 1 --mps 4
