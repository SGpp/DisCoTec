#!/bin/bash -x
### PROCSPERNODE 8
### MAXWT 24
### SUBMITCMD msub

### nodes=number of nodes, ppn=processes per node (max. 8)
#MSUB -l nodes=NODES:ppn=PPN
### maximum walltime currently 24hrs
#MSUB -l walltime=WALLTIME
### jobname
#MSUB -N JOBNAME
###MSUB -M mail@mail.mail     official mail address
#MSUB -m n           
###sendmail: (n)ever, on (a)bort, (b)eginning or (e)nd of job 
###
###Number of OpenMP processes:
#MSUB -v tpt=1
 

cd $PBS_O_WORKDIR
export MKL_SERIAL=YES
export OMP_NUM_THREADS=1

mpiexec -env MKL_SERIAL yes -env OMP_NUM_THREADS 1 -np=NMPIPROCS ./gene_hpcff
