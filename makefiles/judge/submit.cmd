#!/bin/bash 
###MSUB -l nodes=1:ppn=8:gpus=2
#MSUB -l nodes=1:ppn=8
#MSUB -l walltime=1:00:00
#MSUB -N GENE
###MSUB -M mail@mail.mail     official mail address
#MSUB -m n           
###sendmail: (n)ever, on (a)bort, (b)eginning or (e)nd of job 
###
###Number of OpenMP processes:
#MSUB -v tpt=1
#MSUB -j oe

##############################################################
export OMP_NUM_THREADS=1
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
export MKL_SERIAL=YES

### in the next line, the number of processes is calculated 
NP=$(cat $PBS_NODEFILE | wc -l)
mpiexec -np=$NP --exports=MKL_SERIAL,OMP_NUM_THREADS,KMP_AFFINITY ./gene_judge
##############################################################
