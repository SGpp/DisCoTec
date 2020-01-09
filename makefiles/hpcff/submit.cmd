#!/bin/bash -x
### nodes=number of nodes, ppn=processes per node (max. 8)
#MSUB -l nodes=1:ppn=8
### maximum walltime currently 12hrs
#MSUB -l walltime=4:00:00
### jobname
#MSUB -N GENE
###MSUB -M mail@mail.mail     official mail address
#MSUB -m n           
###sendmail: (n)ever, on (a)bort, (b)eginning or (e)nd of job 
###
###Number of OpenMP processes:
#MSUB -v tpt=1
 

cd $PBS_O_WORKDIR
export MKL_SERIAL=YES
export OMP_NUM_THREADS=1

### in the next line, the number of processes is calculated 
NP=$(cat $PBS_NODEFILE | wc -l)

mpiexec -env MKL_SERIAL yes -env OMP_NUM_THREADS 1 -np=$NP ./gene_hpcff

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np $NP --ppn 8 --mps 4
