#maximal wallclock time (the default in the launcher)
### MAXWT 240
#maximal number of MPI processes per node
### PROCSPERNODE 12
#submit command
### SUBMITCMD qsub

# queue
#$ -q std.q
# join stdout and stderr
#$ -j y
# start the job in the current working dir
#$ -cwd
# send notification (on (a)bort, (b)egin, (e)nd, (n)ever)
#$ -m n

# fillup means that each host is filled up to its limit with processes
#$ -pe openmpi-fillup NMPIPROCS

# wallclock limit (inserted by the launcher)
#$ -l h_rt=WALLTIME
# the job's name (inserted by the launcher)
#$ -N JOBNAME

mpiexec -n NMPIPROCS ./gene_leo3_cluster_uibk
