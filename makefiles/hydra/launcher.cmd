# @ shell=/bin/tcsh
#
## MAXWT 24
## PROCSPERNODE 20
## SUBMITCMD llsubmit
# @ job_name   = JOBNAME
# @ error   = $(job_name).$(jobid).err
# @ output  = $(job_name).$(jobid).out
# @ job_type = parallel
# @ node_usage= not_shared
# @ node = NODES
# @ tasks_per_node = PPN
# @ resources = ConsumableCpus(1)
# @ network.MPI = sn_all,not_shared,us
# @ wall_clock_limit = WALLTIME
# @ notification = never
# @ queue

module load fftw petsc-cplx slepc-cplx hdf5-mpi
setenv OMP_NUM_THREADS 1

poe ./gene_hydra

