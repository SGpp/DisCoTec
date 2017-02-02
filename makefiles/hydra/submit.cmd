# @ shell=/bin/tcsh
#
# @ job_name = GENE
# @ output   = $(job_name).$(jobid).out
# @ error    = $(job_name).$(jobid).err
# @ job_type = parallel
# @ node_usage= not_shared
# @ node = 1
## tasks_per_node = 20 (or 40 in SMT mode):
# @ tasks_per_node = 20
# @ resources = ConsumableCpus(1)
# @ network.MPI = sn_all,not_shared,us
# @ wall_clock_limit = 24:00:00
# @ notification = never
# @ queue

module load fftw petsc-cplx slepc-cplx hdf5-mpi

setenv OMP_NUM_THREADS 1
set NTASKS=`wc -l < $LOADL_HOSTFILE`

poe ./gene_hydra


### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np $NTASKS --ppn 16 --mps 4

