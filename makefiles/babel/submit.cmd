# @ job_name = GENE
# @ job_type = BLUEGENE
# Standard output
# @ output = $(job_name).$(jobid)
# Standard error output
# @ error = $(output)
# Time limit
# @ wall_clock_limit = 20:00:00
# Number of compute nodes (=n_mpiprocs/4 if VN mode is used)
# @ bg_size = 64
# @ queue

# Copy executable and input file to TMPDIR
# Warning: if you need to transfer important volumes
# of data, please use a multi-step job
#cp my_code $TMPDIR
#cp data.in $TMPDIR
#cd $TMPDIR

mpirun -mode VN -np 256 -mapfile TXYZ -exe ./gene_babel

# Copy output file to submission directory
# Warning: if you need to transfer important volumes
# of data, please use a multi-step job
# $LOADL_STEP_INITDIR is the submission directory
#cp data.out $LOADL_STEP_INITDIR/

### to submit a parameter scan, comment the mpirun... line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np 256 --ppn 4 --mps 4
