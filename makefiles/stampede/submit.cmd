#!/bin/bash
#SBATCH -J GENE        # Job Name
#SBATCH -o GENE.out%j    # Output and error file name (%j expands to jobID)
#SBATCH -n 32           # Total number of mpi tasks requested
#SBATCH -p development  # Queue (partition) name -- normal, development, etc.
#SBATCH -t 01:30:00     # Run time (hh:mm:ss) - 1.5 hours
###./scanscript --np 32 --ppn 16 --syscall 'ibrun ./gene_stampede'
ibrun ./gene_stampede           # Run the MPI executable named a.out






