#!/bin/bash
# Job Name and Files (also --job-name)
#SBATCH -J third_level_manager_and_broker
#Output and error (also --output, --error):
#SBATCH -o ./%x.%j.out
#SBATCH -e ./%x.%j.err
#Initial working directory (also --chdir):
#SBATCH -D ./
#Notification and type
#SBATCH --mail-type=END
#SBATCH --mail-user=pollinta@ipvs.uni-stuttgart.de
# Wall clock limit:
#SBATCH --time=00:20:00
#SBATCH --no-requeue

#SBATCH --ntasks=2

export LD_LIBRARY_PATH=$(pwd)/../../../lib/sgpp:$(pwd)/../../../glpk/lib:$LD_LIBRARY_PATH

. ./setenv.sh

mpiexec.openmpi -n 1 ./manager_only ctparam_tl_system1 : -n 1 ../../third_level_manager/thirdLevelManager ../../third_level_manager/example.ini
