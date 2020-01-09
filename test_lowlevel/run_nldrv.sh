#!/bin/bash

if [ $# -ge 1 ]; then
    nmpi=$1
else
    nmpi=8
fi

declare -i max_rank=$nmpi-1
export I_MPI_PIN_PROCESSOR_LIST=0-${max_rank}
export I_MPI_DEBUG=5
export PERFLIB_OUTPUT_FORMAT=xml
executable=./nldrv
echo "Running $executable with $nmpi MPI ranks."
mpiexec -n ${nmpi} $executable
