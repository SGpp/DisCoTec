#!/bin/bash

NGROUP=$(grep ngroup ctparam | awk -F"=" '{print $2}')
NPROCS=$(grep nprocs ctparam | awk -F"=" '{print $2}')

mpiexec -n $(($NGROUP*$NPROCS)) ./combi_example_worker_only
