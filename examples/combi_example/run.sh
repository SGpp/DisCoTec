#!/bin/bash

NGROUP=$(grep ngroup ctparam | awk -F"=" '{print $2}')
NPROCS=$(grep nprocs ctparam | awk -F"=" '{print $2}')

mpirun.mpich -n $(($NGROUP*$NPROCS+1)) -l ./combi_example
