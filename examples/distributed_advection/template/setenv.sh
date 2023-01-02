#!/bin/bash
export LD_LIBRARY_PATH=~/epyc/DisCoTec-advection-mass/lib/sgpp:~/epyc/DisCoTec-advection-mass/glpk/lib:$LD_LIBRARY_PATH

. /home/pollinta/epyc/spack/share/spack/setup-env.sh
spack load boost~debug@1.74.0
# spack load mpich@3.3.1
