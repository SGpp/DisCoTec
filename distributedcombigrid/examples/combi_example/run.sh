#!/bin/bash
export LD_LIBRARY_PATH=/path/to/SGpp/lib/sgpp:$LD_LIBRARY_PATH
mpirun.mpich -n 5 ./combi_example
