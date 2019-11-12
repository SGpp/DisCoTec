#!/bin/bash
cd ../../../
./compile.sh
cd distributedcombigrid/examples/combi_example
make clean
make
./run.sh