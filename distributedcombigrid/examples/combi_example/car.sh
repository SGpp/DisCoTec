#!/bin/bash
cd ../../../
./compile.sh
cd distributedcombigrid/examples/combi_example
cp Make_copy Makefile
make clean
make
./run.sh