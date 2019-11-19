#!/bin/bash
cd ../../../
./compile.sh
cd distributedcombigrid/examples/dealii_example
cp Make_copy Makefile
make clean
make
./run.sh