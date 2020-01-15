#!/bin/bash
cd ../../../
./compile.sh
cd distributedcombigrid/examples/dealii_example

make clean
make
. ./run2.sh