#!/bin/bash
cmake -S . -B build
cmake --build build

#old scons
#Release
#scons -j 8 ISGENE=0 VERBOSE=1 COMPILE_BOOST_TESTS=1 RUN_BOOST_TESTS=1 RUN_CPPLINT=0 BUILD_STATICLIB=0 CC=/opt/mpich/bin/mpicc FC=mpifort.mpich CXX=/opt/mpich/bin/mpic++ OPT=1 TIMING=0  USE_VTK=0

#Helium
#scons -j 8 ISGENE=0 VERBOSE=1 COMPILE_BOOST_TESTS=1 RUN_BOOST_TESTS=0 RUN_CPPLINT=0 BUILD_STATICLIB=0 CC=mpicc.mpich FC=mpifort.mpich CXX=mpicxx.mpich OPT=1 TIMING=0 UNIFORMDECOMPOSITION=1 ENABLEFT=0 USE_VTK=0

#HLRS
#scons -j 8 ISGENE=0 VERBOSE=1 COMPILE_BOOST_TESTS=1 RUN_BOOST_TESTS=0 RUN_CPPLINT=0 BUILD_STATICLIB=0 CC=cc FC= CXX=CC OPT=1 TIMING=0 UNIFORMDECOMPOSITION=1 ENABLEFT=0 USE_VTK=0

#LRZ
#scons -j 8 ISGENE=0 VERBOSE=1 COMPILE_BOOST_TESTS=1 RUN_BOOST_TESTS=0 RUN_CPPLINT=0 BUILD_STATICLIB=0 CC=mpicc FC=mpifort.mpich CXX=mpiCC OPT=1 TIMING=0 UNIFORMDECOMPOSITION=1 ENABLEFT=0 USE_VTK=0
#mpiexec.mpich -np 9 .tests/test_distributedcombigrid_boost --run_test=ftolerance
mpiexec.mpich -np 9 .tests/test_distributedcombigrid_boost --run_test=fullgrid
mpiexec.mpich -np 9 .tests/test_distributedcombigrid_boost --run_test=hierarchization
mpiexec.mpich -np 9 .tests/test_distributedcombigrid_boost --run_test=loadmodel
mpiexec.mpich -np 9 .tests/test_distributedcombigrid_boost --run_test=reduce
mpiexec.mpich -np 9 .tests/test_distributedcombigrid_boost --run_test=rescheduling
mpiexec.mpich -np 9 .tests/test_distributedcombigrid_boost --run_test=stats
mpiexec.mpich -np 9 .tests/test_distributedcombigrid_boost --run_test=task
mpiexec.mpich -np 9 .tests/test_distributedcombigrid_boost --run_test=distributedfullgrid
mpiexec.mpich -np 9 .tests/test_distributedcombigrid_boost --run_test=integration
mpiexec.mpich -np 9 .tests/test_distributedcombigrid_boost --run_test=distributedsparsegrid
mpiexec.mpich -np 9 .tests/test_distributedcombigrid_boost --run_test=mpisystem