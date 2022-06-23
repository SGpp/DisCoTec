scons -j 8 ISGENE=0 VERBOSE=1 COMPILE_BOOST_TESTS=1 RUN_BOOST_TESTS=0 RUN_CPPLINT=0 BUILD_STATICLIB=0 CC=mpicc FC=mpifort CXX=mpicxx OPT=1 TIMING=0 UNIFORMDECOMPOSITION=1 ENABLEFT=0 DOC=0 #DEBUG_OUTPUT=1

export LD_LIBRARY_PATH=$(pwd)/lib/sgpp:$(pwd)/glpk/lib:$LD_LIBRARY_PATH
mpiexec -np 9 ./distributedcombigrid/tests/test_distributedcombigrid_boost --run_test=ftolerance
mpiexec -np 9 ./distributedcombigrid/tests/test_distributedcombigrid_boost --run_test=fullgrid
mpiexec -np 9 ./distributedcombigrid/tests/test_distributedcombigrid_boost --run_test=hierarchization
mpiexec -np 9 ./distributedcombigrid/tests/test_distributedcombigrid_boost --run_test=loadmodel
mpiexec -np 9 ./distributedcombigrid/tests/test_distributedcombigrid_boost --run_test=reduce
mpiexec -np 9 ./distributedcombigrid/tests/test_distributedcombigrid_boost --run_test=rescheduling
mpiexec -np 9 ./distributedcombigrid/tests/test_distributedcombigrid_boost --run_test=stats
mpiexec -np 9 ./distributedcombigrid/tests/test_distributedcombigrid_boost --run_test=task
mpiexec -np 9 ./distributedcombigrid/tests/test_distributedcombigrid_boost --run_test=distributedfullgrid
mpiexec -np 9 ./distributedcombigrid/tests/test_distributedcombigrid_boost --run_test=integration
mpiexec -np 9 ./distributedcombigrid/tests/test_distributedcombigrid_boost --run_test=distributedsparsegrid
mpiexec -np 9 ./distributedcombigrid/tests/test_distributedcombigrid_boost --run_test=adaptive
