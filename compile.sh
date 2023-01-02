cmake -S . -B build
cmake --build build

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
