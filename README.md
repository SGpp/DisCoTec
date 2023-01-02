What is this project about?
---------------------------
[![Build Status](https://simsgs.informatik.uni-stuttgart.de/jenkins/buildStatus/icon?job=DisCoTec-master)](https://simsgs.informatik.uni-stuttgart.de/jenkins/view/DisCoTec/job/DisCoTec-master/)

This project contains __DisCoTec__, a code for running the distributed sparse grid combination technique with MPI parallelization. While it originates from the excellent [SGpp project](https://github.com/SGpp/SGpp), all the parallelization makes it a very different code, such that it is its own project now.

Guidelines
---------
*  Please develop new features in a new branch and then create a merge request 
to notify other users about the changes. So everyone can check whether there are 
side-effects to other branches.
* Of course, you are still very welcome to directly fix obvious bugs in the 
master branch.
* Before merging new features to the master branch, please make sure that they
are sufficiently commented and extensively tested.

Requirements
--------------
cmake >= (3.24.2),
libmpich-dev (>= 3.2-6), or other MPI library

Additional (optional) dependencies (can be optained with cmake):
- OpenMP
- libboost-all-dev (>= 1.60)
- HDF5
- HighFive 

Installation instructions: 
--------------------------
#### Complete build:
```bash
git clone https://github.com/SGpp/DisCoTec.git DisCoTec
cmake -S DisCoTec -B DisCoTec/build
cmake --build DisCoTec/build -j
```

#### Advanced build:

All examples and the tools can be build on their own. To build a specific example, run the specific cmake command in the example folder. For example, run the following commands
```bash
git clone https://github.com/SGpp/DisCoTec.git DisCoTec
cd DisCoTec/examples/combi_example
cmake -S . -B build
cmake --build build -j
```
to build the combi_example.


For building only the DisCoTec library, run cmake with the `src` folder as source folder.

#### Optional CMake Options
- `DISCOTEC_TEST=**ON**|OFF` - Build tests if you build the complete project.
- `DISCOTEC_BUILD_MISSING_DEPS=**ON**|OFF`- First order dependencies that are not found are built automatically (glpk is always built).
- `DISCOTEC_TIMING=**ON**|OFF` - Enables internal timing
- `DISCOTEC_USE_HDF5=**ON**|OFF`
- `DISCOTEC_USE_HIGHFIVE=**ON**|OFF` - Enables HDF5 support via HighFive. If `DISCOTEC_USE_HIGHFIVE=ON`, `DISCOTEC_USE_HDF5` has also to be `ON`.
- `DISCOTEC_UNIFORMDECOMPOSITION=**ON **|OFF` - Enables the uniform decomposition of the grid.
- `DISCOTEC_GENE=ON|**OFF**` - Currently GEne is not supported with CMake!
- `DISCOTEC_OPENMP=ON|**OFF**` - Enables OpenMP support.
- `DISCOTEC_ENABLEFT=**ON**|OFF` - Enables the use of the FT library.
- `DISCOTEC_DEBUG=ON|**OFF**` - Enables debug mode. Has to be used with a cmake debug build mode.
- `DISCOTEC_USE_LTO=**ON**|OFF` - Enables link time optimization if the compiler supports it.
- `DISCOTEC_OMITREADYSIGNAL=ON|**OFF**` - Omit the ready signal in the MPI communication. This can be used to reduce the communication overhead.
- `DISCOTEC_USENONBLOCKINGMPICOLLECTIVE=ON|**OFF**` - TODO: Add description


To run the compiled tests, go to folder `tests` and run
```bash
mpiexec -np 9 ./test_distributedcombigrid_boost
```
where you can use all the parameters of the boost test suite.
Or you can run the tests with `ctest` in the build folder.


GENE  submodules as dependencies for GENE examples
----------------
_Warning: The CMake Integration is currently not adappted to use GENE and SeLaLib!_

There are gene versions as submodules: a linear one in the gene_mgr folder, and 
a nonlinear one in gene-non-linear. To get them, you need access to their 
respective repos. Then go into the folder and

``` bash
git submodule init
git submodule update
```
or use the `--recursive` flag when cloning the repo.
and then you just hope you can `make` it by adapting the local library flags.