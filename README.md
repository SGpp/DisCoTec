# What is DisCoTec?
---------------------------
[![Build Status](https://simsgs.informatik.uni-stuttgart.de/jenkins/buildStatus/icon?job=DisCoTec-main-spack)](https://simsgs.informatik.uni-stuttgart.de/jenkins/view/DisCoTec/job/DisCoTec-main-spack/)
[![DOI](https://zenodo.org/badge/226341053.svg)](https://zenodo.org/badge/latestdoi/226341053)

This project contains __DisCoTec__, a code for running the distributed sparse grid combination technique with MPI parallelization. While it originates from the excellent [SGpp project](https://github.com/SGpp/SGpp), all the parallelization makes it a very different code, such that it is its own project now.

It is designed as a framework that can run multiple instances of a (black-box) grid-based solver implementation.
The most basic example we use is a [mass-conserving FDM/FVM constant advection upwinding solver](/examples/distributed_advection/).
An example of a separate, coupled solver is [SeLaLib](/examples/selalib_distributed/).

## Contributing
---------
* We welcome contributions! To find a good place to start, have a look at the currently open issues
* Please develop new features in a new branch (typically on your own fork) and then create a PR
* New features will only be merged to the main branch if they are sufficiently tested: please add unit tests in [/tests] .


## Installation instructions
--------------------------
#### Dependencies
--------------
cmake >= (3.24.2),
libmpich-dev (>= 3.2-6), or other MPI library

Additional (optional) dependencies:
- OpenMP
- libboost-all-dev (>= 1.60)
- HDF5
- HighFive 

For the dependencies, the DisCoTec spack package (currently under review [here](https://github.com/spack/spack/pull/35239)) could be of interest.

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
- `DISCOTEC_ENABLEFT=ON|**OFF**` - Enables the use of the FT library.
- `DISCOTEC_USE_LTO=**ON**|OFF` - Enables link time optimization if the compiler supports it.
- `DISCOTEC_OMITREADYSIGNAL=ON|**OFF**` - Omit the ready signal in the MPI communication. This can be used to reduce the communication overhead.
- `DISCOTEC_USENONBLOCKINGMPICOLLECTIVE=ON|**OFF**` - TODO: Add description
- `DISCOTEC_WITH_SELALIB=ON|**OFF**` - Looks for SeLaLib dependencies and compiles [the matching example](/examples/selalib_distributed/)


To run the compiled tests, go to folder `tests` and run
```bash
mpiexec -np 9 ./test_distributedcombigrid_boost
```
where you can use all the parameters of the boost test suite.
If timings matter, consider the pinning described in the respective section.
Or you can run the tests with `ctest` in the build folder.

## Executing DisCoTec Binaries / Rank and Thread Pinning
The number and size of process groups (in MPI ranks) can be read from the `ctparam` files:
```
NGROUP=$(grep ngroup ctparam | awk -F"=" '{print $2}')
NPROCS=$(grep nprocs ctparam | awk -F"=" '{print $2}')
```

Correct pinning is of utmost importance, especially if DisCoTec was compiled with OpenMP support.
The desired pinning is the simplest one you can imagine: Each node and, within a node, each socket is filled consecutively with MPI rank numbers.
If OpenMP is used, MPI pinning is strided by the number of threads. The threads should be placed close to their MPI rank.
This means for all compilers, if the desired number of OpenMP threads per rank is `N`, one needs
```
export OMP_NUM_THREADS=$N;export OMP_PLACES=cores;export OMP_PROC_BIND=close
```
If OpenMP is used for compilation but should not be used for execution, `N` should be set to 1 to avoid unintended effects.

#### Intel-MPI
Intel-MPI requires some [environment variables](https://software.intel.com/content/www/us/en/develop/documentation/mpi-developer-reference-linux/top/environment-variable-reference/process-pinning/environment-variables-for-process-pinning.html), in particular [for OpenMP](https://www.intel.com/content/www/us/en/docs/mpi-library/developer-guide-linux/2021-6/running-an-mpi-openmp-program.html):
```
mpiexec -n $(($NGROUP*$NPROCS)) -genv I_MPI_PIN_PROCESSOR_LIST="allcores" -genv I_MPI_PIN_DOMAIN=omp $DISCOTEC_EXECUTABLE
```
In SLURM, it is important to set --ntasks-per-node to match the number of desired tasks ($CORES_PER_NODE/$OMP_NUM_THREADS). 
Validate with verbose output: `export I_MPI_DEBUG=4` .

#### OpenMPI
OpenMPI uses command line arguments, pinning may clash with SLURM settings depending on the exact call.

```
mpiexec.openmpi --rank-by core --map-by node:PE=${OMP_NUM_THREADS} -n $(($NGROUP*$NPROCS)) $DISCOTEC_EXECUTABLE
```
Validate with verbose output: `--report-bindings` .

#### MPT
With mpt, one uses the `omplace` wrapper tool to set the correct pinning.
```
mpirun -n $(($NGROUP*$NPROCS)) omplace -v -nt ${OMP_NUM_THREADS} $DISCOTEC_EXECUTABLE
```
Validate with very verbose output: `-vv` .

### GENE  submodules as dependencies for GENE examples
----------------
_Warning: The CMake Integration is currently not adapted to use GENE!_

There are gene versions as submodules: a linear one in the gene_mgr folder, and 
a nonlinear one in gene-non-linear. To get them, you need access to their 
respective repos. Then go into the folder and

``` bash
git submodule init
git submodule update
```
or use the `--recursive` flag when cloning the repo.
and then you just hope you can `make` it by adapting the local library flags.


### Further Readmes
- [widely-distribued simulations](/third_level_manager/README.md)
- [the subspace writer tool](/tools/subspace_writer/README.md)
- [selalib simulations with DisCoTec](/examples/selalib_distributed/README.md)
