
[![Build Status](https://jenkins-sim.informatik.uni-stuttgart.de/buildStatus/icon?job=DisCoTec%2Fmain)](https://jenkins-sim.informatik.uni-stuttgart.de/job/DisCoTec/job/main/)
[![Zenodo DOI](https://zenodo.org/badge/226341053.svg)](https://zenodo.org/badge/latestdoi/226341053)
[![Latest spack version](https://img.shields.io/spack/v/discotec)](https://spack.readthedocs.io/en/latest/package_list.html#discotec)

# What is DisCoTec?
---------------------------

This project contains __DisCoTec__, a code for running the *dis*tributed sparse grid *co*mbination *tec*hnique with MPI parallelization. 
While it originates from the excellent [SGpp project](https://github.com/SGpp/SGpp), all the parallelization makes it a very different code, such that it has become its own project.

DisCoTec is designed as a framework that can run multiple instances of a (black-box) grid-based solver implementation.
The most basic example we use is a [mass-conserving FDM/FVM constant advection upwinding solver](/examples/distributed_advection/).
An example of a separate, coupled solver is [SeLaLib](/examples/selalib_distributed/).


## Contributing
---------
* We welcome contributions! To find a good place to start coding, have a look at the currently open issues.
* Please describe issues in the [issue tracker](https://github.com/SGpp/DisCoTec/issues).
* Please develop new features in a new branch (typically on your own fork) and then create a [pull request](https://github.com/SGpp/DisCoTec/pulls).
* New features will only be merged to the main branch if they are sufficiently tested: please add unit tests in [/tests] .


## Installation instructions: spack
--------------------------

DisCoTec can be installed via Spack, which handles all dependencies.
If you want to develop DisCoTec code or examples, the `spack dev-build` workflow is recommended.

Clone both `spack` and `DisCoTec` to find or build the dependencies and then compile DisCoTec:
```bash
git clone git@github.com:spack/spack.git  # use https if your ssh is not set up on github
./spack/bin/spack external find  # find already-installed packages
./spack/bin/spack compiler find  # find compilers present on system
./spack/bin/spack info discotec@main  # shows DisCoTec's variants 
./spack/bin/spack spec discotec@main  # shows DisCoTec's dependency tree and which parts are already found

git clone git@github.com:SGpp/DisCoTec.git
cd DisCoTec
../spack/bin/spack dev-build -b install discotec@main
```
This will build DisCoTec inside the cloned folder, and the executables are placed in the respective `example` and `test` folders.


## Installation instructions: CMake
--------------------------
#### Dependencies
--------------
cmake >= (3.24.2),
libmpich-dev (>= 3.2-6), or other MPI library
libboost-all-dev (>= 1.60)

Additional (optional) dependencies:
- OpenMP
- HDF5
- HighFive 

#### Complete build:
```bash
git clone https://github.com/SGpp/DisCoTec.git DisCoTec
cmake -S DisCoTec -B DisCoTec/build
cmake --build DisCoTec/build -j
```

#### Advanced build:

All examples and the tools can be built on their own. 
To build a specific example, run the specific cmake command in the example folder. 
For example, run the following commands
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
respective repos at MPCDF. Then go into the folder and

``` bash
git submodule init
git submodule update
```
or use the `--recursive` flag when cloning the repo.


### Further Readmes
- [widely-distribued simulations](/third_level_manager/README.md)
- [the subspace writer tool](/tools/subspace_writer/README.md)
- [selalib simulations with DisCoTec](/examples/selalib_distributed/README.md)
