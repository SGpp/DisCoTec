# DisCoTec: Distributed Combination Technique Framework

[![Build Status](https://jenkins-sim.informatik.uni-stuttgart.de/buildStatus/icon?job=DisCoTec%2Fmain)](https://jenkins-sim.informatik.uni-stuttgart.de/job/DisCoTec/job/main/)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![Zenodo DOI](https://zenodo.org/badge/226341053.svg)](https://zenodo.org/badge/latestdoi/226341053)
[![Latest spack version](https://img.shields.io/spack/v/discotec)](https://spack.readthedocs.io/en/latest/package_list.html#discotec)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/cac5bc0841784657b2bb75ea46e7cf01)](https://app.codacy.com/gh/SGpp/DisCoTec/dashboard)

## What is DisCoTec?

This project contains **DisCoTec**, a code for running the *dis*tributed sparse
grid *co*mbination *tec*hnique with MPI parallelization.
While it originates from the excellent
[SGpp project](https://github.com/SGpp/SGpp), all the parallelization makes it a
very different code, such that it has become its own project.

DisCoTec is designed as a framework that can run multiple instances of a
(black-box) grid-based solver implementation.
The most basic example we use is a [mass-conserving FDM/FVM constant advection
upwinding solver](/examples/distributed_advection/).
An example of a separate, coupled solver is [SeLaLib](/examples/selalib_distributed/).

### Sparse Grid Combination Technique with Time Stepping

The sparse grid combination technique (Griebel et al.
[1992](https://ins.uni-bonn.de/media/public/publication-media/griesiam.ps.gz),
Garcke [2013](https://link.springer.com/chapter/10.1007/978-3-642-31703-3_3),
Harding [2016](https://link.springer.com/chapter/10.1007/978-3-319-28262-6_4))
can be used to alleviate the curse of dimensionality encountered in
high-dimensional simulations.
Instead of using your solver on a single structured full grid (where every
dimension is finely resolved), you would use it on many different structured
full grids (each of them differently resolved).
We call these coarsely-resolved grids component grids.
Taken together, all component grids form a sparse grid approximation, which can
be explicitly obtained by a linear superposition of the individual grid
functions, with the so-called combination coefficients.

![schematic of a combination scheme in 2D](gfx/combischeme-2d.pdf)
In this two-dimensional combination scheme, all combination coefficients are 1
and -1, respectively.
Figure originally published in (Pollinger [2024](https://elib.uni-stuttgart.de/handle/11682/14229)).

Between time steps, the grids exchange data through a multi-scale approach,
which is summarized as the "combination" step in DisCoTec.
Assuming a certain smoothness in the solution, this allows for a good
approximation of the finely-resolved function, while achieving drastic
reductions in compute and memory requirements.

### Parallelism in DisCoTec

The DisCoTec framework can work with existing MPI parallelized solver codes
operating on structured grids.
In addition to the parallelism provided by the solver, it adds the combination
technique's parallelism:

![schematic of MPI ranks in DisCoTec](gfx/discotec-ranks.pdf)

## Contributing

We welcome contributions! To find a good place to start coding, have a look at
the currently open issues.

- Please describe issues and intended changes in the [issue tracker](https://github.com/SGpp/DisCoTec/issues).
- Please develop new features in a new branch (typically on your own fork) and
  then create a [pull request](https://github.com/SGpp/DisCoTec/pulls).
- New features will only be merged to the main branch if they are sufficiently
  tested: please add unit tests in [/tests](/tests).

## Installing

### Installation instructions: spack

DisCoTec can be installed via spack, which handles all dependencies.
If you want to develop DisCoTec code or examples, the `spack dev-build` workflow
is recommended as follows.

Clone both `spack` and `DisCoTec` to find or build the dependencies and then
compile DisCoTec:

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

This will first build all dependencies, and then build DisCoTec inside the
cloned folder.
The executables are placed in the respective `example` and `test` folders.

### Installation instructions: CMake

#### Dependencies

cmake >= (3.24.2),
libmpich-dev (>= 3.2-6), or other MPI library
libboost-all-dev (>= 1.60)

Additional (optional) dependencies:

- OpenMP
- HDF5
- HighFive

#### Complete build

```bash
git clone https://github.com/SGpp/DisCoTec.git DisCoTec
cmake -S DisCoTec -B DisCoTec/build
cmake --build DisCoTec/build -j
```

#### Advanced build

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

For building only the DisCoTec library, run cmake with the `src` folder as
source folder.

#### Optional CMake Options

- `DISCOTEC_TEST=**ON**|OFF` - Build tests if you build the complete project.
- `DISCOTEC_BUILD_MISSING_DEPS=**ON**|OFF`- First order dependencies that are
  not found are built automatically (glpk is always built).
- `DISCOTEC_TIMING=**ON**|OFF` - Enables internal timing
- `DISCOTEC_USE_HDF5=**ON**|OFF`
- `DISCOTEC_USE_HIGHFIVE=**ON**|OFF` - Enables HDF5 support via HighFive. If
  `DISCOTEC_USE_HIGHFIVE=ON`, `DISCOTEC_USE_HDF5` has also to be `ON`.
- `DISCOTEC_UNIFORMDECOMPOSITION=**ON **|OFF` - Enables the uniform
  decomposition of the grid.
- `DISCOTEC_GENE=ON|**OFF**` - Currently GEne is not supported with CMake!
- `DISCOTEC_OPENMP=ON|**OFF**` - Enables OpenMP support.
- `DISCOTEC_ENABLEFT=ON|**OFF**` - Enables the use of the FT library.
- `DISCOTEC_USE_LTO=**ON**|OFF` - Enables link time optimization if the compiler
  supports it.
- `DISCOTEC_OMITREADYSIGNAL=ON|**OFF**` - Omit the ready signal in the MPI
  communication. This can be used to reduce the communication overhead.
- `DISCOTEC_USENONBLOCKINGMPICOLLECTIVE=ON|**OFF**` - Flag currently unused
- `DISCOTEC_WITH_SELALIB=ON|**OFF**` - Looks for SeLaLib dependencies and
  compiles [the matching example](/examples/selalib_distributed/)

To run the compiled tests, go to folder `tests` and run

```bash
mpiexec -np 9 ./test_distributedcombigrid_boost
```

where you can use all the parameters of the boost test suite.
If timings matter, consider the pinning described in the respective section.
Or you can run the tests with `ctest` in the build folder.

## Executing DisCoTec Binaries

TODO what's in ctparam

## Rank and Thread Pinning

The number and size of process groups (in MPI ranks) can be read from the
`ctparam` files:

```bash
NGROUP=$(grep ngroup ctparam | awk -F"=" '{print $2}')
NPROCS=$(grep nprocs ctparam | awk -F"=" '{print $2}')
```

Correct pinning is of utmost importance to performance, especially if DisCoTec
was compiled with OpenMP support.
The desired pinning is the simplest one you can imagine: Each node and, within a
node, each socket is filled consecutively with MPI rank numbers.
If OpenMP is used, MPI pinning is strided by the number of threads. The threads
should be placed close to their MPI rank.
This means for all compilers, if the desired number of OpenMP threads per rank
is `N`, one needs

```bash
export OMP_NUM_THREADS=$N;export OMP_PLACES=cores;export OMP_PROC_BIND=close
```

If OpenMP is used for compilation but should not be used for execution, `N`
should be set to 1 to avoid unintended effects.

### Intel-MPI

Intel-MPI requires some [environment variables](https://software.intel.com/content/www/us/en/develop/documentation/mpi-developer-reference-linux/top/environment-variable-reference/process-pinning/environment-variables-for-process-pinning.html),
in particular [for OpenMP](https://www.intel.com/content/www/us/en/docs/mpi-library/developer-guide-linux/2021-6/running-an-mpi-openmp-program.html):

```bash
mpiexec -n $(($NGROUP*$NPROCS)) -genv I_MPI_PIN_PROCESSOR_LIST="allcores" -genv I_MPI_PIN_DOMAIN=omp $DISCOTEC_EXECUTABLE
```

In SLURM, it is important to set --ntasks-per-node to match the number of
desired tasks (`$CORES_PER_NODE`/`$OMP_NUM_THREADS`).
Validate with verbose output: `export I_MPI_DEBUG=4`

### OpenMPI

OpenMPI uses command line arguments, pinning may clash with SLURM settings
depending on the exact call.

```bash
mpiexec.openmpi --rank-by core --map-by node:PE=${OMP_NUM_THREADS} -n $(($NGROUP*$NPROCS)) $DISCOTEC_EXECUTABLE
```

Validate with verbose output: `--report-bindings` .

### MPT

With mpt, one uses the `omplace` wrapper tool to set the correct pinning.

```bash
mpirun -n $(($NGROUP*$NPROCS)) omplace -v -nt ${OMP_NUM_THREADS} $DISCOTEC_EXECUTABLE
```

Validate with very verbose output: `-vv` .

## GENE submodules as dependencies for GENE examples

*Warning: The CMake Integration is currently not adapted to use GENE!*

There are gene versions as submodules: a linear one in the gene_mgr folder, and
a nonlinear one in gene-non-linear. To get them, you need access to their
respective repos at MPCDF. Then go into the folder and

```bash
git submodule init
git submodule update
```

or use the `--recursive` flag when cloning the repo.

## Further Readmes

- [widely-distribued simulations](/third_level_manager/README.md)
- [the subspace writer tool](/tools/subspace_writer/README.md)
- [selalib simulations with DisCoTec](/examples/selalib_distributed/README.md)
