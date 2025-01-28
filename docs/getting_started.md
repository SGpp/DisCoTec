# Getting Started with DisCoTec

## Installing DisCoTec

There are two established ways of installing: Spack and CMake.
After installing all dependencies, Spack actually delegates to
CMake to build DisCoTec.

### ... with Spack

DisCoTec can be installed via spack, which handles all dependencies.

If you want to use DisCoTec as a library only, follow the usual
[Spack installation](https://spack.readthedocs.io/en/latest/features.html)
instructions until the command `spack install discotec@main`
(or an extended version of it) succeeds.

If you want to develop DisCoTec code or examples, the `spack dev-build` workflow
is recommended as follows.

Clone both `spack` and `DisCoTec` to find or build the dependencies and then
compile DisCoTec:

```bash
git clone https://github.com/spack/spack.git
./spack/bin/spack external find  # find already-installed packages
./spack/bin/spack compiler find  # find compilers present on system
./spack/bin/spack info discotec@main  # shows DisCoTec's variants

 # shows DisCoTec's dependency tree and which parts are already found
./spack/bin/spack spec discotec@main
# optional: configuring external packages according to Spack documentation

git clone https://github.com/SGpp/DisCoTec.git
cd DisCoTec
../spack/bin/spack dev-build -b install discotec@main
```

This will first build all dependencies, and then build DisCoTec inside the
cloned folder.
The executables are placed in the respective `example` and `test` folders.

If you encounter problems with the Spack installation, check out
the [Spack command documentation](https://spack.readthedocs.io/en/latest/getting_started.html),
and see if other users had the same
[spack issues](https://github.com/spack/spack/issues?q=is%3Aissue) or
[DisCoTec issues](https://github.com/SGpp/DisCoTec/issues) before.

### ... with CMake

#### Dependencies

cmake >= (3.24.2),
libmpich-dev (>= 3.2-6), or other MPI library
libboost-all-dev (>= 1.60) (or the subset test, serialization, filesystem, system,
program_options, date_time Boost libraries)

Additional (optional) dependencies:

- OpenMP
- HDF5
- HighFive
- lz4
- glpk (used version bundled in the DisCoTec repo)
- vtk

You can also install the dependencies with Spack, and `spack load`
them, before executing the next steps.

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

to build the `combi_example`.

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
  compiles [the matching example](https://github.com/SGpp/DisCoTec/blob/main/examples/selalib_distributed/)

## Run the tests

Beware: Most of the tests in DisCoTec need significant resources
(eight to nine MPI processes,
file I/O, some performance tests when built in Release mode).
Make sure you run them on appropriate hardware.

To run the compiled tests, go to folder `tests` and make sure the correct MPI runtime
is loaded to your path, e.g. with

```bash
. ../spack/share/spack/setup-env.sh
spack load openmpi # <- use the one you compiled with
mpiexec --version # <- expected output?
```

Then, run

```bash
mpiexec -np 9 ./test_distributedcombigrid_boost
```

where you can use all the [parameters of the boost test suite](https://beta.boost.org/doc/libs/1_60_0/libs/test/doc/html/boost_test/utf_reference/rt_param_reference.html).
If timings matter, consider [pinning](#pinning-with-various-mpi-implementations).

## Run an Example

DisCoTec executables are typically configured through `ctparam` files, which are
parsed on startup.
The `ctparam` file will contain the combination technique parameters (dimension,
minimum and maximum level) as well as parallelization parameters (number and
size of process groups, parallelism within process groups) in `.ini` file format.

Any DisCoTec executable must be run through MPI (either `mpirun` or `mpiexec`),
and if no argument to the binary is specified, it will use the file called
`ctparam` in the current working directory.
Make sure that it exists and describes the parameters you want to run.

As with the tests, make sure the correct MPI is loaded to your path.

The exact format and naming in `ctparam` is not (yet) standardized, to allow
adaptation for different PDE solver applications.
Please refer to existing parameter files and example implementations.

The different examples intend to illustrate various aspects of DisCoTec:

- [`combi_example`](https://github.com/SGpp/DisCoTec/blob/main/examples/combi_example/):
  manager-worker setup with a generic `TaskExample`
  that interpolates a cosine function on the component grids
- [`combi_example_faults`](https://github.com/SGpp/DisCoTec/blob/main/examples/combi_example_faults/):
  like `combi_example`, but with fault tolerance functionality
  (requires that DisCoTec is built with the `DISCOTEC_ENABLEFT` flag)
- [`distributed_advection`](https://github.com/SGpp/DisCoTec/blob/main/examples/distributed_advection/):
  manager-worker setup using a first-order upwinding advection solver,
  illustrates how to evaluate the Monte-Carlo norms
  (and error norms, as an analytical solution is given for this PDE problem)
- [`combi_workers_only`](https://github.com/SGpp/DisCoTec/blob/main/examples/combi_workers_only/):
  like `distributed_advection`, but with worker-only setup;
  illustrates the static task assignment required for worker-only scenarios
- [`distributed_third_level`](https://github.com/SGpp/DisCoTec/blob/main/examples/distributed_third_level/):
  like `distributed_advection` and `combi_workers_only`;
  adds logic for widely-distributed simulations
- GENE examples, starting with `gene_`: show the interaction with
  the GENE gyrokinetic solver (using a complicated build/link setup
  that we don't recommend to use any more, see next example instead)
- [`selalib_distributed`](https://github.com/SGpp/DisCoTec/blob/main/examples/selalib_distributed/):
  shows the interaction with a Semi-Lagrangian Vlasov solver in SeLaLib:
  a static library is built as part of SeLaLib, and linked through DisCoTec's
  `Task` interface in `SelalibTask`.

(pinning-with-various-mpi-implementations)=
### Pinning With Various MPI Implementations

Correct pinning is of utmost importance to performance, especially if DisCoTec
was compiled with OpenMP support.
The desired pinning is the simplest one you can imagine: Each node and, within a
node, each socket is filled consecutively with MPI rank numbers.
If OpenMP is used, MPI pinning is strided by the number of threads.
The threads should be placed close to their MPI rank.
This means for all compilers, if the desired number of OpenMP threads per rank
is `N`, one needs

```bash
export OMP_NUM_THREADS=$N;export OMP_PLACES=cores;export OMP_PROC_BIND=close
```

If OpenMP is used for compilation but should not be used for execution, `N`
should be set to 1 to avoid unintended effects.

The number and size of process groups (in MPI ranks) can be read from the
`ctparam` files:

```bash
NGROUP=$(grep ngroup ctparam | awk -F"=" '{print $2}')
NPROCS=$(grep nprocs ctparam | awk -F"=" '{print $2}')
```

#### Intel-MPI

Intel-MPI requires some [environment variables](https://software.intel.com/content/www/us/en/develop/documentation/mpi-developer-reference-linux/top/environment-variable-reference/process-pinning/environment-variables-for-process-pinning.html),
in particular [for OpenMP](https://www.intel.com/content/www/us/en/docs/mpi-library/developer-guide-linux/2021-6/running-an-mpi-openmp-program.html):

```bash
mpiexec -n $(($NGROUP*$NPROCS)) -genv I_MPI_PIN_PROCESSOR_LIST="allcores" -genv I_MPI_PIN_DOMAIN=omp $DISCOTEC_EXECUTABLE
```

In SLURM, it is important to set `--ntasks-per-node` to match the number of
desired tasks (`$CORES_PER_NODE`/`$OMP_NUM_THREADS`).
Validate with verbose output: `export I_MPI_DEBUG=4`

#### OpenMPI

OpenMPI uses command line arguments, pinning may clash with SLURM settings
depending on the exact call.

```bash
mpiexec.openmpi --rank-by core --map-by node:PE=${OMP_NUM_THREADS} -n $(($NGROUP*$NPROCS)) $DISCOTEC_EXECUTABLE
```

Validate with verbose output: `--report-bindings`

#### MPT

With mpt, one uses the `omplace` wrapper tool to set the correct pinning.

```bash
mpirun -n $(($NGROUP*$NPROCS)) omplace -v -nt ${OMP_NUM_THREADS} $DISCOTEC_EXECUTABLE
```

Validate with very verbose output: `-vv` .
