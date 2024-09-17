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

![schematic of a combination scheme in 2D](gfx/combischeme-2d.svg)

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
technique's parallelism.
This is achieved through *process groups* (pgs):
`MPI_COMM_WORLD` is subdivided into equal-sized process groups
(and optionally, a manager rank).

![schematic of MPI ranks in DisCoTec](gfx/discotec-ranks.svg)

The image describes the two ways of scaling up:
One can either increase the size or the number of process groups.
Figure originally published in (Pollinger [2024](https://elib.uni-stuttgart.de/handle/11682/14229)).

All component grids are distributed to process groups (either statically, or
dynamically through the manager rank).
During the solver time step and most of the combination, MPI communication only
happens within the process groups.
Conversely, for the sparse grid reduction using the combination coefficients,
MPI communication only happens between a rank and its colleagues in the other
process groups, e.g., rank 0 in group 0 will only talk to rank 0 in all other groups.
Thus, major bottlenecks arising from global communication can be avoided altogether.

Combining the two ways of scaling up, DisCoTec's scalability was demonstrated on
several machines, with the experiments comprising up to 524288 cores:

![timings for advection solver step on HAWK at various
parallelizations](gfx/times-solver-on-hawk.svg)![timings for combination step on
HAWK at various parallelizations](gfx/times-combination-on-hawk.svg)

We see the timings (in seconds) for the advection solver step and the
combination step, respectively.
This weak scaling experiment used four OpenMP threads per rank, and starts with
one pg of four processes in the upper left corner.
The largest parallelization is 64 pgs of 2048 processes each.
Figure originally published in (Pollinger [2024](https://elib.uni-stuttgart.de/handle/11682/14229)).

This image also describes the challenges in large-scale experiments with DisCoTec:
If the process groups become too large, the MPI communication of the multiscale
transform starts to dominate the combination time.
If there are too many pgs, the combination reduction will dominate the
combination time.
However, the times required for the solver stay relatively constant;
they are determined by the solver's own scaling and the load balancing quality.

There are only few codes that allow weak scaling up to this problem size:
a size that uses most of the available main memory of the entire system.

## Contributing

We welcome contributions! To find a good place to start coding, have a look at
the currently open issues.

- Please describe issues and intended changes in the [issue tracker](https://github.com/SGpp/DisCoTec/issues).
- Please develop new features in a new branch (typically on your own fork) and
  then create a [pull request](https://github.com/SGpp/DisCoTec/pulls).
- New features will only be merged to the main branch if they are sufficiently
  tested: please add unit tests in [/tests](/tests).

## Installing

DisCoTec can be installed via spack, which handles all dependencies.
We recommend the `spack dev-build` workflow:

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

[Here are the Docs](https://discotec.readthedocs.io/en/latest/getting_started.html#installation-with-spack) for CMake options and further Spack customization hints.


