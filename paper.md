---
title: 'DisCoTec: Distributed higher-dimensional HPC simulations with the sparse grid combination technique'
tags:
  - C++
  - MPI
  - structured grid-based simulations
  - sparse grids
  - black-box solvers
  - Vlasov solvers
  - massively parallel
authors:
  - name: Theresa Pollinger
    affiliation: 1 
    corresponding: true
    orcid: 0000-0002-0186-4340
  - name: Marcel Hurler
    affiliation: 1
  - name: Alexander Van Craen
    affiliation: 1 
    orcid: 0000-0002-3336-7226
  - name: Michael Obersteiner
    affiliation: 2
  - name: Dirk Pflüger
    affiliation: 1 
    orcid: 0000-0002-4360-0212
affiliations:
 - name: University of Stuttgart, Scientific Computing, Stuttgart, Germany
   index: 1
 - name: Technical University of Munich, Chair of Scientific Computing, Munich, Germany
   index: 2
nocite: |
  @obersteiner_spatially_2021
date: 16 November 2023
bibliography: paper.bib

---

# Summary

`DisCoTec` is a C++ framework for the sparse grid combination technique [@griebelCombinationTechniqueSolution1992],
designed for massively parallel settings.
It is implemented with shared-memory parallelism via OpenMP and
distributed-memory parallelism via MPI, and is intended to be used in
conjunction with existing simulation codes.
For simulation codes that can handle nested structured grids, little to no
adaptation work is needed for use with the `DisCoTec` framework.
The combination technique with `DisCoTec` demonstrates its superiority in
memory-per-precision for higher-dimensional time-dependent simulations, such
as high-fidelity plasma turbulence simulations
in four to six dimensions
and even for simulations in two dimensions, improvements can be observed [@pollingerStableMassconservingSparse2023].

A central part of the combination technique at scale is the transformation of
grid coefficients into a multi-scale basis.
`DisCoTec` provides a selection of three different lifting wavelets for this
purpose: hierachical hat basis, biorthogonal, and fullweighting basis.
In addition, any code that can operate on nested structured grids can benefit
from the model order reduction provided by the underlying sparse grid approach
used by `DisCoTec`, without requiring any multi-scale operations.
An additional feature of `DisCoTec` is the possibility of performing
widely-distributed simulations of higher-dimensional problems, where multiple
High-Performance Computing (HPC) systems collaborate to solve a joint simulation,
as demonstrated in [@pollingerRealizingJointExtremeScale2024].
Thus, `DisCoTec` can leverage the compute power and main memory of multiple HPC
systems, with comparatively low and manageable transfer costs due to the
combination technique.

# Statement of need

Higher-dimensional problems (by which we mean more than three space
dimensions and one time dimension) quickly require infeasible amounts of
computational resources such as memory and core-hours as the problem size
increases---they are haunted by the so-called 'curse of dimensionality'.
An example of this are high-fidelity plasma turbulence simulations in the field
of confined fusion research.
Currently employed approaches to this problem include dimensionally-reduced models,
such as gyrokinetics [@brizardFoundationsNonlinearGyrokinetic2007]
(which may not always be applicable),
particle-in-cell methods (which suffer from inherent noise [@verboncoeurParticleSimulationPlasmas2005]),
and restricting computations to a very limited resolution.
A further---still developing but very promising---approach to the problem are
low-rank methods [@einkemmerMassMomentumEnergy2021].
Multi-scale (hierarchical) methods, such as the sparse grid combination
technique (CT) that `DisCoTec` employs,
provide an alternative approach to addressing the curse of dimensionality by
considering only those resolutions where the highest amount of information is expected
[@bungartzSparseGrids2004].
While some implementations of the CT are
available, there is currently no other implementation for
parallel simulations that require distributed computing.
`DisCoTec` is a C++ framework for massively-parallel time-dependent problems
with the CT, which fills this gap.


# Methods: Sparse grid combination technique and implementation

The sparse grid combination technique (with time-stepping) is a multi-scale
approach for solving higher-dimensional problems.
Instead of solving the problem on one grid that is very finely resolved in all dimensions,
the problem is solved on the so-called 'component grids' which are all rather
coarsely resolved---each of them differently in the different dimensions.
For instance, the following schematic shows a two-dimensional combination scheme, 
consisting of seven component grids.

![Combination scheme in two dimensions with $\vec{l}_{min} = (2,1)$ and $\vec{l}_{max} = (5,4)$, periodic boundary conditions. Figure first published in  [@pollingerStableMassconservingHighdimensional2024]. \label{fig:combischeme-2d}](gfx/combischeme-2d.pdf)

By updating each other's information throughout the simulation, the component grids
still obtain an accurate solution of the overall problem [@griebelCombinationTechniqueSolution1992].
This is enabled by an intermedate transformation into a multi-scale (hierarchical)
basis, and application of the combination formula
$$ f^{(\text{s})} = \sum_{\vec{\ell} \in \mathcal{I} } c_{\vec{\ell}} f_{\vec{\ell}} $$
where $f^{(\text{s})}$ is the sparse grid approximation, and $f_{\vec{\ell}}$ are
the component grid functions.
The set of all used levels $\vec{\ell}$ is often called a combination scheme $\mathcal{I}$.
In \autoref{fig:combischeme-2d}, the coefficients $c_{\vec{\ell}}$ are $-1$ for
the coarser component grids (red background) and $1$ for the finer component grids
(orange background).
In summary, each of the grids will run (one or more) time steps of the simulation,
then exchange information with the other grids, and repeat this process until
the simulation is finished.

`DisCoTec` provides the necessary infrastructure for the combination technique
with a black-box approach, enabling massive parallelism---suitable for existing
distributed solvers that use structured grids.
An important feature is the usage of 'process groups', where multiple MPI ranks
will collaborate on a set of component grids, and the solver's existing
parallelism can be re-used.
The process groups are displayed as $pg_i$ in \autoref{fig:discotec-ranks}.

![`DisCoTec` process groups: Each black square denotes one MPI rank. The ranks are grouped into the so-called 'process groups'. Distributed operations in `DisCoTec` require either communication in the process group, or perpendicular to it---there is no need for global communication or synchronization, which avoids a major scaling bottleneck. The manager rank is optional. Figure first published in  [@pollingerStableMassconservingHighdimensional2024]. \label{fig:discotec-ranks}](gfx/discotec-ranks.pdf)

In addition, the number of process groups can be increased to leverage the
combination technique's embarrassing parallelism in the solver time steps.
In \autoref{fig:discotec-ranks}, this would be equivalent to adding more and more 
process groups to the right.

Using `DisCoTec`, kinetic simulations were demonstrated to scale up to hundreds
of thousands of CPU cores [@pollingerStableMassconservingHighdimensional2024].
By putting a special focus on saving memory, most of the memory is available for
use by the black-box solver, even at high core counts.
In addition, OpenMP parallelism can be used to further increase parallelism while
being more lightweight than MPI in terms of memory.

Through highly parallel I/O operations, `DisCoTec` can be used to perform
simulations on multiple High Performance Computing (HPC) systems simultaneously,
if there exists a tool for
sufficiently fast file transfer between the systems [@pollingerStableMassconservingHighdimensional2024].
The `DisCoTec` repository contains example scripts and documentation for
utilizing UFTP as an example of a transfer tool, but the approach is not limited
to UFTP.

`DisCoTec` provides a conveniently automated way of installation using a
[`spack` package](https://github.com/spack/spack/blob/develop/var/spack/repos/builtin/packages/discotec/package.py)
[@gamblinSpackPackageManager2015],
which can be used to install `DisCoTec` and its whole dependency tree
in an automated manner optimized for HPC hardware.

# State of the field

Besides `DisCoTec` there exist other frameworks that allow the usage of
sparse grids and the combination technique.
We will give a brief overview and outline the differences and
application areas of the codes.

The C++ code `SG++`[@SGppSGpp2024] provides a direct interface to
sparse grids and applying them to a variety of different tasks such as interpolation,
quadrature, optimization, PDEs,  regression, and classification.
With the help of wrappers, the framework can be used from various other programming
languages such as Python and Matlab.
The code targets direct implementations within sparse grids and provides a basic
implementation of the combination technique.
Although offering parallelization for some of the tasks, the code
mainly targets single-node computations.

The `Sparse Grids Matlab Kit`[@tamelliniLorenzotamelliniSparsegridsmatlabkit2024]
by Piazzola and Tamellini was originally designed for teaching purposes and
uncertainty quantification with the combination technique [@piazzolaSparseGridsMatlab2022].
It offers a user friendly MATLAB interface for the combination technique.
In addition, dimensional adaptivity is available for nested and non-nested sequences
of component grid collocation points.
The code is designed for usage on a single node which limits the parallelism
to shared memory.

The `sparseSpACE` [@obersteinerSparseSpACESparseGrid2023] project offers
different variants of the combination technique including a spatially adaptive
combination technique.
It provides implementations for various applications such as numerical integration,
interpolation, uncertainty quantification, sparse grid density estimation (for
classification and clustering), regression, and PDE calculations.
The code is completely written in Python and is mostly sequential.
The main novelty of this project is the possibility to add spatial adaptivity
to the combination technique.

This demonstrates that there exist multiple codes for sparse grids and the
combination technique.
However, `DisCoTec` is the only code that offers distributed parallelization with
the combination technique and has demonstrated that it can scale up to full
supercomputers and beyond.
In addition, `DisCoTec` uses the most sophisticated approach to utilize the
combination technique with time-dependent PDEs by employing recombinations,
which increases the overall numerical accuracy [@pollingerStableMassconservingSparse2023].

# Acknowledgements

We acknowledge contributions from Mario Heene, Christoph Kowitz, Alfredo Parra
Hinojosa, Johannes Rentrop, Keerthi Gaddameedi, Marvin Dostal,
Marcel Breyer, Christoph Niethammer, Philipp Offenhäuser,
and support from HLRS, LRZ, JSC, and NHR@FAU, where we would like to highlight
the long-standing support by Martin Bernreuther and Martin Ohlerich in particular.

# References
