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
  - name: Alexander Van Craen
    affiliation: 1 
    orcid: 0000-0002-3336-7226
  - name: Dirk Pflüger
    affiliation: 1 
    orcid: 0000-0002-4360-0212
affiliations:
 - name: University of Stuttgart, Scientific Computing, Stuttgart, Germany
   index: 1
nocite: |
  @obersteiner_spatially_2021
date: 16 November 2023
bibliography: paper.bib

---

# Summary

`DisCoTec` is a C++ framework for the sparse grid combination technique, targeting massively parallel settings.
It provides shared-memory parallelism via OpenMP and distributed-memory parallelism via MPI, 
and is designed to be used in combination with existing simulation codes.

Any code that can operated on nested structured grids can employ the model order reduction 
provided by the underlying sparse grid approach, without considering any multi-scale operations; this part is provided by DisCoTec.
Although already 2D applications can see significant benefits, the higher-dimensional (4- to 6-dimensional) 
grids employed in high-fidelity plasma simulations benefit even more from the combination technique [@pollingerStableMassconservingSparse2023].

Further features include the widely-distributed simulation of higher-dimensional problems,
in which multiple HPC systems cooperate to solve a joint simulation [@pollingerLeveragingComputePower2023].
Thus, `DisCoTec` can leverage the compute power and main memory of multiple HPC systems.
The transfer cost is relatively low due to the multi-scale approach in the combination technique 
-- much less than with a traditional domain decomposition.

# Statement of need

Higher-dimensional problems (by which we typically mean more than three space 
dimensions and one time dimension) quickly require infeasible amounts of computational resources 
such as memory and core-h---they are haunted by the so-called curse of dimensionality.
An example of this are high-fidelity plasma simulations in the field of confined fusion research.
Current approaches to this problem include dimensionally-reduced models (which may not always be applicable),
and restricting oneself to a very limited resolution.
Multi-scale (hierarchical) methods, such as the sparse grid combination technique, 
provide an alternative approach to addressing the curse of dimensionality.
While some implementations of the sparse grid combination technique are available in the context of UQ,
there is currently no implementation for parallel simulations that require distributed computing---apart from `DisCoTec`.

`DisCoTec` is a C++ framework for the sparse grid combination technique.
Targeted at HPC systems, it is used for parallel simulations,
drawing on distributed-memory parallelism via MPI [@heeneMassivelyParallelCombination2018] 
and shared-memory parallelism via OpenMP.
It is designed to be used in combination with existing simulation codes,
which can be used with `DisCoTec` in a black-box fashion.


# Method: Sparse grid combination technique

The sparse grid combination technique (with time-stepping) is a multi-scale approach for solving higher-dimensional problems.
Instead of solving the problem on one grid that is very finely resolved in all dimensions,
the problem is solved on the so-called component grids which are all rather coarsely resolved --
each of them differently in the different dimensions.

![Combination scheme in two dimensions with $\vec{l}_{min} = (1,1)$ and $\vec{l}_{max} = (3,3)$, periodic boundary conditions](gfx/combi-2d-small-periodic.pdf)

By updating each other's information throughout the simulation, the component grids
still obtain an accurate solution of the overall problem. 
This is enabled by an intermedate transformation into a multi-scale (hierarchical) basis, and application of the combination formula
$$ f^{(\text{s})} = \sum_{\vect{l} \in \mathcal{I} } c_{\vect{l}} f_{\vect{l}} $$
where $f^{(\text{s})}$ is the sparse grid approximation, and $f_{\vect{l}}$ are the component grid functions.
In summary, each of the grids will run (one or more) time steps of the simulation, 
then exchange information with the other grids, and repeat this process until the simulation is finished.


# Implementation

`DisCoTec` provides the necessary infrastructure for the combination technique with a black-box approach, 
enabling massive parallelism---suitable for existing solvers that use MPI and structured grids.
An important feature is the usage of process groups, where multiple MPI ranks will collaborate on a set of component grids, 
and the solver's existing parallelism can be re-used.
In addition, the number of process groups can be increased to leverage the 
combination technique's embarrassing parallelism in the solver time steps.
![DisCoTec process groups](gfx/discotec-ranks.pdf)

Using DisCoTec, kinetic simulations could be demonstrated to scale up to hundreds of thousands of cores.
By putting a special focus on saving memory, most of the memory is available for use by the black-box solver, even at high core counts. 
In addition, OpenMP parallelism can be used to further increase parallelism and decrease main memory usage.

Through highly parallel I/O operations, `DisCoTec` can be used to perform simulations on multiple HPC systems simultaneously, 
if there exists a tool for fast file transfer between the systems[@pollingerLeveragingComputePower2023].
The communication between different systems is enabled by file transfer. 
The software repository contains example scripts and documentation for utilizing UFTP as an example of a transfer tool,
but the approach is not limited to UFTP.


# Acknowledgements

We acknowledge contributions from Mario Heene, Christoph Kowitz, Alfredo Parra Hinojosa, Michael Obersteiner, 
Marcel Hurler, Johannes Rentrop, Keerthi Gaddameedi, Marvin Dostal, Marcel Breyer, Christoph Niethammer, Philipp Offenhäuser, 
and support from HLRS, LRZ, JSC, and NHR@FAU, where we would like to particularly highlight the long-standing support by Martin Bernreuther and Martin Ohlerich.

# References
