---
title: 'DisCoTec: MPI-based code for distributed HPC simulations with the sparse grid combination technique'
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

Higher-dimensional problems (by which we typically mean more than three space 
dimensions and one time dimension) quickly require huge amounts of computational resources such as memory and core-h.
An example of this are plasma simulations in the field of confined fusion research.
The sparse grid combination technique addresses this problem:
Instead of solving the problem on one grid that is very finely resolved in all dimensions,
the problem is solved on a combination of grids which are all rather coarsely resolved --
each of them differently in the different dimensions.
By updating each other's information throughout the simulation, the grids still solve the problem accurately.


# Statement of need

`DisCoTec` is a C++ framework for the sparse grid combination technique.
Targeted at HPC systems, it is used for parallel simulations,
drawing on distributed-memory parallelism via MPI [@heeneMassivelyParallelCombination2018] and shared-memory parallelism via OpenMP.
It is designed to be used in combination with existing simulation codes,
which can be used with `DisCoTec` in a black-box fashion.

A further application includes the widely-distributed simulation of higher-dimensional problems,
in which multiple HPC systems cooperate to solve a joint simulation [@pollingerLeveragingComputePower2023].
The transfer cost is relatively low due to the multi-scale approach in the combination technique 
-- much less than with a traditional domain-decomposition.   
This feature is enabled by file transfer through tools like UFTP.

Basically, any code that can operated on nested structured grids can employ the model order reduction 
provided by the underlying sparse grid approach without considering any multiscale operations; this is provided by DisCoTec.
Although already 2D applications can see significant benefits, the higher-dimensional (4- to 6-dimensional) 
grids employed in high-fidelity plasma simulations benefit even more [@pollingerStableMassconservingSparse2023].


# Acknowledgements

We acknowledge contributions from Mario Heene, Christoph Kowitz, Alfredo Parra Hinojosa, Michael Obersteiner, Marcel Hurler, Johannes Rentrop, Keerthi Gaddameedi, Marvin Dostal, Marcel Breyer, Christoph Niethammer, Philipp Offenhäuser, and support from HLRS, LRZ, JSC, and NHR@FAU, where we would like to highlight the long-standing support by Martin Bernreuther and Martin Ohlerich.

# References
