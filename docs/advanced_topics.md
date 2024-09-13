# Advanced Topics and Further Resources

## General Reading

### Sparse Grids
- H.-J. Bungartz, M. Griebel. ‘Sparse grids’. In: Acta Numerica 13 (May 2004). Publisher: Cambridge University Press, pp. 147–269. url: https://www.cambridge.org/core/journals/acta-numerica/article/sparse-grids/47EA2993DB84C9D231BB96ECB26F615C
- J. Garcke. ‘Sparse Grids in a Nutshell’. In: Sparse Grids and Applications. Ed. by J. Garcke, M. Griebel. Lecture Notes in Computational Science and Engineering. Berlin, Heidelberg: Springer, 2013, pp. 57–80.
- D. Pflüger. Spatially Adaptive Sparse Grids for High-Dimensional Problems. Verlag Dr. Hut, 2010.
- W. Guo, Y. Cheng. ‘A Sparse Grid Discontinuous Galerkin Method for High-Dimensional Transport Equations and Its Application to Kinetic Simulations’. In: SIAM Journal on Scientific Computing 38.6 (Jan. 1, 2016), A3381–A3409. url: https://epubs.siam.org/doi/10.1137/16M1060017
- S. Schnake et al. ‘Sparse-grid discontinuous Galerkin methods for the Vlasov–Poisson–Lenard–Bernstein model‘. In: Journal of Computational Physics, vol. 510, p. 113053, Aug. 2024, doi: 10.1016/j.jcp.2024.113053.
- A. Zeiser. ‘Sparse Grid Time-Discontinuous Galerkin Method with Streamline Diffusion for Transport Equations’. In: Partial Differential Equations and Applications 4.4 (Aug. 1, 2023), p. 38. url: https://doi.org/10.1007/s42985-023-00250-2 .


### Combination Technique
- M. Griebel, M. Schneider, C. Zenger. ‘A combination technique for the solution of sparse grid problems’. In: Proceedings of the IMACS International Symposium on Iterative Methods in Linear Algebra: Brussels, Belgium, 2 - 4 April, 1991. Ed. by P. de Groen, R. Beauwens. North Holland, 1992.
- M. Heene. ‘A Massively Parallel Combination Technique for the Solution of High-Dimensional PDEs’. PhD thesis. Institut für Parallele und Verteilte Systeme der Universität Stuttgart, 2018.


## Variants of the Combination Technique in DisCoTec

### Truncated Combination Technique
DisCoTec by default uses the *truncated combination technique*, which sets a minimum level for all component grids used.

- J. Benk, D. Pflüger. ‘Hybrid parallel solutions of the Black-Scholes PDE with the truncated combination technique’. In: 2012 International Conference on High Performance Computing Simulation (HPCS). 2012 International Conference on High Performance Computing Simulation (HPCS). July 2012.

### Adaptivity

Of the two ways of bringing adaptivity into the combination technique, DisCoTec currently supports only one: 
Static dimensional adaptivity can be achieved by providing the combination scheme as a .json input file.
These files can be generated with https://github.com/SGpp/DisCoTec-combischeme-utilities . 
To get the more advanced combination schemes as discussed in C. Kowitz' dissertation, you may have to adapt the utilities code.
If your application needs dynamic (= during runtime) or spatial adaptivity (= block-structured grids), please get in touch with the developers.

- T. Gerstner, M. Griebel. ‘Dimension–adaptive tensor–product quadrature’. In: Computing 71.1 (2003), pp. 65–87.
- C. Kowitz. ‘Applying the Sparse Grid Combination Technique in Linear Gyrokinetics’. Dissertation. München: Technische Universität München, 2016.
- M. Obersteiner. ‘A spatially adaptive and massively parallel implementation of the fault-tolerant combination technique’. Dissertation. Technische Universität München, 2021. url: https://mediatum.ub.tum.de/doc/1613369/1613369.pdf .


### Load Balancing

The static load balancing strategy used in the [tutorial](./simple_tutorial.md) is solely based on equally distributing the number of DOF, to allow for maximum usage of the main memory.
This ignores that some grids may take a lot longer on certain operations due to their high anisotropy and the resulting cache effects.
The grids can also be assigned to groups through the files generated with https://github.com/SGpp/DisCoTec-combischeme-utilities .
To achieve anisotropy-based and dynamic load balancing, DisCoTec uses `LoadModel`s.

- M. Heene, C. Kowitz, D. Pflüger. ‘Load Balancing for Massively Parallel Computations with the Sparse Grid Combination Technique.’ In: Parallel Computing: Accelerating Computational Science and Engineering (CSE). Ed. by M. Bader, A. Bode, H.-J. Bungartz, M. Gerndt, G. R. Joubert, F. Peters. Vol. 25. Advances in Parallel Computing. 2014, pp. 574–583.
- M. Dostal. ‘Lastbalancierung durch dynamische Aufgaben-Umverteilung mit der Dünngitter-Kombinationstechnik’. Publisher: Universität Stuttgart. BSc Thesis. 2020. url: http://elib.uni-stuttgart.de/handle/11682/10956 .


### Combination Variants

DisCoTec can use several `CombinationVariant`s for the reduction in the combination: `sparseGridReduce`, `subspaceReduce`, `outgroupSparseGridReduce`, `chunkedOutgroupSparseGridReduce`.
Their benefits depend on the assignment of component grids to process groups, which can be competing with load balance.
We currently recommend to use `subspaceReduce` for small problems and `chunkedOutgroupSparseGridReduce` for larger problems, unless you need to collect all data on a single rank (-> `sparseGridReduce`).

- P. Hupp, M. Heene, R. Jacob, D. Pflüger. ‘Global Communication Schemes for the Numerical Solution of High-dimensional PDEs’. In: Parallel Computing 52.C (Feb. 2016), pp. 78–105. url: http://dx.doi.org/10.1016/j.parco.2015.12.006 .
- T. Pollinger. ‘Stable and mass-conserving high-dimensional simulations with the sparse grid combination technique for full HPC systems and beyond.‘ Dissertation, 2024. doi: [10.18419/opus-14210](http://elib.uni-stuttgart.de/handle/11682/14229) .

### Hierarchical Basis Functions / Biorthogonal Wavelets and Boundary Treatment
<!-- TODO -->

- T. Pollinger, J. Rentrop, D. Pflüger, K. Kormann. ‘A Stable and Mass-Conserving Sparse Grid Combination Technique with Biorthogonal Hierarchical Basis Functions for Kinetic Simulations’. In: Journal of Computational Physics (July 7, 2023), p. 112338. url: https://www.sciencedirect.com/science/article/pii/S0021999123004333 .

### Timers
<!-- TODO -->

### File IO

Currently, every DisCoTec application uses a custom variant of `.ini` parameter file parsing from a `ctparam` file. 
This may be standardized in future revisions.

[Interpolation coordinates](../examples/combi_workers_only/combi_example_worker_only.cpp#L141) and other input data can be read from hdf5 files.

The custom combination schemes generated with https://github.com/SGpp/DisCoTec-combischeme-utilities can be used with [`CombiMinMaxSchemeFromFile`](https://github.com/SGpp/DisCoTec/blob/main/examples/combi_workers_only/combi_example_worker_only.cpp#L97).

Timers are written with `.json` files.

The subspace size input and sparse grid output files are only needed for widely-distributed simulations and are discussed below.

#### Application Code I/O

The GENE and SeLaLib examples use a separate folder for each component grid, and generate the input parameter files at the beginning of the main program.
The task then changes the directory at initialization and for the solver update, so that outputs will be placed there.
The derived quantities like energy can then be [combined as a postprocessing step](../examples/selalib_distributed/postprocessing/combine_selalib_diagnostics.cpp#L38).


### Widely-Distributed Simulations 
<!-- TODO -->
-- combination scheme input files
discotec-combischeme-utilities

-- subspace size input files 
--sparse grid output files 

- T. Pollinger, A. Van Craen, C. Niethammer, M. Breyer, D. Pflüger. ‘Leveraging the Compute Power of Two HPC Systems for Higher-Dimensional Grid-Based Simulations with the Widely-Distributed Sparse Grid Combination Technique’. In: SC ’23. Association for Computing Machinery, Nov. 11, 2023. url: https://dl.acm.org/ doi/10.1145/3581784.3607036
  
### Using Solvers Written In Other Programming Languages

Your functions need the same described interface and need to somehow expose it to the C++ compiler.
For Fortran, the SelalibTask shows how the [relevant functions](https://github.com/selalib/selalib/blob/main/simulations/parallel/bsl_vp_3d3v_cart_dd/sll_m_sim_bsl_vp_3d3v_cart_dd_slim_interface.F90) can be [made accessible with `extern "C"`](https://github.com/SGpp/DisCoTec/blob/main/examples/selalib_distributed/src/SelalibTask.hpp).

