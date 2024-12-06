# Parallelism in DisCoTec

The DisCoTec framework can work with existing MPI-parallelized PDE solver codes
operating on structured grids.
In addition to the parallelism provided by the PDE solver, it adds the combination
technique's parallelism.
This is achieved mainly through *process groups* (pgs):
`MPI_COMM_WORLD` is subdivided into equal-sized process groups
(and optionally, a manager rank).

![schematic of MPI ranks in DisCoTec](../gfx/discotec-ranks.svg)

The image describes the two ways of scaling up:
One can either increase the size or the number of process groups.
Figure originally published in (Pollinger [2024](https://elib.uni-stuttgart.de/handle/11682/14229)).

All component grids are distributed to process groups (either statically, or
dynamically through the manager rank).
Each process group will be responsible for a set of component grids and the
solver instance belonging to each.
Within a process group and each grid, one can observe the exact same domain 
decomposition parallelism as when using the solver directly.
(The domain decomposition needs to be the same for all solver instances.)

During a time step, the PDE solver step applies one or multiple time step updates
to the values in each component grid.
During the PDE solver time step and most of the combination step, MPI communication
only happens within the process groups.
Conversely, for the sparse grid reduction using the combination coefficients $c_\vec{\ell}^c$,
MPI communication only happens between a rank and its colleagues in the other
process groups, e.g., rank 0 in group 0 will only talk to rank 0 in all other groups.
Thus, major bottlenecks arising from global communication can be avoided altogether.

In practice, the parallelization in DisCoTec is set at runtime, by calling

```cpp
  size_t ngroup = 8;
  size_t nprocs = 32;
  bool withWorldManager = false;
  combigrid::theMPISystem()->initWorldReusable(MPI_COMM_WORLD, ngroup, nprocs, withWorldManager);
```

for example.
If needed, the exact distribution of coordinates to ranks is set in the
`CombiParameters` class.

Combining the two ways of scaling up, DisCoTec's scalability was demonstrated on
several machines, with the experiments comprising up to 524288 cores:

![timings for advection solver step on HAWK at various
parallelizations](../gfx/times-solver-on-hawk.svg)
![timings for combination step on
HAWK at various parallelizations](../gfx/times-combination-on-hawk.svg)

We see the timings (in seconds) for the advection PDE solver step and the
combination step, respectively.
This weak scaling experiment used four OpenMP threads per rank, and starts with
one pg of four processes in the upper left corner.
The largest parallelization is 64 pgs of 2048 processes each.
Figure originally published in (Pollinger [2024](https://elib.uni-stuttgart.de/handle/11682/14229)).

This image also describes the challenges in large-scale experiments with DisCoTec:
Generally, the process group size can be chosen so that it matches the solver
and the desired resolutions well.
But if the process groups become too large, the MPI communication of the multiscale
transform starts to dominate the combination time.
Conversely, the number of process groups can be freely chosen, as long as there
is sufficient work for each group to do (one grid per group at a bare minimum).
If there are too many pgs, the combination reduction will dominate the
combination time.
However, the times required for the PDE solver stay relatively constant;
they are determined by the PDE solver's own scaling and the load balancing quality.

There are only few codes that allow weak scaling up to this problem size:
a size that uses most of the available main memory of the entire system.
