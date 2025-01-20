# Basic Example

This example illustrates the flow of information in a DisCoTec simulation with
a very simple task: A time-dependent cosine function is interpolated on every
combination technique component grid.

You can adapt the example by editing the `ctparam` file:

```ini
[ct]
dim = 2
# minimum and maximum level, result in different numbers of grids
# and different resolutions
lmin = 3 3 
lmax = 12 12
# level at which the solution will be interpolated and written
leval = 5 5 
# parallelization vector, needs to contain powers of 2
p = 1 2
# number of combinations
ncombi = 10

# parameters for the task
[application]
# time step size
dt = 1e-3
# number of task time steps per combination
nsteps = 100

[manager]
# number of groups (must be smaller than the number of component grids)
ngroup = 2
# number of processes per group, must be the product of p's entries
nprocs = 2

```

You can then run the example by

```bash
mpirun -n $N ./combi_example
```

with the number of MPI processes `$N`
as `ngroup` * `nprocs` + 1 (for the manager process).

The example will write the combined solution after every time step
to the `out/` folder.
