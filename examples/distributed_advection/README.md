# Advection Example

This example illustrates how to use DisCoTec to run an advection solver
wrapped in `TaskAdvection`. It moves a Gaussian diagonally through a domain
of arbitrary dimensionality.

Details on the initial condition and the solver can be found in
[this preprint](http://arxiv.org/abs/2209.14064).

You can adapt the example by editing the `ctparam` file, exactly like in
[the basic example](../combi_example/README.md).

You can run the singular example by

```bash
mpirun -n $N ./advection_example
```

with the number of MPI processes `$N`
as `ngroup` * `nprocs` + 1 (for the manager process)
-> `3` in the current parameter files.

In order to write `hdf5` output, DisCoTec needs to be compiled with `DISCOTEC_USE_HIGHFIVE`.
In addition, the output will be written in raw format if a folder named `out`
was created before.

The raw data can be visualized, for example, with

```bash
../../tools/plot_raw_slice.py out/solution_hat_2D_3-10_10_0.raw
```

which will put a corresponding `.pdf` file in the `out/` folder.
