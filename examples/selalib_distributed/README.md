# Selalib Example

* runs the Semi-Lagrangian BSL code
* public part available here <https://github.com/selalib/selalib>
* "our" currently used branch is <https://gitlab.mpcdf.mpg.de/clapp/selalib/-/tree/exahd_communicator>

## Build Selalib Dependency

* clone and check out right version (see above)
* depends on `fftw` and `hdf5+fortran@1.10.5` or higher (in spack notation)
* create build directory, load dependencies
* in the build directory, run
  
```sh
cmake -DCMAKE_BUILD_TYPE=Release -DOPENMP_ENABLED=ON -DHDF5_PARALLEL_ENABLED=ON  -DUSE_FMEMPOOOL=OFF -DCMAKE_INSTALL_PREFIX=$(pwd)/install $SELALIB_DIR
make test_cpp_interface
make sll_m_sim_bsl_vp_3d3v_cart_dd_slim_interface
make sim_bsl_vp_3d3v_cart_dd_slim
make install
```

where `test_cpp_interface` can be used to test the gerenal C/Fortran interface,
`sim_bsl_vp_3d3v_cart_dd_slim` is the mononlithic PDE solver for our test case, and
`sll_m_sim_bsl_vp_3d3v_cart_dd_slim_interface` builds the libraries needed for
this `selalib_distributed` example.
(This may take some time and usually fails if tried in parallel with `make -j`).

### Build Selalib+DisCoTec example

* run DisCoTec's cmake with `-DDISCOTEC_WITH_SELALIB=1 -DSELALIB_DIR=` flags
  (`SELALIB_DIR` should be set to selalib's `CMAKE_INSTALL_PREFIX` folder, appended
  with `/cmake`, or wherever you find `SELALIBConfig.cmake` accompanied by
  `libselalib.a` and `libsll_m_sim_bsl_vp_3d3v_cart_dd_slim_interface.a`)
* then `make selalib_distributed`

### Run Selalib+DisCoTec example

* update combination technique parameters in `ctparam`, selalib parameters in
  `template/param.nml` (be careful with `haveDiagnosticsTask` -- it interpolates
  onto `lmax` to do the diagnostics, which can be extremely costly)
* load dependencies and run with mpi, like shown in `sbatch.sh` -- the number of
  MPI tasks has to match the process group number and size in `ctparam`.
  Without modification, this should look similar to
  `mpiexec -n 9 ./selalib_distributed` or
  `mpiexec -n 8 ./selalib_distributed_workers_only` ;
  Make sure that OpenMP runs with 1 process per MPI rank.
  The examples will create a separate output folder called `seladisco$n`
  for every component grid (with index `$n`) in the current directory.
  You should see some output like this for the first time step

```bash
 Entering main loop ... 
 Time   0.010 of   0.010 :: step run time [s] =   0.021
 Leaving main loop.  Main loop run time [s] =      0.021
 Running 6D Vlasov simulation with MPI domain decomposition (slim) ...
 git: main-0-g591da44
 compiler: gfortran 13.2.0 (-std=f2008 -ffree-line-length-none -fstack-arrays -O3 -fPIC  -w  -march=native -fallow-argument-mismatch  -fopenmp)
parallelism :: MPI procs:     8, OMP threads:  18
               process grid:    1   1   1   2   2   2
         process-local grid:    8   8   8   4   4   4
```

  and like this for all remaining time steps

```bash
 Entering main loop ... 
 Time   0.020 of   0.010 :: step run time [s] =   0.048
 Leaving main loop.  Main loop run time [s] =      0.049
```

  for every component grid in the simulation (-> can be a lot of output!)

### Postprocessing for Selalib+DisCoTec example

* run `./combine_selalib_diagnostics ctparam` or
  `./plot_selalib_diagnostics.py --ctparam_to_use ctparam`
  (which also invokes`./combine_selalib_diagnostics`
  but then plots the quantities of interest)
* analyze the sparse grid "surplusses" (maximum hierarchical coefficients per
  level) by running
  `../../tools/visualize_sg_minmax.py plot.dat_selalib_sg_${TIMESTEP}_0.txt`.
  It shows a correlation matrix and plots projections of maxima on the different
  level combinations.
