# SeLaLib Example

* runs the Semi-Lagrangian BSL code
* the example can be used with <https://github.com/selalib/selalib>;
  it uses the SeLaLib target `sll_m_sim_bsl_vp_3d3v_cart_dd_slim_interface`
* the experiments on mass conservation presented
  [here](https://doi.org/10.1016/j.jcp.2023.112338)
  used the branch <https://gitlab.mpcdf.mpg.de/clapp/selalib/-/tree/exahd_communicator> .

## Build SeLaLib Dependency

* clone and check out right version (see above)
* depends on `fftw` and `hdf5+fortran@1.10.5` or higher (in spack notation)
* create build directory, load dependencies
* in the build directory, run
  
```sh
cmake -DCMAKE_BUILD_TYPE=Release -DOPENMP_ENABLED=ON -DHDF5_PARALLEL_ENABLED=ON -DUSE_FMEMPOOOL=OFF -DCMAKE_INSTALL_PREFIX=$(pwd)/install $SELALIB_DIR
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

### Build SeLaLib+DisCoTec example

* run DisCoTec's cmake with `-DDISCOTEC_WITH_SELALIB=1 -DSELALIB_DIR=` flags
  (`SELALIB_DIR` should be set to selalib's `CMAKE_INSTALL_PREFIX` folder, appended
  with `/cmake`, or wherever you find `SELALIBConfig.cmake` accompanied by
  `libselalib.a` and `libsll_m_sim_bsl_vp_3d3v_cart_dd_slim_interface.a`)
* then `make selalib_distributed`

### Run SeLaLib+DisCoTec example

* update combination technique parameters in `ctparam`, selalib parameters in
  `template/param.nml` (be careful with `haveDiagnosticsTask` -- it interpolates
  onto `lmax` to do the diagnostics, which can be extremely costly)
* load dependencies and run with mpi, like shown in `sbatch.sh` -- the number of
  MPI tasks has to match the process group number and size in `ctparam`.
  Make sure that OpenMP runs with 1 process per MPI rank.

### Postprocessing for SeLaLib+DisCoTec example

* run `./combine_selalib_diagnostics ctparam` or update paths in
  `./plot_landau.py` and run (which also invokes`./combine_selalib_diagnostics`
  but then plots the quantities of interest)
* analyze the sparse grid "surplusses" (maximum hierarchical coefficients per
  level) by running
  `../../tools/visualize_sg_minmax.py plot.dat_selalib_sg_${TIMESTEP}_0.txt`.
  It shows a correlation matrix and plots projections of maxima on the different
  level combinations.
