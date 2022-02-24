## Selalib Example

* runs the Semi-Lagrangian BSL code
* public part available here https://github.com/selalib/selalib
* "our" currently used branch is https://gitlab.mpcdf.mpg.de/clapp/selalib/-/tree/exahd_communicator

### Build Selalib Dependency
* clone and check out right version (see above)
* depends on `fftw` and `hdf5+fortran@1.10.5` or higher (in spack notation)
* create build directory, load dependencies
* in the build directory, run
```
cmake -DCMAKE_BUILD_TYPE=Debug -DOPENMP_ENABLED=ON -DHDF5_PARALLEL_ENABLED=ON $SELALIB_DIR
make test_cpp
make sll_m_sim_bsl_vp_3d3v_cart_dd_slim_movingB
make sll_m_sim_bsl_vp_3d3v_cart_dd_slim_movingB_interface
```
where `test_cpp` can be used to test the gerenal C/Fortran interface, `sll_m_sim_bsl_vp_3d3v_cart_dd_slim_movingB` is the mononlithic solver for our test case, and `sll_m_sim_bsl_vp_3d3v_cart_dd_slim_movingB_interface` builds the libraries needed for this `selalib_distributed` example. (This may take some time and usually fails if tried in parallel with `make -j`).


### Build Selalib+DisCoTec example
* update Makefile.template in this folder (could not figure out an automatic way, please update mpi and selalib paths in there)
* install DisCoTec as described in the general README; this will generate a Makefile from Makefile.template in this folder
* run `make` here


### Run Selalib+DisCoTec example
* update combination technique parameters in `ctparam`, selalib parameters in `template/param.nml` (be careful with `haveDiagnosticsTask` -- it interpolates onto `lmax` to do the diagnostics, which can be extremely costly)
* load dependencies and run with mpi, like shown in `sbatch.sh` -- the number of MPI tasks has to match the process group number and size in `ctparam`

### Postprocessing for Selalib+DisCoTec example
* run `./combine_selalib_diagnostics ctparam` or update paths in `./plot_landau.py` and run (which also invokes`./combine_selalib_diagnostics` but then plots the quantities of interest)
* analyze the sparse grid "surplusses" (maximum hierarchical coefficients per level) by running `../../tools/visualize_sg_minmax.py plot.dat_selalib_sg_${TIMESTEP}_0.txt`. It shows a correlation matrix and plots projections of maxima on the different level combinations.