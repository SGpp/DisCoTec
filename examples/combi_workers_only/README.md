# Example Workers-Only

This example

* shows how the simulation can be set up without a manager rank (and fault tolerance)
* illustrates how to use json-scheme reading: The file `scheme_6.06_GiB.json`
  was generated with [the combischeme utilities](https://github.com/SGpp/DisCoTec-combischeme-utilities).
  The file does not assign process group numbers, so the tasks will be distributed
  round-robin to all process groups.
  The GiB in the file name is the memory volume that will be needed for
  the component grids -> expect more for sparse grid and MPI data structures
* uses different-size chunks for chunked reduction algorithms
  (it is large and runs relatively long so you can see a difference in run time)
* also shows how to evaluate the Monte-Carlo error integral (on 100,000 points)

## Set up and run

* adapt `ctparam`: the `chunkSize` parameter is the buffer size in MiB
  used for the chunked reduction;
  Larger values need more memory but speed up the combination.
* call (with default parameters): `mpirun -n 4 ./combi_example_worker_only`
  
  sample output (observe from the time stamps that this took around 50min):

    ```bash
    $ mpirun -n 4 ./combi_example_worker_only
    MPI: 2 groups with 2 ranks each without world manager
    [20250127T104016] initialized communicators
    [20250127T104016]  Process group 0 will run 14 of 28 tasks.
    Combination scheme DOF : 411041792 i.e. 3.28833 GB 
    [20250127T104016]  Process group 1 will run 14 of 28 tasks.
    Combination scheme DOF : 402653184 i.e. 3.22123 GB 
    [20250127T104016] generated parameters
    [20250127T104343] worker: initialized tasks 
    [20250127T104411] worker: initialized SG
    [20250127T104411] group 1: set sparse grid sizes, will allocate 234.881 MB (but only 16.7772 MB at once)
    [20250127T104411] initialization took: 234 seconds
    [20250127T104411] start simulation loop
    [20250127T104411] group 0: set sparse grid sizes, will allocate 234.881 MB (but only 16.7772 MB at once)
    [20250127T105258] calculation 0 took: 527.57 seconds
    [20250127T105359] interpolation 0 took: 61.179 seconds
    [20250127T105451] combination 0 took: 51.659 seconds
    [20250127T110328] calculation 1 took: 516.501 seconds
    [20250127T110502] interpolation 1 took: 94.305 seconds
    [20250127T110554] combination 1 took: 51.673 seconds
    [20250127T111430] calculation 2 took: 515.722 seconds
    [20250127T111538] interpolation 2 took: 67.719 seconds
    [20250127T111630] combination 2 took: 52.007 seconds
    [20250127T112508] last calculation 3 took: 518.559 seconds
    [20250127T112617] last interpolation 3 took: 68.31 seconds
    ```

## Inspect results

* The example will write interpolated values as well as their interpolation
  coordinates as `.hdf5` files. They can be postprocessed e.g. with
  `../../tools/hdf5_interpolation_norms.py`:

  ```bash
  $ ../../tools/hdf5_interpolation_norms.py --solution=advection --coordinates=interpolation_coords_6D_100000.h5 worker_interpolated_values_3.h5 
  ic| values_file: 'worker_interpolated_values_3.h5'
      reference_solution: ['advection']
  ic| values_dataset: <HDF5 dataset "interpolated_3": shape (100000,), type "<f8">
  ic| simulation_time: np.float64(0.004)
  ic| coordinates_dataset: <HDF5 dataset "only": shape (100000, 6), type "<f8">
  Monte Carlo error integral of order inf : 0.006725605769273502
  Monte Carlo error integral of order 1 : 0.00012767923912647587
  Monte Carlo error integral of order 2 : 1.122707832581847e-06
  ```
  
  This output shows the maximum norm of the error, as well as the L_1 and L_2
  error norms, w.r.t. the analytical (Gaussian) solution.

* The `timers.json` can be visualized just like in
  [the `combi_example`](../combi_example/README.md)
