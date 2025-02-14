# Third-level / Widely-distributed Combination Technique Example

There are two main different ways of doing the widely-distributed combination technique:

1. Using a direct connection (TCP or TCP-over-SSH) from the manager rank
   to a separate process which then communicates the data to the other system
2. Using a loosely-coupled file exchange mechanism to the other system

In this folder, you find both flavors of widely-distributed combination technique.
See also the [README for third-level manager / file exchange scripts](../../third_level_manager/README.md)
on how to configure the non-DisCoTec component in such simulations.

We are going to use different `ctparam` files for the different systems and setups;
when you adapt them you always need to make sure they match each other:
they can only differ in the input files they process, and the number of process groups;
they must have different system numbers.

## 1. Connection-based

For testing on a single system, we need to use a total of three different terminals.

First, set up the third-level manager process:

```bash
cd ../../third_level_manager
./thirdLevelManager example.ini
```

This will have the manager listen at the specified port.
If something goes wrong, you may have to switch to another port in this ini file
and the corresponding `ctparam` files (e.g. `9999` -> `9998`).

Then, in two other shells, run

```bash
mpirun -n 2 ./distributed_third_level_combi_example ctparam_tl_system0
```

and

```bash
mpirun -n 2 ./distributed_third_level_combi_example ctparam_tl_system1
```

each invocation of `./distributed_third_level_combi_example` will write their
`timers-$s.json` file (with `$s` the system number).

### Sample command line outputs

For the `thirdLevelManager`:

```bash
$ ./thirdLevelManager example.ini
thirdLevelManager running on separate system
broker = pollinta-ThinkPad-P16s-Gen-3
Third level manager running on port: 9999
 Waiting for system (1/2)
 Waiting for system (2/2)
All systems connected successfully
Simulation took:                     58841.9sec
Num combinations:                    10
Total transfer during combination:   2621460B
Transfer per combination:            262146B
Total transfer during size exchange: 482B
```

For one system:

```bash
$ mpirun -n 2 ./distributed_third_level_combi_example ctparam_tl_system1 
Using third-level parallelism
broker running on same system
Using third-level parallelism
broker running on same system
9 tasks to distribute.
manager = pollinta-ThinkPad-P16s-Gen-3
lmin = [3 3 ]
lmax = [10 10 ]
manager: generated parameters
manager = pollinta-ThinkPad-P16s-Gen-3
Connecting to third level manager at host localhost on port 9999
Connected.
manager: updated parameters in 0 seconds
manager: ran solver in 3.842 seconds, of which SG init were 0
manager: unify sparse grid data structures w/ remote
manager: unified SG in 0 seconds
combination 1 took: 0.134 seconds
write 1 took: 0 seconds
read 1 took: 0 seconds
calculation 1 took: 4.44 seconds
combination 2 took: 0.013 seconds
[...]
combination 9 took: 0.017 seconds
write 9 took: 0 seconds
read 9 took: 0 seconds
calculation 9 took: 6.227 seconds
manager exit
```

## 2. File-based

For testing on a single system, we can just run the executable in two different
terminals, but in the same folder:

```bash
mpirun -n 6 ./distributed_third_level_workers_only ctparam_2D
```

and

```bash
mpirun -n 1 ./distributed_third_level_workers_only ctparam_2D_system1
```

(in case you get an error due to missing interpolation coordinates file,
re-compile with HighFive or obtain such a file from the experiment's data repo
which you are trying to replicate.)

This should write a set of output files to your working directory:

* the interpolated values are written as `worker_interpolated_$s_values_$t.h5`
  with `$s` the system number and `$t` the time step number.
  They can be combined into the total interpolated values with, e.g.,
  `$ ../../tools/hdf5_reduce_datasets.py worker_interpolated_0_values_2.h5 worker_interpolated_1_values_2.h5`,
  which writes a file `worker_interpolated_0_values_2_and_worker_interpolated_1_values_2.h5`
  to be further processed as in [the workers-only example](../combi_workers_only/README.md)
* the timers are written as slices for each system and for each combination interval
  (files named `$t_stats_worker_$s_group$g.json`) as well as once per group
  with all measurements at the end (files named `timers_system$s_group$g.json`).
  
  To get the average timings e.g. for system 0, run

  ```bash
  ../../tools/evaluate.py --input_files *_stats_worker_0_group*.json
  ```

  to use the slices (useful if the simulation aborts prematurely), or

  ```bash
  ../../tools/evaluate.py --input_files timers_system0_group*.json
  ```

  to use the final timers files. This tool writes output statistics to `.csv` files
  in your working folder, and prints two lines that can also be used for csv output.
