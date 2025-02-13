# Example for Fault Tolerance

In addition to the parameters in the [basic combi example](../combi_example/README.md),
you can set the following to control the behavior of the random faults:

```ini
[faults]
num_faults = 2
iteration_faults = 1 2
global_rank_faults = 3 1
```

Here `num_faults` specifies the number of faults that should occur. After that the combination iteration of the fault is specified for each fault via `iteration_faults` and in addition the MPI rank that fails is specified in `global_rank_faults`. In the above example, we specify two faults: rank 3 fails in the first combination interval and rank 1 in the second combination interval. It is important to specify the same amount of iterations and ranks as `num_faults indicates.

It is also possible to use a random distribution for the failure of ranks by using a negative value for `num_faults`. This implies that faults will occur non-deterministic. In that case the absolute value of `num_faults` is used to specify the lambda parameter of the Weibull distribution. A small `lambda` value relates to a higher failure probability (it is related to the mean-time between failures). See [1] for more information on this random process and the Weibull distribution. The values of `ìteration_fault` and `global_rank_faults` are both ignored for negative values of `num_faults`. An example with a random failure behaviour and a `lambda = 10'000` could be:

```ini
[faults]
num_faults = -10000
iteration_faults = 1 2
global_rank_faults = 3 1
```

Important: To use the fault-tolerant mechanisms it is necessary to compile the code using the `-DENABLEFT` flag!

Using the default settings, run with

```bash
mpirun -n 5 ./combi_example_faults
```

[1] Obersteiner, Michael; Parra Hinojosa, Alfredo; Heene, Mario; Bungartz, Hans-Joachim; Pflüger, Dirk: A Highly Scalable, Algorithm-Based Fault-Tolerant Solver for Gyrokinetic Plasma Simulations. ScalA '17: Proceedings of the 8th Workshop on Latest Advances in Scalable Algorithms for Large-Scale Systems, 2017 (https://dl.acm.org/doi/10.1145/3148226.3148229)
