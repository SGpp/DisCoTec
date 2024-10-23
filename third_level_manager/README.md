# Widely-distributed Combination Technique Simulations

In this folder, and throughout DisCoTec, there are parts that allow you to
distribute combination technique simulations across HPC compute systems.
In particular, the executables in `/examples/distributed_third_level` give an
overview on the different flavours of DisCoTec applications that are distributed
with the "third level of parallelism".

Notably, to pre-compute the desired splitting between systems, please refer to <https://github.com/SGpp/DisCoTec-combischeme-utilities>.
With these scripts, you can generate the `scheme_*.json` files that the DisCoTec
instance on each system can read.

Generally, there are two ways of performing widely-distributed DisCoTec simulations.

1. The first relies on serializing the sparse grid data through the manager ranks
   on each system, which connect to the `thirdLevelManager` executable (which
   may also run on the compute nodes, but also on an entirely different machine
   that is accessible from the systems).

2. The second relies on file-based combinations, where each system writes their
   partially-reduced sparse grid to the file system, and the data exchage needs
   to be handled separately, for instance with the UFTP scripts in this folder.
   If the other file is received, it is used to further reduce the own sparse grid
   data and computation can continue.

## Third-Level-Manager

The third-level-manager handles the synchronization and data transfer between
the systems.

1. To compile, `make thirdLevelManager`, or build all of DisCoTec on the machine
   where you want the manager to run.

2. The manager takes an `.ini` file as input parameter.
   This folder holds the example file `example.ini` where the number of systems and
   the port on which the manager listens can be adjusted.
   Currently, only 2 systems are supported.

3. To run the manager on its own separate system, use the `run.sh` script to
   execute the manager.
   To run the manager within the MPI call with one of the HPC systems,
   use the `:` syntax for `mpirun` or `mpiexec`, e.g.
   `mpirun -n $nprocs ./distributed_third_level : -n 1 ./thirdLevelManager`
   and set the flag `thirdLevel.brokerOnSameSystem` to true in the `ctparam`
   file on that system.

## File-Based Recombination

The file-based recombination omits the serialization and the `thirdLevelManager`;
It even allows to omit having any manager rank in DisCoTec.

Currently, the subspace sizes need to be precomputed so that DisCoTec knows how
much data it should allocate for the sparse grid reduce mechanism in each process
group, and how large the files containing the conjoint sparse grid should be.
To this end, use the `tools/subspace_writer` ([README here](/tools/subspace_writer/README.md))
to generate the files called `scheme_*.sizes` .

For UFTP based file transfers have a look at the example scripts
`copy_file_from_hawk|ng_to_hawk|no_[parallel].sh` in this folder.  
example call:
`PATHLRZ=$WORK/scenarios/widely USERHAWK=ipvsavcr PATHHAWK=$WORK/widely copy_file_from_hawk_to_ng.sh`
