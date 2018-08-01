#!/ bin / bash
export LD_LIBRARY_PATH = / path / to / SGpp / lib / sgpp : $LD_LIBRARY_PATH

                                                               NGROUP =
                             $(grep ngroup ctparam | awk - F "=" '{print $2}') NPROCS =
                                 $(grep nprocs ctparam | awk - F "=" '{print $2}')

                                     mpirun.mpich -
                                 n $(($NGROUP * $NPROCS + 1))./ combi_example
