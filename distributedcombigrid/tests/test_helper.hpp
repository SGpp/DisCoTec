#ifndef TEST_HELPER_HPP
#define TEST_HELPER_HPP

#include <mpi.h>
#include <complex>
#include <boost/test/floating_point_comparison.hpp>

namespace TestHelper{
  static constexpr double tolerance = 1e-12;
  static constexpr double higherTolerance = 1e-5;

  static inline bool checkNumMPIProcsAvailable(int nprocs) {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size >= nprocs;
  }

  static MPI_Comm getComm(int nprocs) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int color = rank < nprocs ? 0 : 1;
    MPI_Comm lcomm;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &lcomm);
    if (rank < nprocs) {
      return lcomm;
    } else {
      return MPI_COMM_NULL;
    }
  }

  static inline int getRank(MPI_Comm comm) {
    int rank;
    MPI_Comm_rank(comm, &rank);
    return rank;
  }
}

#endif  // TEST_HELPER_HPP
