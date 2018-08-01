#ifndef TEST_HELPER_HPP
#define TEST_HELPER_HPP

#include <mpi.h>
#include <complex>

class TestHelper {
  static constexpr double tolerance = 1e-12;

 public:
  static bool checkNumProcs(int nprocs) {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size >= nprocs;
  }

  template <class T>
  static bool equals(const T& a, const T& b) {
    return std::abs(a - b) <= tolerance;
  }

  template <class T>
  static bool equals(const T& a, const T& b, double tol) {
    return std::abs(a - b) <= tol;
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

  static int getRank(MPI_Comm comm) {
    int rank;
    MPI_Comm_rank(comm, &rank);
    return rank;
  }
};

#endif  // TEST_HELPER_HPP
