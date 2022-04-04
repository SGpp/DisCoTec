#ifndef TEST_HELPER_HPP
#define TEST_HELPER_HPP

// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include <complex>
// #include <boost/test/floating_point_comparison.hpp>
// new header for boost >= 1.59
#include <boost/test/tools/floating_point_comparison.hpp>

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

  static void testStrayMessages(MPI_Comm comm = MPI_COMM_WORLD){
    // general test for stray messages
    int flag;
    MPI_Status status;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &flag, &status);
    if (flag){
      int number_amount;
      MPI_Get_count(&status, MPI_CHAR, &number_amount);
      std::cout << getRank(MPI_COMM_WORLD) << " received " << number_amount << " from "
                << status.MPI_SOURCE << " with tag " << status.MPI_TAG << std::endl;
    }
    BOOST_CHECK(flag == false);
  }

  struct BarrierAtEnd {
    BarrierAtEnd() = default;
    ~BarrierAtEnd() {
      MPI_Barrier(MPI_COMM_WORLD);
      TestHelper::testStrayMessages();
    }
  };
}  // namespace TestHelper

#endif  // TEST_HELPER_HPP
