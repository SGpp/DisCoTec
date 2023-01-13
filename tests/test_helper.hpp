#ifndef TEST_HELPER_HPP
#define TEST_HELPER_HPP

// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include <algorithm>
#include <complex>
#include <numeric>
#include <vector>
// new header for boost >= 1.59
#include <boost/test/tools/floating_point_comparison.hpp>
#include "utils/Stats.hpp"

namespace TestHelper{
  static constexpr double tolerance = 1e-12;
  static constexpr double higherTolerance = 1e-5;

  static inline bool checkNumMPIProcsAvailable(int nprocs) {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size >= nprocs;
  }

  static inline bool checkNumMPIProcsAvailable(size_t nprocs) {
    return checkNumMPIProcsAvailable(static_cast<int>(nprocs));
  }

  static inline MPI_Comm getComm(int nprocs) {
    BOOST_CHECK(TestHelper::checkNumMPIProcsAvailable(nprocs));
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

  /**
   * @brief Get a cartesian communicator of specified extents
   *
   * @param procs a vector of the extents
   * @return MPI_Comm the cartesian communicator (or MPI_COMM_NULL)
   */
  static inline MPI_Comm getComm(std::vector<int> procs, std::vector<int> periods = {}) {
    auto comm = getComm(std::accumulate(procs.begin(), procs.end(), 1, std::multiplies<int>()));
    if (comm == MPI_COMM_NULL) {
      return comm;
    } else {
      if (procs.size() != periods.size()) {
        // Make all dimensions not periodic
        periods.resize(procs.size(), 0);
      }
      // let MPI assign arbitrary ranks?
      int reorder = false;
      // Create a communicator given the topology
      MPI_Comm new_communicator;
      MPI_Cart_create(comm, static_cast<int>(procs.size()), procs.data(), periods.data(), reorder,
                      &new_communicator);
      return new_communicator;
    }
  }

  static inline int getRank(MPI_Comm comm) {
    int rank;
    MPI_Comm_rank(comm, &rank);
    return rank;
  }

  static bool testStrayMessages(MPI_Comm comm = MPI_COMM_WORLD) {
    // general test for stray messages
    int flag;
    MPI_Status status;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &flag, &status);
    int number_amount;
    if (flag) {
      MPI_Get_count(&status, MPI_CHAR, &number_amount);
      std::cout << getRank(MPI_COMM_WORLD) << " received " << number_amount << " bytes from "
                << status.MPI_SOURCE << " with tag " << status.MPI_TAG << std::endl;
    }
    // BOOST_CHECK(flag == false);
    if (flag) {
      std::vector<char> buffer(number_amount);
      MPI_Recv(buffer.data(), number_amount, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, comm,
               MPI_STATUS_IGNORE);
      std::cout << " content " ;
      for (const auto& c: buffer)
        std::cout << std::to_string(c) << " ";
      std::cout << std::endl;
    }
    return flag;
  }

  struct BarrierAtEnd {
    BarrierAtEnd() = default;
    ~BarrierAtEnd() {
      BOOST_CHECK(!combigrid::Stats::isInitialized());
      MPI_Barrier(MPI_COMM_WORLD);
      BOOST_CHECK(!TestHelper::testStrayMessages());
    }
  };
}  // namespace TestHelper

#endif  // TEST_HELPER_HPP
