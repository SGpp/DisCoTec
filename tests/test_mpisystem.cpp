#define BOOST_TEST_DYN_LINK
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <boost/serialization/export.hpp>
#include <boost/test/unit_test.hpp>
#include <cstdio>

#include "mpi/MPIMemory.hpp"
#include "mpi/MPISystem.hpp"
#include "test_helper.hpp"
#include "utils/Stats.hpp"

using namespace combigrid;

unsigned long checkMPIMemory(size_t ngroup, size_t nprocs) {
  size_t size = ngroup * nprocs + 1;
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(size));

  CommunicatorType comm = TestHelper::getComm(size);
  if (comm == MPI_COMM_NULL) {
    BOOST_TEST_CHECKPOINT("drop out of test comm");
    return 0;
  }

  combigrid::Stats::initialize();
  theMPISystem()->initWorldReusable(comm, ngroup, nprocs);
  // theMPISystem()->init(ngroup, nprocs);

  BOOST_CHECK_EQUAL(theMPISystem()->getWorldSize(), size);
  BOOST_CHECK_EQUAL(getCommSize(theMPISystem()->getWorldComm()), size);

  unsigned long local_vmrss, local_vmsize = 0;
  WORLD_MANAGER_EXCLUSIVE_SECTION {}
  else {
    BOOST_TEST_CHECKPOINT("after MPI init ");
    // mpimemory::print_memory_usage_local();
    mpimemory::get_memory_usage_local_kb(&local_vmrss, &local_vmsize);
  }

  combigrid::Stats::finalize();
  MPI_Barrier(comm);
  BOOST_CHECK(!TestHelper::testStrayMessages(comm));
  return local_vmsize;
}

BOOST_FIXTURE_TEST_SUITE(mpisystem, TestHelper::BarrierAtEnd, *boost::unit_test::timeout(60))

BOOST_AUTO_TEST_CASE(test_1, *boost::unit_test::timeout(60)) {
  std::vector<size_t> groupSizes{1, 2, 4, 8};
  std::vector<long unsigned int> vmSizes(groupSizes.size());
  for (size_t i = 0; i < groupSizes.size(); ++i) {
    vmSizes[i] = checkMPIMemory(1, groupSizes[i]);
    if (i > 0) {
      // compare allocated memory sizes
      // check for linear scaling (grace 20%)
      if (TestHelper::getRank(MPI_COMM_WORLD) == 0) {
        BOOST_TEST(static_cast<double>(vmSizes[i]) <=
                   (vmSizes[0] * groupSizes[i] / groupSizes[0] * 1.2));
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  BOOST_CHECK(!TestHelper::testStrayMessages());
}

BOOST_AUTO_TEST_CASE(test_with_manager) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(6));
  for (const auto& procs : std::vector<std::vector<int>>{{2, 3}, {3, 2}, {1, 6}, {6, 1}}) {
    CommunicatorType comm = TestHelper::getComm(procs[0] * procs[1] + 1);
    if (comm != MPI_COMM_NULL) {
      Stats::initialize();
      theMPISystem()->initWorldReusable(comm, procs[0], procs[1]);
      BOOST_CHECK_EQUAL(theMPISystem()->getWorldSize(), procs[0] * procs[1] + 1);
      BOOST_CHECK_EQUAL(theMPISystem()->getWorldRank(), getCommRank(comm));

      WORLD_MANAGER_EXCLUSIVE_SECTION {
        BOOST_CHECK_EQUAL(theMPISystem()->getWorldRank(), procs[0] * procs[1]);
        BOOST_CHECK_EQUAL(theMPISystem()->getLocalComm(), MPI_COMM_NULL);
        BOOST_CHECK_LT(theMPISystem()->getLocalRank(), 0);
        BOOST_CHECK_EQUAL(theMPISystem()->getGlobalReduceComm(), MPI_COMM_NULL);
        BOOST_CHECK_LT(theMPISystem()->getGlobalReduceRank(), 0);
        BOOST_CHECK_NE(theMPISystem()->getGlobalComm(), MPI_COMM_NULL);
        BOOST_CHECK_EQUAL(theMPISystem()->getGlobalRank(), procs[0]);
      }
      else {
        BOOST_CHECK_NE(theMPISystem()->getLocalComm(), MPI_COMM_NULL);
        BOOST_CHECK_LT(theMPISystem()->getLocalRank(), procs[1]);
        BOOST_CHECK_NE(theMPISystem()->getGlobalReduceComm(), MPI_COMM_NULL);
        BOOST_CHECK_LT(theMPISystem()->getGlobalReduceRank(), procs[0]);
        MASTER_EXCLUSIVE_SECTION {
          BOOST_CHECK_NE(theMPISystem()->getGlobalComm(), MPI_COMM_NULL);
          BOOST_CHECK_LT(theMPISystem()->getGlobalRank(), procs[0]);
        }
        else {
          BOOST_CHECK_EQUAL(theMPISystem()->getGlobalComm(), MPI_COMM_NULL);
          BOOST_CHECK_LT(theMPISystem()->getGlobalRank(), 0);
        }
      }
      Stats::finalize();
    }
  }
}

BOOST_AUTO_TEST_CASE(test_manager_only) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(1));
  auto procs = std::vector<int>{0, 0};
  CommunicatorType comm = TestHelper::getComm(1);
  if (comm != MPI_COMM_NULL) {
    Stats::initialize();
    theMPISystem()->initWorldReusable(comm, procs[0], procs[1]);
    BOOST_CHECK_EQUAL(theMPISystem()->getWorldSize(), procs[0] * procs[1] + 1);
    BOOST_CHECK_EQUAL(theMPISystem()->getWorldRank(), getCommRank(comm));

    WORLD_MANAGER_EXCLUSIVE_SECTION {
      BOOST_CHECK_EQUAL(theMPISystem()->getWorldRank(), procs[0] * procs[1]);
      BOOST_CHECK_EQUAL(theMPISystem()->getLocalComm(), MPI_COMM_NULL);
      BOOST_CHECK_LT(theMPISystem()->getLocalRank(), 0);
      BOOST_CHECK_EQUAL(theMPISystem()->getGlobalReduceComm(), MPI_COMM_NULL);
      BOOST_CHECK_LT(theMPISystem()->getGlobalReduceRank(), 0);
      // BOOST_CHECK_NE(theMPISystem()->getGlobalComm(), MPI_COMM_NULL);
      // BOOST_CHECK_EQUAL(theMPISystem()->getGlobalRank(), procs[0]);
    }
    else {
      BOOST_CHECK(false);
    }
    Stats::finalize();
  }
}


BOOST_AUTO_TEST_SUITE_END()
