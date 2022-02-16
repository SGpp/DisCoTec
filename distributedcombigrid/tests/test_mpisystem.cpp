#define BOOST_TEST_DYN_LINK
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <boost/serialization/export.hpp>
#include <boost/test/unit_test.hpp>
#include <cstdio>

// #include "sgpp/distributedcombigrid/loadmodel/LearningLoadModel.hpp"
// #include "sgpp/distributedcombigrid/loadmodel/LinearLoadModel.hpp"
// #include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
// #include "sgpp/distributedcombigrid/manager/ProcessGroupManager.hpp"
// #include "sgpp/distributedcombigrid/manager/ProcessGroupWorker.hpp"
// #include "sgpp/distributedcombigrid/manager/ProcessManager.hpp"
#include "sgpp/distributedcombigrid/mpi/MPIMemory.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
// #include "sgpp/distributedcombigrid/task/Task.hpp"
// #include "sgpp/distributedcombigrid/utils/Config.hpp"
#include "sgpp/distributedcombigrid/utils/Stats.hpp"
#include "test_helper.hpp"

using namespace combigrid;

void checkMPIRanksAndCommunication(size_t ngroup, size_t nprocs) {
  // size_t size = ngroup * nprocs + 1;
  // BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(size));

  // CommunicatorType comm = TestHelper::getComm(size);
  // if (comm == MPI_COMM_NULL) {
  //   BOOST_TEST_CHECKPOINT("drop out of test comm");
  //   return;
  // }

  // theMPISystem()->initWorldReusable(comm, ngroup, nprocs);
  // theMPISystem()->init(ngroup, nprocs);

  // WORLD_MANAGER_EXCLUSIVE_SECTION {
  //   // make sure the manager's ranks are set right
  //   BOOST_CHECK_EQUAL(getCommRank(theMPISystem()->getWorldComm()), size - 1);
  //   BOOST_CHECK_EQUAL(getCommSize(theMPISystem()->getGlobalComm()), ngroup + 1);
  //   BOOST_CHECK_EQUAL(getCommRank(theMPISystem()->getGlobalComm()), ngroup);

  //   DimType dim = 6;
  //   LevelVector lmin(dim, 2);
  //   LevelVector lmax(dim, 3);

  //   ProcessGroupManagerContainer pgroups;
  //   for (int i = 0; i < ngroup; ++i) {
  //     int pgroupRootID(i);
  //     pgroups.emplace_back(std::make_shared<ProcessGroupManager>(pgroupRootID));
  //   }

  //   auto loadmodel = std::unique_ptr<LoadModel>(new LinearLoadModel());

  //   std::vector<bool> boundary(dim, boundaryV);

  //   CombiMinMaxScheme combischeme(dim, lmin, lmax);
  //   combischeme.createAdaptiveCombischeme();

  //   std::vector<LevelVector> levels = combischeme.getCombiSpaces();
  //   std::vector<combigrid::real> coeffs = combischeme.getCoeffs();

  //   // create Tasks
  //   TaskContainer tasks;
  //   std::vector<size_t> taskIDs;
  //   for (size_t i = 0; i < levels.size(); i++) {
  //     Task* t = new TaskCount(dim, levels[i], boundary, coeffs[i], loadmodel.get());
  //     tasks.push_back(t);
  //     taskIDs.push_back(t->getID());
  //   }

  //   // create combiparameters
  //   CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs, ncombi);
  //   params.setParallelization({static_cast<IndexType>(nprocs), 1});

  //   // create abstraction for Manager
  //   ProcessManager manager{pgroups, tasks, params, std::move(loadmodel)};

  //   // the combiparameters are sent to all process groups before the
  //   // computations start
  //   manager.updateCombiParameters();

  //   manager.exit();
  //   TestHelper::testStrayMessages(theMPISystem()->getGlobalComm());
  // }
  // else {
  //   BOOST_CHECK_EQUAL(getCommSize(theMPISystem()->getLocalComm()), nprocs);
  //   if (nprocs == 1) {
  //     BOOST_CHECK(theMPISystem()->isMaster());
  //   }
  //   BOOST_TEST_CHECKPOINT("Worker starts");
  //   ProcessGroupWorker pgroup;
  //   MASTER_EXCLUSIVE_SECTION { std::cout << "at worker construction " << std::flush; }
  //   mpimemory::print_memory_usage_local();
  //   SignalType signal = -1;
  //   while (signal != EXIT) {
  //     BOOST_TEST_CHECKPOINT(signal);
  //     MASTER_EXCLUSIVE_SECTION { std::cout << "at signal " << signal << std::flush; }
  //     mpimemory::print_memory_usage_local();
  //     signal = pgroup.wait();
  //   }
  // }
  // MPI_Barrier(comm);
  // TestHelper::testStrayMessages(comm);
}

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
  TestHelper::testStrayMessages(comm);
  return local_vmsize;
}

BOOST_AUTO_TEST_SUITE(mpisystem)

BOOST_AUTO_TEST_CASE(test_1, *boost::unit_test::timeout(60)) {
  std::vector<size_t> groupSizes {1,2,4,8};
  std::vector<long unsigned int> vmSizes(groupSizes.size());
  for (size_t i = 0; i < groupSizes.size(); ++i) {
    vmSizes[i] = checkMPIMemory(1, groupSizes[i]);
    if (i > 0) {
      // compare allocated memory sizes
      // check for linear scaling (grace 20%)
      if (TestHelper::getRank(MPI_COMM_WORLD) == 0) {
        BOOST_TEST(static_cast<double>(vmSizes[i]) <= (vmSizes[0] * groupSizes[i] / groupSizes[0] * 1.2));
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  TestHelper::testStrayMessages();
}

//TODO make test for checkMPIRanksAndCommunication

BOOST_AUTO_TEST_SUITE_END()
