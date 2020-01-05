#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <chrono>
#include <limits>
#include <random>
#include <thread>

#include "sgpp/distributedcombigrid/loadmodel/LearningLoadModel.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupSignals.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/mpi/MPIUtils.hpp"
#include "sgpp/distributedcombigrid/utils/Stats.hpp"
#include "test_helper.hpp"

using namespace combigrid;

void testDataSave(int size) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(size));

  CommunicatorType comm = TestHelper::getComm(size);
  if (comm == MPI_COMM_NULL) {
    return;
  }

  combigrid::Stats::initialize();

  size_t ngroup = size - 1;
  size_t nprocs = 1;
  theMPISystem()->initWorldReusable(comm, ngroup, nprocs);

  WORLD_MANAGER_EXCLUSIVE_SECTION {
    std::vector<LevelVector> lvv;
    for (long int i = 0; i < ngroup; ++i) {
      lvv.push_back({i});
    }
    auto loadModel = std::unique_ptr<LoadModel>(new LearningLoadModel(lvv));
    for (size_t j = 0; j < 600; ++j) {
      for (size_t i = 0; i < ngroup; ++i) {
        durationInformation recvbuf;
        MPI_Status stat;
        MPIUtils::receiveClass(&recvbuf, MPI_ANY_SOURCE, theMPISystem()->getGlobalComm());
        if (LearningLoadModel* llm = dynamic_cast<LearningLoadModel*>(loadModel.get())) {
          llm->addDataPoint(recvbuf, lvv.at(recvbuf.task_id)); 
        }
      }
    }
    // test loadmodel
    for (long int i = 0; i < ngroup; ++i) {
      // std::cout << "llm eval " << loadModel->eval({i}) << std::endl;
      BOOST_TEST(loadModel->eval({i}) == 1000000 * i);
    }
  }
  else {
    // the other processes act as the process managers in the combigrid setup
    uint nprocs = 64000;

    long unsigned int d = 1000000 * getCommRank(comm);

    Stats::Event e = Stats::Event();
    e.end = e.start + std::chrono::microseconds(d);
    durationInformation info = {TestHelper::getRank(comm), Stats::getEventDurationInUsec(e), 12.34, 0.00001, 1234, nprocs};
    // MPI_Request request;
    // send durationInfo to manager
    for (size_t i = 0; i < 600; ++i) { 
      MPIUtils::sendClass(&info, theMPISystem()->getManagerRank(), theMPISystem()->getGlobalComm());
      // TODO see if we can send asynchronously
    }
  }

  //TODO currently writing data only; need to implement eval properly
  // // see if reading the data works, too
  // WORLD_MANAGER_EXCLUSIVE_SECTION {
  //   std::cout << "reading data as if for next run" << std::endl;
  //   std::vector<LevelVector> lvv;
  //   for (long int i = 1; i < ngroup; ++i) {
  //     lvv.push_back({i});
  //   }
  //   auto loadModel = std::unique_ptr<LoadModel>(new LearningLoadModel(lvv));

  //   // test loadmodel
  //   for (long int i = 1; i < ngroup; ++i) {
  //     // std::cout << "llm eval " << loadModel->eval({i}) << std::endl;
  //     BOOST_TEST(loadModel->eval({i}) == 1000000 * i);
  //   }
  // }
  combigrid::Stats::finalize();
}

BOOST_AUTO_TEST_SUITE(loadmodel)

BOOST_AUTO_TEST_CASE(test_2) {
    testDataSave(2);
}

BOOST_AUTO_TEST_CASE(test_9, *boost::unit_test::timeout(120)) {
  testDataSave(9);
  MPI_Barrier(MPI_COMM_WORLD);
}

BOOST_AUTO_TEST_SUITE_END()
