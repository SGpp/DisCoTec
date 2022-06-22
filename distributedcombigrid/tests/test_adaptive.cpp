#define BOOST_TEST_DYN_LINK
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <boost/serialization/export.hpp>
#include <boost/test/unit_test.hpp>
#include <cstdio>

#include "TaskInterpolateWave.hpp"
#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LearningLoadModel.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LinearLoadModel.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupManager.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupWorker.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessManager.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/DistributedSparseGridUniform.hpp"
#include "sgpp/distributedcombigrid/task/Task.hpp"
#include "sgpp/distributedcombigrid/utils/Config.hpp"
#include "sgpp/distributedcombigrid/utils/MonteCarlo.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "stdlib.h"
#include "test_helper.hpp"

using namespace combigrid;

BOOST_CLASS_EXPORT(TaskInterpolateWave)

void checkAdaptivity(size_t ngroup = 1, size_t nprocs = 1, bool boundaryV = true) {
  size_t size = ngroup * nprocs + 1;
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(size));

  CommunicatorType comm = TestHelper::getComm(size);
  if (comm == MPI_COMM_NULL) {
    BOOST_TEST_CHECKPOINT("drop out of test comm");
    return;
  }

  combigrid::Stats::initialize();

  theMPISystem()->initWorldReusable(comm, ngroup, nprocs);

  DimType dim = 2;
  LevelVector lmin(dim, 2);
  LevelVector lmax(dim, 5);

  size_t ncombi = 10;

  WORLD_MANAGER_EXCLUSIVE_SECTION {
    ProcessGroupManagerContainer pgroups;
    for (int i = 0; i < ngroup; ++i) {
      int pgroupRootID(i);
      pgroups.emplace_back(std::make_shared<ProcessGroupManager>(pgroupRootID));
    }

    auto loadmodel = std::unique_ptr<LoadModel>(new LinearLoadModel());

    std::vector<bool> boundary(dim, boundaryV);

    CombiMinMaxScheme combischeme(dim, lmin, lmax);
    combischeme.createAdaptiveCombischeme();

    std::vector<LevelVector> levels = combischeme.getCombiSpaces();
    std::vector<combigrid::real> coeffs = combischeme.getCoeffs();

    // create Tasks
    double timeStepSize =
        1. / ncombi;  // try ncombi time steps -> end should be the same as initial state
    TaskContainer tasks;
    std::vector<size_t> taskIDs;
    for (size_t i = 0; i < levels.size(); i++) {
      Task* t = new TaskInterpolateWave(dim, levels[i], boundary, coeffs[i], loadmodel.get(),
                                        timeStepSize);
      tasks.push_back(t);
      taskIDs.push_back(t->getID());
    }

    // why are there duplicate task IDs over different calls to this function??
    // ah its because different ranks will be master, so there are different Task::count instances
    // std::cout << "taskIDs " << taskIDs << std::endl;

    // create combiparameters
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs, ncombi);
    params.setParallelization({static_cast<IndexType>(nprocs), 1});

    BOOST_TEST_CHECKPOINT("Manager starts");
    // create abstraction for Manager
    ProcessManager manager{pgroups, tasks, params, std::move(loadmodel)};

    // the combiparameters are sent to all process groups before the
    // computations start
    manager.updateCombiParameters();

    /* distribute task according to load model and start computation for
     * the first time */
    BOOST_TEST_CHECKPOINT("run first");
    manager.runfirst(true);

    for (size_t it = 0; it < ncombi - 1; ++it) {
      BOOST_TEST_CHECKPOINT("combine");
      manager.combine();

      /**
       * TODO adaptation
       *
       */

      std::string filename("adaptive_" + std::to_string(it) + ".raw");
      BOOST_TEST_CHECKPOINT("write solution " + filename);
      Stats::startEvent("manager write solution");
      // this writes .raw files, which can be plotted in 2-D e.g. with ../tools/raw_to_image.py
      manager.parallelEval(lmax, filename, 0);
      // this writes human-readable files, which can be plotted in 3-D with
      // ../tools/visualize_sg_minmax.py
      manager.writeSparseGridMinMaxCoefficients("adaptive_" + std::to_string(it) +
                                                "_sparse_minmax");
      Stats::stopEvent("manager write solution");

      BOOST_TEST_CHECKPOINT("run next");
      manager.runnext();
    }
    manager.combine();
    std::string filename("adaptive_" + std::to_string(ncombi) + ".raw");
    BOOST_TEST_CHECKPOINT("write solution " + filename);
    Stats::startEvent("manager write solution");
    manager.parallelEval(lmax, filename, 0);
    manager.writeSparseGridMinMaxCoefficients("adaptive_" + std::to_string(ncombi) +
                                              "_sparse_minmax");
    Stats::stopEvent("manager write solution");
    manager.exit();

    TestHelper::testStrayMessages(theMPISystem()->getGlobalComm());
  }
  else {
    BOOST_TEST_CHECKPOINT("Worker starts");
    ProcessGroupWorker pgroup;
    SignalType signal = -1;
    // omitting to count RUN_FIRST signal, as it is executed once for every task
    while (signal != EXIT) {
      BOOST_TEST_CHECKPOINT("Last Successful Worker Signal " + std::to_string(signal));
      signal = pgroup.wait();
    }
    TestHelper::testStrayMessages(theMPISystem()->getLocalComm());
    MASTER_EXCLUSIVE_SECTION { TestHelper::testStrayMessages(theMPISystem()->getGlobalComm()); }
  }

  combigrid::Stats::finalize();

  MPI_Barrier(comm);
  TestHelper::testStrayMessages(comm);
}

BOOST_FIXTURE_TEST_SUITE(adaptive, TestHelper::BarrierAtEnd, *boost::unit_test::timeout(90))

BOOST_AUTO_TEST_CASE(test_1) {
  for (bool boundary : {true}) {
    for (size_t ngroup : {2}) {
      for (size_t nprocs : {4}) {
        checkAdaptivity(ngroup, nprocs, boundary);
        MPI_Barrier(MPI_COMM_WORLD);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    TestHelper::testStrayMessages();
  }
}

BOOST_AUTO_TEST_SUITE_END()
