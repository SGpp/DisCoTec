#define BOOST_TEST_DYN_LINK
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <boost/serialization/export.hpp>
#include <boost/test/unit_test.hpp>
#include <cstdio>

#include "TaskCount.hpp"
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
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "stdlib.h"
#include "test_helper.hpp"

using namespace combigrid;

// omitted here, because already defined in test_thirdLevel.cpp
// BOOST_CLASS_EXPORT(TaskCount)

bool checkReducedFullGridIntegration(ProcessGroupWorker& worker, int nrun) {
  auto tasks = worker.getTasks();
  int numGrids = (int)worker.getCombiParameters().getNumGrids();

  BOOST_CHECK(tasks.size() > 0);
  BOOST_CHECK(numGrids > 0);

  // to check if any data was actually compared
  bool any = false;

  for (Task* t : tasks) {
    for (int g = 0; g < numGrids; g++) {
      DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid(g);
      for (auto b : dfg.returnBoundaryFlags()) {
        BOOST_CHECK(b);
      }

      TestFnCount<CombiDataType> initialFunction;
      for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
        std::vector<double> coords(dfg.getDimension());
        dfg.getCoordsLocal(li, coords);
        CombiDataType expected = initialFunction(coords, nrun);
        CombiDataType occuring = dfg.getData()[li];
        for (auto& c : coords){
          BOOST_CHECK(c >= 0.);
          BOOST_CHECK(c <= 1.);
        }
        BOOST_CHECK_CLOSE(expected, occuring, TestHelper::tolerance);
        // BOOST_REQUIRE_CLOSE(expected, occuring, TestHelper::tolerance);
        any = true;
      }
    }
  }
  BOOST_CHECK(any);
  return any;
}

void checkIntegration(size_t ngroup = 1, size_t nprocs = 1, bool boundaryV = true) {
  size_t size = ngroup * nprocs + 1;
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(size));

  CommunicatorType comm = TestHelper::getComm(size);
  if (comm == MPI_COMM_NULL) {
    BOOST_TEST_CHECKPOINT("drop out of test comm");
    return;
  }

  combigrid::Stats::initialize();

  theMPISystem()->initWorldReusable(comm, ngroup, nprocs);
  // theMPISystem()->init(ngroup, nprocs);

  DimType dim = 2;
  LevelVector lmin(dim, 2);
  LevelVector lmax(dim, 5);

  size_t ncombi = 4;

  BOOST_CHECK_EQUAL(theMPISystem()->getWorldSize(), size);
  BOOST_CHECK_EQUAL(getCommSize(theMPISystem()->getWorldComm()), size);

  WORLD_MANAGER_EXCLUSIVE_SECTION {
    // make sure the manager's ranks are set right
    BOOST_CHECK_EQUAL(getCommRank(theMPISystem()->getWorldComm()), size - 1);
    BOOST_CHECK_EQUAL(getCommSize(theMPISystem()->getGlobalComm()), ngroup + 1);
    BOOST_CHECK_EQUAL(getCommRank(theMPISystem()->getGlobalComm()), ngroup);

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
    TaskContainer tasks;
    std::vector<int> taskIDs;
    for (size_t i = 0; i < levels.size(); i++) {
      Task* t = new TaskCount(dim, levels[i], boundary, coeffs[i], loadmodel.get());
      tasks.push_back(t);
      taskIDs.push_back(t->getID());
    }

    // create combiparameters
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs, ncombi);
    params.setParallelization({nprocs, 1});

    // create abstraction for Manager
    ProcessManager manager{pgroups, tasks, params, std::move(loadmodel)};

    // the combiparameters are sent to all process groups before the
    // computations start
    manager.updateCombiParameters();

    /* distribute task according to load model and start computation for
     * the first time */
    BOOST_TEST_CHECKPOINT("run first");
    manager.runfirst();

    for (size_t it = 0; it < ncombi - 1; ++it) {
      BOOST_TEST_CHECKPOINT("combine");
      manager.combine();

      BOOST_TEST_CHECKPOINT("run next");
      manager.runnext();
    }
    manager.combine();

    std::string filename("integration_" + std::to_string(ncombi) + ".raw");
    Stats::startEvent("manager write solution");
    manager.parallelEval( lmax, filename, 0 );
    Stats::stopEvent("manager write solution");

    manager.exit();

    // if output files are not needed, remove them right away
    remove(("integration_" + std::to_string(ncombi) + "_0.raw").c_str());
    remove(("integration_" + std::to_string(ncombi) + "_0.raw_header").c_str());
  }
  else {
    BOOST_CHECK_EQUAL(getCommSize(theMPISystem()->getLocalComm()), nprocs);
    if (nprocs == 1){
      BOOST_CHECK(theMPISystem()->isMaster());
    }
    BOOST_TEST_CHECKPOINT("Worker starts");
    ProcessGroupWorker pgroup;
    SignalType signal = -1;
    // omitting to count RUN_FIRST signal, as it is executed once for every task
    int nrun = 1;
    while (signal != EXIT) {
      BOOST_TEST_CHECKPOINT(signal);
      signal = pgroup.wait();

      if (signal == RUN_NEXT) {
        BOOST_CHECK(pgroup.getTasks().size() > 0);
        ++nrun;
      }
      if (signal == COMBINE) {
        // after combination check workers' grids
        BOOST_CHECK(checkReducedFullGridIntegration(pgroup, nrun));
      }
    }
    BOOST_CHECK_EQUAL(nrun, ncombi);
  }

  combigrid::Stats::finalize();
  Stats::write("integration_" + std::to_string(ngroup) + "_" + std::to_string(nprocs) + ".json");

  MPI_Barrier(comm);
  TestHelper::testStrayMessages(comm);
}

BOOST_AUTO_TEST_SUITE(integration)

BOOST_AUTO_TEST_CASE(test_1, *boost::unit_test::tolerance(TestHelper::higherTolerance) *
                                 boost::unit_test::timeout(60)) {
  for (bool boundary : {true}) {
    for (size_t ngroup : {1, 2, 3, 4}) {
      for (size_t nprocs : {1, 2}) {
        std::cout << "integration/test_1 " << ngroup << " " << nprocs << std::endl;
        checkIntegration(ngroup, nprocs, boundary);
        MPI_Barrier(MPI_COMM_WORLD);
      }
    }
    for (size_t ngroup : {1, 2}) {
      for (size_t nprocs : {4}) { //TODO currently fails for non-power-of-2-decompositions
        std::cout << "integration/test_1 " << ngroup << " " << nprocs << std::endl;
        checkIntegration(ngroup, nprocs, boundary);
        MPI_Barrier(MPI_COMM_WORLD);
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  TestHelper::testStrayMessages();
}

BOOST_AUTO_TEST_SUITE_END()
