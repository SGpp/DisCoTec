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
#include "sgpp/distributedcombigrid/utils/MonteCarlo.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "stdlib.h"
#include "test_helper.hpp"

using namespace combigrid;

BOOST_CLASS_EXPORT(TaskCount)

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
	// checking absolute and real value since comparing std::complex may be tricky
        BOOST_CHECK_CLOSE(std::abs(expected), std::abs(occuring), TestHelper::tolerance);
	BOOST_CHECK_CLOSE(std::real(expected), std::real(occuring), TestHelper::tolerance);
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
    params.setParallelization({static_cast<IndexType>(nprocs), 1});

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
    std::cout << "combined " << ngroup << " " << nprocs << std::endl;

    std::string filename("integration_" + std::to_string(ncombi) + ".raw");
    BOOST_TEST_CHECKPOINT("write solution " + filename);
    Stats::startEvent("manager write solution");
    manager.parallelEval( lmax, filename, 0 );
    Stats::stopEvent("manager write solution");
    std::cout << "wrote solution  " << ngroup << " " << nprocs << std::endl;

    // test Monte-Carlo interpolation
    BOOST_TEST_CHECKPOINT("MC interpolation coordinates");
    auto interpolationCoords = montecarlo::getRandomCoordinates(1000, dim);
    BOOST_TEST_CHECKPOINT("MC interpolation");
    Stats::startEvent("manager interpolate");
    auto values = manager.interpolateValues(interpolationCoords);
    Stats::stopEvent("manager interpolate");
    std::cout << "did interpolation " << ngroup << " " << nprocs << std::endl;

    TestFnCount<CombiDataType> initialFunction;
    for (size_t i = 0; i < interpolationCoords.size(); ++i) {
      if (std::abs(initialFunction(interpolationCoords[i], ncombi) - values[i]) > TestHelper::tolerance) {
        std::cout << "err " << interpolationCoords.size() <<interpolationCoords[i] << " " << i << std::endl;
      }
      auto ref = initialFunction(interpolationCoords[i], ncombi);
      BOOST_CHECK_CLOSE(std::abs(ref), std::abs(values[i]), TestHelper::tolerance);
      BOOST_CHECK_CLOSE(std::real(ref), std::real(values[i]), TestHelper::tolerance);
    }

    manager.exit();

    // if output files are not needed, remove them right away
    remove(("integration_" + std::to_string(ncombi) + "_0.raw").c_str());
    remove(("integration_" + std::to_string(ncombi) + "_0.raw_header").c_str());

    TestHelper::testStrayMessages(theMPISystem()->getGlobalComm());
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
      BOOST_TEST_CHECKPOINT("Last Successful Worker Signal " + std::to_string(signal));
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
    TestHelper::testStrayMessages(theMPISystem()->getLocalComm());
    MASTER_EXCLUSIVE_SECTION{
      TestHelper::testStrayMessages(theMPISystem()->getGlobalComm());
    }
  }

  combigrid::Stats::finalize();
  Stats::write("integration_" + std::to_string(ngroup) + "_" + std::to_string(nprocs) + ".json");

  MPI_Barrier(comm);
  TestHelper::testStrayMessages(comm);
}

/**
 * @brief Test for integrated passing of the hierarchical basis type
 *        (needs a lot of boilerplate code to set up manager etc, but its really only the
 *        `setCombiParametersHierarchicalBasesUniform<T>(params);` and
 *        `BOOST_TEST(dynamic_cast<T*>(b) != nullptr)` parts that are interesting)
 * @tparam T the hierarchical basis class to test (needs to be derived from BasisFunctionBasis,
 * serializable, exported (cf. BoostExports.hpp))
 */
template <typename T>
void checkPassingHierarchicalBases(size_t ngroup = 1, size_t nprocs = 1) {
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

  WORLD_MANAGER_EXCLUSIVE_SECTION {
    ProcessGroupManagerContainer pgroups;
    for (int i = 0; i < ngroup; ++i) {
      int pgroupRootID(i);
      pgroups.emplace_back(std::make_shared<ProcessGroupManager>(pgroupRootID));
    }

    auto loadmodel = std::unique_ptr<LoadModel>(new LinearLoadModel());
    std::vector<bool> boundary(dim, true);
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
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs, 2);
    params.setParallelization({static_cast<IndexType>(nprocs), 1});
    setCombiParametersHierarchicalBasesUniform<T>(params);

    // create abstraction for Manager
    ProcessManager manager{pgroups, tasks, params, std::move(loadmodel)};

    // the combiparameters are sent to all process groups before the
    // computations start
    manager.updateCombiParameters();

    manager.runfirst();
    manager.combine();
    manager.runnext();
    manager.combine();

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
    const auto& bases = pgroup.getCombiParameters().getHierarchicalBases();
    for (const auto& b : bases) {
      BOOST_TEST(dynamic_cast<T*>(b) != nullptr);
    }
    TestHelper::testStrayMessages(theMPISystem()->getLocalComm());
    MASTER_EXCLUSIVE_SECTION { TestHelper::testStrayMessages(theMPISystem()->getGlobalComm()); }
  }
  combigrid::Stats::finalize();
  MPI_Barrier(comm);
  TestHelper::testStrayMessages(comm);
}

#ifndef ISGENE // integration tests won't work with ISGENE because of worker magic
BOOST_FIXTURE_TEST_SUITE(integration, TestHelper::BarrierAtEnd, *boost::unit_test::timeout(120))

BOOST_AUTO_TEST_CASE(test_1, *boost::unit_test::tolerance(TestHelper::higherTolerance) ) {
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
}

BOOST_AUTO_TEST_SUITE_END()
#endif

BOOST_AUTO_TEST_CASE(test_2) { checkPassingHierarchicalBases<HierarchicalHatBasisFunction>(1, 1); }

BOOST_AUTO_TEST_CASE(test_3) { checkPassingHierarchicalBases<FullWeightingBasisFunction>(1, 2); }

BOOST_AUTO_TEST_CASE(test_4) {
  checkPassingHierarchicalBases<FullWeightingPeriodicBasisFunction>(2, 2);
}

BOOST_AUTO_TEST_CASE(test_5) { checkPassingHierarchicalBases<BiorthogonalBasisFunction>(1, 4); }

BOOST_AUTO_TEST_CASE(test_6) {
  checkPassingHierarchicalBases<BiorthogonalPeriodicBasisFunction>(4, 2);
}
