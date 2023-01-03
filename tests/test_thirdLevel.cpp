#define BOOST_TEST_DYN_LINK
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <boost/serialization/export.hpp>
#include <boost/test/unit_test.hpp>
#include <cstdio>

#include "TaskConstParaboloid.hpp"
#include "TaskCount.hpp"
#include "combischeme/CombiMinMaxScheme.hpp"
#include "combischeme/CombiThirdLevelScheme.hpp"
#include "loadmodel/LearningLoadModel.hpp"
#include "loadmodel/LinearLoadModel.hpp"
#include "manager/CombiParameters.hpp"
#include "manager/ProcessGroupManager.hpp"
#include "manager/ProcessGroupWorker.hpp"
#include "manager/ProcessManager.hpp"
#include "sparsegrid/DistributedSparseGridUniform.hpp"
#include "task/Task.hpp"
#include "utils/Config.hpp"
#include "utils/MonteCarlo.hpp"
#include "utils/Types.hpp"
#include "stdlib.h"
#include "test_helper.hpp"

using namespace combigrid;

BOOST_CLASS_EXPORT(TaskConstParaboloid)
BOOST_CLASS_EXPORT(TaskCount)

class TestParams {
 public:
  DimType dim = 2;
  LevelVector lmin;
  LevelVector lmax;
  BoundaryType boundary;
  unsigned int ngroup = 1;
  unsigned int nprocs = 1;
  unsigned int ncombi = 1;
  unsigned int sysNum = 0;
  const CommunicatorType& comm;
  std::string host = "localhost";
  unsigned short port = 9999;

  TestParams(DimType dim, LevelVector& lmin, LevelVector& lmax, BoundaryType boundary, unsigned int ngroup,
             unsigned int nprocs, unsigned int ncombi, unsigned int sysNum,
             const CommunicatorType& comm, const std::string& host = "localhost",
#ifdef NDEBUG
             unsigned short dataPort = 9999)
#else
             unsigned short dataPort = 7777)
#endif  // NDEBUG
      : dim(dim),
        lmin(lmin),
        lmax(lmax),
        boundary(boundary),
        ngroup(ngroup),
        nprocs(nprocs),
        ncombi(ncombi),
        sysNum(sysNum),
        comm(comm),
        host(host),
        port(dataPort) {
  }
};

/**
 * Checks if combination was successful.
 * Since the tasks don't evolve over time the expected result should match the
 * initial function values.
 */
bool checkReducedFullGrid(ProcessGroupWorker& worker, int nrun) {
  TaskContainer& tasks = worker.getTasks();
  int numGrids = (int)worker.getCombiParameters().getNumGrids();

  BOOST_CHECK(tasks.size() > 0);
  BOOST_CHECK(numGrids > 0);

  // to check if any data was actually compared
  bool any = false;

  for (Task* t : tasks) {
    for (int g = 0; g < numGrids; g++) {
      DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid(g);
      // dfg.print(std::cout);
      // std::cout << std::endl;
      // TestFnCount<CombiDataType> initialFunction;
      ParaboloidFn<CombiDataType> initialFunction;
      for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
        std::vector<double> coords(dfg.getDimension());
        dfg.getCoordsLocal(li, coords);
        // CombiDataType expected = initialFunction(coords, nrun);
        CombiDataType expected = initialFunction(coords);
        CombiDataType occuring = dfg.getData()[li];
        if (expected == 0.) {
          BOOST_CHECK_SMALL(occuring, 1e-300);
        } else {
          BOOST_REQUIRE_CLOSE(occuring, expected, TestHelper::tolerance);
        }
        any = true;
      }
    }
  }
  BOOST_CHECK(any);
  return any;
}

void assignProcsToSystems(unsigned int ngroup, unsigned int nprocs, unsigned int numSystems,
                          unsigned int& sysNum, CommunicatorType& newcomm) {


  int procsPerSys = ngroup * nprocs + 1;
  int totalProcs = numSystems * procsPerSys;

  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(totalProcs));
  auto smallerComm = TestHelper::getComm(totalProcs);
  if (smallerComm == MPI_COMM_NULL) {
    newcomm = MPI_COMM_NULL;
    sysNum = -1;
    return;
  } else {
    int size = 0;
    MPI_Comm_size(smallerComm, &size);
    BOOST_CHECK_EQUAL(size, totalProcs);
    auto rank = TestHelper::getRank(smallerComm);
    BOOST_CHECK_LT(rank, totalProcs);
    // assign procs to systems
    sysNum = unsigned(rank) / procsPerSys;
    int color = int(sysNum);
    BOOST_CHECK_LT(color, numSystems);
    int key = rank % (int)procsPerSys;
    BOOST_CHECK_LT(key, procsPerSys);

    MPI_Comm_split(smallerComm, color, key, &newcomm);
    MPI_Comm_size(newcomm, &size);
    BOOST_CHECK_EQUAL(size, procsPerSys);
  }
}

/** Runs the third level manager in the background as a forked child process */
void runThirdLevelManager(unsigned short port) {
  std::cout << "starting thirdLevelManager..." << std::endl;
  std::string command = "../third_level_manager/thirdLevelManager --port=" +
                        std::to_string(port) + " &";
  auto status = system(command.c_str());
  BOOST_WARN_GE(status, 0);
}

/** Runs the tl manager*/
#ifdef NDEBUG
void startInfrastructure(unsigned short port = 9999) {
#else
void startInfrastructure(unsigned short port = 7777) {
#endif // NDEBUG
  // give former infrastructure some time to shut down
  sleep(3);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    BOOST_TEST_CHECKPOINT("starting broker");
    runThirdLevelManager(port);
  }
  // give infrastructure some time to set up
  sleep(10);
}

void testCombineThirdLevel(TestParams& testParams, bool thirdLevelExtraSparseGrid = false) {
  BOOST_CHECK(testParams.comm != MPI_COMM_NULL);

  combigrid::Stats::initialize();

  theMPISystem()->initWorldReusable(testParams.comm, testParams.ngroup, testParams.nprocs);

  WORLD_MANAGER_EXCLUSIVE_SECTION {
    ProcessGroupManagerContainer pgroups;
    for (size_t i = 0; i < testParams.ngroup; ++i) {
      int pgroupRootID((int)i);
      pgroups.emplace_back(std::make_shared<ProcessGroupManager>(pgroupRootID));
    }

    auto loadmodel = std::unique_ptr<LoadModel>(new LinearLoadModel());
    std::vector<BoundaryType> boundary(testParams.dim, testParams.boundary);

    // create third level specific scheme
    CombiMinMaxScheme combischeme(testParams.dim, testParams.lmin, testParams.lmax);
    combischeme.createClassicalCombischeme();
    // combischeme.createAdaptiveCombischeme();
    // get full scheme first
    std::vector<LevelVector> fullLevels = combischeme.getCombiSpaces();
    std::vector<combigrid::real> fullCoeffs = combischeme.getCoeffs();
    // split scheme and assign each half to a system
    std::vector<LevelVector> levels;
    std::vector<combigrid::real> coeffs;
    CombiThirdLevelScheme::createThirdLevelScheme(fullLevels, fullCoeffs, testParams.sysNum, 2,
                                                  levels, coeffs);

    BOOST_REQUIRE_EQUAL(levels.size(), coeffs.size());
    // std::cout << "Combischeme " << testParams.sysNum << ":" << std::endl;
    // for (size_t i =0; i < levels.size(); ++i)
    //  std::cout << toString(levels[i]) << " " << coeffs[i]<< std::endl;

    // create Tasks
    TaskContainer tasks;
    std::vector<size_t> taskIDs;
    for (size_t i = 0; i < levels.size(); i++) {
      Task* t = new TaskConstParaboloid(levels[i], boundary, coeffs[i], loadmodel.get());
      // Task* t = new TaskCount(2, levels[i], boundary, coeffs[i], loadmodel.get());

      tasks.push_back(t);
      taskIDs.push_back(t->getID());
    }

    // create combiparameters
    std::vector<int> parallelization(testParams.dim, 1);
    parallelization[1] = static_cast<int>(testParams.nprocs);
    CombiParameters combiParams(testParams.dim, testParams.lmin, testParams.lmax, boundary, levels,
                                coeffs, taskIDs, testParams.ncombi, 1, parallelization,
                                LevelVector(testParams.dim, 0), LevelVector(testParams.dim, 1),
                                false, testParams.host, testParams.port, 0);

    // create abstraction for Manager
    ProcessManager manager(pgroups, tasks, combiParams, std::move(loadmodel));

    // the combiparameters are sent to all process groups before the
    // computations start
    manager.updateCombiParameters();

    for (unsigned int i = 0; i < testParams.ncombi; i++) {
      if (i == 0) {
        Stats::startEvent("manager run");
        manager.runfirst();
        Stats::stopEvent("manager run");

        // exchange subspace sizes to unify the dsgs with the remote system
        Stats::startEvent("manager unify subspace sizes with remote");
        auto newSGsize = manager.unifySubspaceSizesThirdLevel(thirdLevelExtraSparseGrid);
        Stats::stopEvent("manager unify subspace sizes with remote");
        {
          std::vector<LevelVector> created;
          auto lmaxReduced = testParams.lmax;
          for (auto& l : lmaxReduced) {
            l -= 1;
          }
          combigrid::createTruncatedHierarchicalLevels(lmaxReduced, testParams.lmin, created);
          auto numDofSG = combigrid::getNumDofHierarchical(created, boundary);
          if (thirdLevelExtraSparseGrid) {
            //TODO more accurate check
            BOOST_CHECK_LT(newSGsize, numDofSG);
          } else {
            BOOST_CHECK_EQUAL(newSGsize, numDofSG);
          }
        }
      } else {
        Stats::startEvent("manager run");
        manager.runnext();
        Stats::stopEvent("manager run");
      }
      // combine grids
      Stats::startEvent("manager combine third level");
      // do two TCP-based communications and one file-based one
      if (i < 2) {
        BOOST_CHECK_NO_THROW(manager.combineThirdLevel());
      } else {
        std::string filenamePrefixToWrite = "dsgu_combine_" + std::to_string(testParams.sysNum);
        std::string writeCompleteTokenFileName = filenamePrefixToWrite + "_complete.txt";
        std::string filenamePrefixToRead =
            "dsgu_combine_" + std::to_string((testParams.sysNum + 1) % 2);
        std::string startReadingTokenFileName = filenamePrefixToRead + "_complete.txt";
        manager.combineThirdLevelFileBased(filenamePrefixToWrite, writeCompleteTokenFileName,
                                           filenamePrefixToRead, startReadingTokenFileName);
      }
      // manager.combine();
      Stats::stopEvent("manager combine third level");
    }

    // test Monte-Carlo interpolation
    std::vector<std::vector<real>> interpolationCoords;
    size_t numMCValues = (testParams.dim > 2) ? 1000 : 10000;
    std::vector<CombiDataType> values(numMCValues, 0.0);
    real l2ErrorSingle = 0.;
    // TestFnCount<CombiDataType> initialFunction;
    ParaboloidFn<CombiDataType> initialFunction;

    // compare to third-level monte carlo interpolation
    manager.monteCarloThirdLevel(numMCValues, interpolationCoords, values);
    real l2ErrorTwoSystems = 0.;
    for (size_t i = 0; i < interpolationCoords.size(); ++i) {
      // l2ErrorTwoSystems += std::pow(initialFunction(interpolationCoords[i], testParams.ncombi) - values[i], 2);
      l2ErrorTwoSystems += std::pow(initialFunction(interpolationCoords[i]) - values[i], 2);
    }

    Stats::startEvent("manager interpolate");
    values = manager.interpolateValues(interpolationCoords);
    Stats::stopEvent("manager interpolate");

    for (size_t i = 0; i < interpolationCoords.size(); ++i) {
      // l2ErrorSingle += std::pow(initialFunction(interpolationCoords[i], testParams.ncombi) - values[i], 2);
      l2ErrorSingle += std::pow(initialFunction(interpolationCoords[i]) - values[i], 2);
    }

    std::cout << "Monte carlo errors are " << l2ErrorSingle << " on this system and " <<
      l2ErrorTwoSystems << " in total. boundary: " << boundary << std::endl;
    BOOST_CHECK_LE(l2ErrorTwoSystems, l2ErrorSingle);

    std::string filename("thirdLevel_" + std::to_string(testParams.ncombi) + ".raw");
    Stats::startEvent("manager write solution");
    manager.parallelEval(testParams.lmax, filename, 0);
    manager.writeDSGsToDisk("thirdLevel_" + std::to_string(testParams.ncombi) + "_dsgs");
    manager.readDSGsFromDisk("thirdLevel_" + std::to_string(testParams.ncombi) + "_dsgs");
    Stats::stopEvent("manager write solution");

    manager.exit();

    // if output files are not needed, remove them right away
    remove(("thirdLevel_" + std::to_string(testParams.ncombi) + "_0.raw").c_str());
    remove(("thirdLevel_" + std::to_string(testParams.ncombi) + "_0.raw_header").c_str());
  }
  else {
    ProcessGroupWorker pgroup;
    SignalType signal = -1;
    signal = pgroup.wait();
    // omitting to count RUN_FIRST signal, as it is executed once for every task
    int nrun = 1;
    while (signal != EXIT) {
      BOOST_TEST_CHECKPOINT("Last Successful Worker Signal " + std::to_string(signal));
      BOOST_REQUIRE_NO_THROW(signal = pgroup.wait());
      if (signal == RUN_NEXT) {
        ++nrun;
      }
      // std::cout << "Worker with rank " << theMPISystem()->getLocalRank() << " processed signal "
      //           << signal << std::endl;
      if (signal == COMBINE_THIRD_LEVEL || signal == WAIT_FOR_TL_COMBI_RESULT) {
        // after combination check workers' grids
        BOOST_CHECK(checkReducedFullGrid(pgroup, nrun));
      }
      if (signal == REDUCE_SUBSPACE_SIZES_TL) {
        std::cout << "reduce ";
        for (auto& dsg : pgroup.getCombinedUniDSGVector()) {
          for (size_t size : dsg->getSubspaceDataSizes()) {
            std::cout << size << " ";
          }
        }
        std::cout << std::endl;
      }
      if (signal == REDUCE_SUBSPACE_SIZES_TL_AND_ALLOCATE_EXTRA_SG) {
        std::cout << "reduce extra ";
        for (auto& dsg : pgroup.getExtraUniDSGVector()) {
          for (size_t size : dsg->getSubspaceDataSizes()) {
            std::cout << size << " ";
          }
        }
        std::cout << std::endl;
      }
      if(signal == INIT_DSGUS){
        std::cout << "INIT DSGUS ";
        for (auto& dsg : pgroup.getCombinedUniDSGVector()) {
          for (size_t size : dsg->getSubspaceDataSizes()) {
            std::cout << size << " ";
          }
        }
        std::cout << std::endl;
      }
      TestHelper::testStrayMessages(theMPISystem()->getLocalComm());
    }
    BOOST_CHECK_EQUAL(pgroup.getCurrentNumberOfCombinations(), testParams.ncombi);
    for (const auto& b : pgroup.getCombiParameters().getBoundary())
      BOOST_CHECK_EQUAL(b, testParams.boundary);
    for (const auto& r : pgroup.getCombiParameters().getLMaxReductionVector())
      BOOST_CHECK_EQUAL(r, 1);
  }

  combigrid::Stats::finalize();
  combigrid::Stats::write("stats_thirdLevel_" + std::to_string(testParams.sysNum) + ".json");
  MPI_Barrier(testParams.comm);
  BOOST_CHECK(!TestHelper::testStrayMessages(testParams.comm));
}


/**
 * @brief test for the static task assignment mechanism, both systems read their assignment from file `test_scheme.json`
 */
void testCombineThirdLevelStaticTaskAssignment(TestParams& testParams, bool thirdLevelExtraSparseGrid = false) {
  BOOST_CHECK(testParams.comm != MPI_COMM_NULL);

  combigrid::Stats::initialize();
  theMPISystem()->initWorldReusable(testParams.comm, testParams.ngroup, testParams.nprocs);

  auto loadmodel = std::unique_ptr<LoadModel>(new LinearLoadModel());
  std::vector<BoundaryType> boundary(testParams.dim, testParams.boundary);

  std::vector<LevelVector> levels;
  std::vector<combigrid::real> coeffs;
  std::vector<size_t> taskNumbers; // only used in case of static task assignment
  bool useStaticTaskAssignment = false;
  {
  // read in CT scheme
    std::unique_ptr<CombiMinMaxSchemeFromFile> scheme(
        new CombiMinMaxSchemeFromFile(testParams.dim, testParams.lmin, testParams.lmax, "test_scheme.json"));
    const auto& pgNumbers = scheme->getProcessGroupNumbers();
    if (pgNumbers.size() > 0) {
      useStaticTaskAssignment = true;
      const auto& allCoeffs = scheme->getCoeffs();
      const auto& allLevels = scheme->getCombiSpaces();
      const auto [itMin, itMax] = std::minmax_element(pgNumbers.begin(), pgNumbers.end());
      assert(*itMin == 0);  // make sure it starts with 0
      // filter out only those tasks that belong to "our" process group
      const auto& pgroupNumber = theMPISystem()->getProcessGroupNumber();
      for (size_t taskNo = 0; taskNo < pgNumbers.size(); ++taskNo) {
        if (pgNumbers[taskNo] == pgroupNumber) {
          taskNumbers.push_back(taskNo);
          coeffs.push_back(allCoeffs[taskNo]);
          levels.push_back(allLevels[taskNo]);
        }
      }
      MASTER_EXCLUSIVE_SECTION {
        std::cout << " Process group " << pgroupNumber << " will run " << levels.size() << " of "
                  << pgNumbers.size() << " tasks." << std::endl;
      }
    } else {
      // levels and coeffs are only used in manager
      WORLD_MANAGER_EXCLUSIVE_SECTION {
        coeffs = scheme->getCoeffs();
        levels = scheme->getCombiSpaces();
        std::cout << levels.size() << " tasks to distribute." << std::endl;
      }
    }
  }

  BOOST_REQUIRE(useStaticTaskAssignment);
  BOOST_REQUIRE_EQUAL(levels.size(), coeffs.size());

  WORLD_MANAGER_EXCLUSIVE_SECTION {
    ProcessGroupManagerContainer pgroups;
    for (size_t i = 0; i < testParams.ngroup; ++i) {
      int pgroupRootID((int)i);
      pgroups.emplace_back(std::make_shared<ProcessGroupManager>(pgroupRootID));
    }

    TaskContainer tasks;
    std::vector<size_t> taskIDs;
    for (size_t i = 0; i < levels.size(); i++) {
      Task* t = new TaskConstParaboloid(levels[i], boundary, coeffs[i], loadmodel.get());

      tasks.push_back(t);
      taskIDs.push_back(t->getID());
    }

    // create combiparameters
    std::vector<int> parallelization = {static_cast<int>(testParams.nprocs), 1};
    CombiParameters combiParams(testParams.dim, testParams.lmin, testParams.lmax, boundary, levels,
                                coeffs, taskIDs, testParams.ncombi, 1, parallelization,
                                LevelVector(testParams.dim, 0), LevelVector(testParams.dim, 1),
                                false, testParams.host, testParams.port, 0);

    // create abstraction for Manager
    ProcessManager manager(pgroups, tasks, combiParams, std::move(loadmodel));

    if (useStaticTaskAssignment) {
      // read in CT scheme -- again!
      std::unique_ptr<CombiMinMaxSchemeFromFile> scheme(new CombiMinMaxSchemeFromFile(
          testParams.dim, testParams.lmin, testParams.lmax, "test_scheme.json"));
      const auto& pgNumbers = scheme->getProcessGroupNumbers();
      for (size_t taskNo = 0; taskNo < tasks.size(); ++taskNo) {
        pgroups[pgNumbers[taskNo]]->storeTaskReference(tasks[taskNo]);
      }
    }

    // the combiparameters are sent to all process groups before the
    // computations start
    manager.updateCombiParameters();

    for (unsigned int i = 0; i < testParams.ncombi; i++) {
      if (i == 0) {
        Stats::startEvent("manager no run first");
        BOOST_TEST_CHECKPOINT("manager first runnext");
        manager.runnext();
        BOOST_TEST_CHECKPOINT("manager init dsgus");
        manager.initDsgus();
        Stats::stopEvent("manager no run first");

        // exchange subspace sizes to unify the dsgs with the remote system
        Stats::startEvent("manager unify subspace sizes with remote");
        BOOST_TEST_CHECKPOINT("manager unifySubspaceSizesThirdLevel");
        manager.unifySubspaceSizesThirdLevel(thirdLevelExtraSparseGrid),
        Stats::stopEvent("manager unify subspace sizes with remote");
      } else {
        Stats::startEvent("manager run");
        BOOST_TEST_CHECKPOINT("manager runnext " + std::to_string(i));
        manager.runnext();
        Stats::stopEvent("manager run");
      }
      // combine grids
      Stats::startEvent("manager combine third level");
      BOOST_TEST_CHECKPOINT("manager combineThirdLevel " + std::to_string(i));
      manager.combineThirdLevel();
      Stats::stopEvent("manager combine third level");
    }

    BOOST_TEST_CHECKPOINT("manager exit");
    manager.exit();
  }
  else {
    ProcessGroupWorker pgroup;
    SignalType signal = -1;
    while (signal != EXIT) {
      BOOST_TEST_CHECKPOINT("Last Successful Worker Signal " + std::to_string(signal));
      signal = pgroup.wait();
      // if using static task assignment, we initialize all tasks after the combi parameters are
      // updated
      if (useStaticTaskAssignment) {
        if (signal == UPDATE_COMBI_PARAMETERS) {
          BOOST_TEST_CHECKPOINT("Update for static task assignment " +
                                std::to_string(taskNumbers.size()));
          // initialize all "our" tasks
          for (size_t taskIndex = 0; taskIndex < taskNumbers.size(); ++taskIndex) {
            auto task = new TaskConstParaboloid(levels[taskIndex], boundary, coeffs[taskIndex],
                                                loadmodel.get());
            task->setID(taskNumbers[taskIndex]);
            pgroup.initializeTaskAndFaults(task);
          }
        }
        if (signal == RUN_FIRST) {
          BOOST_CHECK(false);
        }
      }
    }
    BOOST_TEST_CHECKPOINT("exited worker main loop");
  }
  combigrid::Stats::finalize();
  combigrid::Stats::write("stats_thirdLevel_static_" + std::to_string(testParams.sysNum) + ".json");
  MPI_Barrier(testParams.comm);
  BOOST_CHECK(!TestHelper::testStrayMessages(testParams.comm));
}

void testPretendThirdLevel(TestParams& testParams) {
  BOOST_CHECK(testParams.comm != MPI_COMM_NULL);

  combigrid::Stats::initialize();

  size_t procsPerSys = testParams.ngroup * testParams.nprocs + 1;

  theMPISystem()->initWorldReusable(testParams.comm, testParams.ngroup, testParams.nprocs);

  WORLD_MANAGER_EXCLUSIVE_SECTION {
    ProcessGroupManagerContainer pgroups;

    auto loadmodel = std::unique_ptr<LoadModel>(new LinearLoadModel());
    std::vector<BoundaryType> boundary(testParams.dim, testParams.boundary);

    // create third level specific scheme
    CombiMinMaxScheme combischeme(testParams.dim, testParams.lmin, testParams.lmax);
    combischeme.createClassicalCombischeme();
    // combischeme.createAdaptiveCombischeme();
    // get full scheme first
    std::vector<LevelVector> fullLevels = combischeme.getCombiSpaces();
    std::vector<combigrid::real> fullCoeffs = combischeme.getCoeffs();
    // split scheme and assign each half to a system
    std::vector<LevelVector> levels;
    std::vector<combigrid::real> coeffs;
    CombiThirdLevelScheme::createThirdLevelScheme(fullLevels, fullCoeffs, testParams.sysNum, 2,
                                                  levels, coeffs);

    BOOST_REQUIRE_EQUAL(levels.size(), coeffs.size());
    TaskContainer tasks;
    std::vector<size_t> taskIDs;
    auto reduceCombinationDimsLmax = LevelVector(testParams.dim, 1);

    // create combiparameters
    std::vector<int> parallelization = {static_cast<int>(testParams.nprocs), 1};
    CombiParameters combiParams(
        testParams.dim, testParams.lmin, testParams.lmax, boundary, levels, coeffs, taskIDs,
        testParams.ncombi, 1, parallelization, LevelVector(testParams.dim, 0),
        reduceCombinationDimsLmax, true, testParams.host, testParams.port, 0);

    auto decomposition = combigrid::getStandardDecomposition(testParams.lmax, parallelization);
    combiParams.setDecomposition(decomposition);
    // create abstraction for Manager
    ProcessManager manager(pgroups, tasks, combiParams, std::move(loadmodel));

    // the combiparameters are sent to all process groups before the
    // computations start
    manager.updateCombiParameters();
    auto dsguSizes =
        getPartitionedNumDOFSGAdaptive(testParams.lmin, testParams.lmax - reduceCombinationDimsLmax,
                                       testParams.lmax, decomposition);

    BOOST_TEST_CHECKPOINT("partitioned num dof " + std::to_string(dsguSizes[0]));
    for (unsigned int i = 0; i < testParams.ncombi; i++) {
      BOOST_TEST_CHECKPOINT("pretend combine for broker " + std::to_string(i));
      auto numErrors = manager.pretendCombineThirdLevelForBroker(dsguSizes, true);
      BOOST_CHECK_EQUAL(numErrors, 0);
      // BOOST_TEST_CHECKPOINT("pretend combine for workers " + std::to_string(i));
      // manager.pretendCombineThirdLevelForWorkers(dsguSizes, true);
      // BOOST_CHECK_EQUAL(numErrors, 0);
    }
    manager.exit();
  }
  else {
    // do nothing
    BOOST_TEST_CHECKPOINT("exited worker main loop");
  }
  combigrid::Stats::finalize();
  MPI_Barrier(testParams.comm);
  BOOST_CHECK(!TestHelper::testStrayMessages(testParams.comm));
}

BOOST_FIXTURE_TEST_SUITE(thirdLevel, TestHelper::BarrierAtEnd, *boost::unit_test::timeout(600))

BOOST_AUTO_TEST_CASE(test_0, *boost::unit_test::tolerance(TestHelper::tolerance)) {
  unsigned int numSystems = 2;
  unsigned int ngroup = 1;
  unsigned int nprocs = 1;
  unsigned int ncombi = 3;
  DimType dim = 2;
  LevelVector lmin(dim, 1);
  LevelVector lmax(dim, 2);

  unsigned int sysNum;
  CommunicatorType newcomm;

  for (auto boundary : std::vector<BoundaryType>({0, 1, 2})) {
    assignProcsToSystems(ngroup, nprocs, numSystems, sysNum, newcomm);

    if (newcomm != MPI_COMM_NULL) {  // remove unnecessary procs
      TestParams testParams(dim, lmin, lmax, boundary, ngroup, nprocs, ncombi, sysNum, newcomm);
      startInfrastructure();
      testCombineThirdLevel(testParams, false);
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }
}

BOOST_AUTO_TEST_CASE(test_2, *boost::unit_test::tolerance(TestHelper::tolerance)) {
  unsigned int numSystems = 2;
  unsigned int ngroup = 1;
  unsigned int nprocs = 1;
  unsigned int ncombi = 10;
  DimType dim = 2;
  LevelVector lmin(dim, 2);
  LevelVector lmax(dim, 3);

  unsigned int sysNum;
  CommunicatorType newcomm;

  for (auto boundary : std::vector<BoundaryType>({2})) {
    assignProcsToSystems(ngroup, nprocs, numSystems, sysNum, newcomm);

    if (newcomm != MPI_COMM_NULL) {  // remove unnecessary procs
      TestParams testParams(dim, lmin, lmax, boundary, ngroup, nprocs, ncombi, sysNum, newcomm);
      startInfrastructure();
      testCombineThirdLevel(testParams, false);
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }
}

BOOST_AUTO_TEST_CASE(test_3, *boost::unit_test::tolerance(TestHelper::tolerance)) {
  unsigned int numSystems = 2;
  unsigned int ngroup = 1;
  unsigned int nprocs = 1;
  unsigned int ncombi = 10;
  DimType dim = 2;
  LevelVector lmin(dim, 4);
  LevelVector lmax(dim, 7);

  unsigned int sysNum;
  CommunicatorType newcomm;

  for (auto boundary : std::vector<BoundaryType>({0, 1, 2})) {
    for (bool extraSparseGrid : {false, true}) {
      assignProcsToSystems(ngroup, nprocs, numSystems, sysNum, newcomm);

      if (newcomm != MPI_COMM_NULL) {  // remove unnecessary procs
        TestParams testParams(dim, lmin, lmax, boundary, ngroup, nprocs, ncombi, sysNum, newcomm);
        startInfrastructure();
        testCombineThirdLevel(testParams, extraSparseGrid);
      }

      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
}

BOOST_AUTO_TEST_CASE(test_4, *boost::unit_test::tolerance(TestHelper::tolerance)) {
  unsigned int numSystems = 2;
  unsigned int ngroup = 2;
  unsigned int nprocs = 1;
  unsigned int ncombi = 10;
  DimType dim = 2;
  LevelVector lmin(dim, 4);
  LevelVector lmax(dim, 7);

  unsigned int sysNum;
  CommunicatorType newcomm;

  for (auto boundary : std::vector<BoundaryType>({0, 1, 2})) {
    assignProcsToSystems(ngroup, nprocs, numSystems, sysNum, newcomm);

    if (newcomm != MPI_COMM_NULL) {  // remove unnecessary procs
      TestParams testParams(dim, lmin, lmax, boundary, ngroup, nprocs, ncombi, sysNum, newcomm);
      startInfrastructure();
      testCombineThirdLevel(testParams, false);
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }
}

BOOST_AUTO_TEST_CASE(test_5, *boost::unit_test::tolerance(TestHelper::tolerance)) {
  unsigned int numSystems = 2;
  unsigned int ngroup = 1;
  unsigned int nprocs = 2;
  unsigned int ncombi = 10;
  DimType dim = 2;
  LevelVector lmin(dim, 4);
  LevelVector lmax(dim, 7);

  unsigned int sysNum;
  CommunicatorType newcomm;

  for (auto boundary : std::vector<BoundaryType>({0, 1, 2})) {
    assignProcsToSystems(ngroup, nprocs, numSystems, sysNum, newcomm);

    if (newcomm != MPI_COMM_NULL) {  // remove unnecessary procs
      TestParams testParams(dim, lmin, lmax, boundary, ngroup, nprocs, ncombi, sysNum, newcomm);
      startInfrastructure();
      testCombineThirdLevel(testParams, false);
      MPI_Barrier(newcomm);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

// like test_5, but with static group assignment
BOOST_AUTO_TEST_CASE(test_6, *boost::unit_test::tolerance(TestHelper::tolerance)) {
  unsigned int numSystems = 2;
  unsigned int ngroup = 3;
  unsigned int nprocs = 1;
  unsigned int ncombi = 10;
  DimType dim = 2;
  LevelVector lmin = {3,6};
  LevelVector lmax = {7,10};

  unsigned int sysNum;
  CommunicatorType newcomm;

  for (auto boundary : std::vector<BoundaryType>({0, 1, 2})) {
    assignProcsToSystems(ngroup, nprocs, numSystems, sysNum, newcomm);
    BOOST_TEST_CHECKPOINT("static group assignment. sysNum: " + std::to_string(sysNum));
    if (newcomm != MPI_COMM_NULL) {  // remove unnecessary procs
      BOOST_TEST_CHECKPOINT("static sysNum: " + std::to_string(sysNum));
      for (bool extraSparseGrid : {false, true}) {
        TestParams testParams(dim, lmin, lmax, boundary, ngroup, nprocs, ncombi, sysNum, newcomm);
        startInfrastructure();
        testCombineThirdLevelStaticTaskAssignment(testParams, extraSparseGrid);
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }
}

// like test_5, but send only dummy data between managers
BOOST_AUTO_TEST_CASE(test_7, *boost::unit_test::tolerance(TestHelper::tolerance)) {
  unsigned int numSystems = 2;
  unsigned int ngroup = 1;
  unsigned int ncombi = 10;
  DimType dim = 2;
  LevelVector lmin = {3, 6};
  LevelVector lmax = {7, 10};

  unsigned int sysNum;
  CommunicatorType newcomm;

  for (auto boundary : std::vector<BoundaryType>({2})) {
    for (unsigned int nprocs : {1, 2, 3}) {
      assignProcsToSystems(ngroup, nprocs, numSystems, sysNum, newcomm);

      if (newcomm != MPI_COMM_NULL) {  // remove unnecessary procs
        TestParams testParams(dim, lmin, lmax, boundary, ngroup, nprocs, ncombi, sysNum, newcomm);
        startInfrastructure();
        testPretendThirdLevel(testParams);
      }

      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
}

// "target" case: 6-dim, one-sided boundary, extra sparse grid, many process groups
BOOST_AUTO_TEST_CASE(test_8, *boost::unit_test::tolerance(TestHelper::tolerance)) {
  unsigned int numSystems = 2;
  unsigned int nprocs = 1;
  unsigned int ncombi = 10;
  DimType dim = 6;
  BoundaryType boundary = 2;
  // LevelVector lmin(dim, 4);
  // LevelVector lmax(dim, 7);
  LevelVector lmin(dim, 1);
  LevelVector lmax(dim, 4);

  // the fact that values are wrong when calling this for ngroup = 3
  // but not for other values is a bug that I dont understand yet
  for (unsigned int ngroup : std::vector<unsigned int>({1, 2, 3})) {
    unsigned int sysNum;
    CommunicatorType newcomm = MPI_COMM_NULL;
    assignProcsToSystems(ngroup, nprocs, numSystems, sysNum, newcomm);

    if (newcomm != MPI_COMM_NULL) {  // remove unnecessary procs
      TestParams testParams(dim, lmin, lmax, boundary, ngroup, nprocs, ncombi, sysNum, newcomm);
      startInfrastructure();
      BOOST_CHECK_NO_THROW(testCombineThirdLevel(testParams, true));
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }
}

BOOST_AUTO_TEST_SUITE_END()
