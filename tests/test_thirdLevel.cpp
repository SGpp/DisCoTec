#define BOOST_TEST_DYN_LINK
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <boost/serialization/export.hpp>
#include <boost/test/unit_test.hpp>
#include <cstdio>
#include <filesystem>
#include <thread>

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
#include "utils/DecompositionUtils.hpp"
#include "utils/MonteCarlo.hpp"
#include "utils/Types.hpp"
#include "stdlib.h"
#include "test_helper.hpp"

using namespace combigrid;

BOOST_CLASS_EXPORT(TaskConstParaboloid)
BOOST_CLASS_EXPORT(TaskCount)

#define CHECK_NO_THROW_LOUD(statement)                \
  BOOST_CHECK_NO_THROW(                               \
      try { statement; } catch (std::exception & e) { \
        std::cout << e.what() << std::endl;           \
        throw e;                                      \
      })

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
  unsigned int numSystems = 2;
  const CommunicatorType& comm;
  std::string host = "localhost";
  unsigned short port = 9999;

  TestParams(DimType dim, LevelVector& lmin, LevelVector& lmax, BoundaryType boundary,
             unsigned int ngroup, unsigned int nprocs, unsigned int ncombi, unsigned int sysNum,
             const CommunicatorType& comm, unsigned int numSystems = 2,
             const std::string& host = "localhost",
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
        numSystems(numSystems),
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
  const auto& tasks = worker.getTasks();
  int numGrids = (int)worker.getCombiParameters().getNumGrids();

  BOOST_CHECK(tasks.size() > 0);
  BOOST_CHECK(numGrids > 0);

  // to check if any data was actually compared
  bool any = false;

  for (const auto& t : tasks) {
    for (int g = 0; g < numGrids; g++) {
      const DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid(g);
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

void assignProcsToSystems(unsigned int procsPerSys, unsigned int numSystems, unsigned int& sysNum,
                          CommunicatorType& newcomm) {
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
  sleep(5);
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
  BOOST_CHECK(thirdLevelExtraSparseGrid == false);
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
                                coeffs, taskIDs, testParams.ncombi, 1,
                                CombinationVariant::sparseGridReduce, parallelization,
                                LevelVector(testParams.dim, 0), LevelVector(testParams.dim, 1), 32,
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
            // TODO more accurate check
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
      // // do two TCP-based communications and one or more file-based ones
      // if (i < 2) {
        BOOST_CHECK_NO_THROW(manager.combineThirdLevel());
//       } else {
//         // remove previously interpolated files
//         auto status =
//             system(("rm tl_group" + std::to_string(testParams.sysNum) + "_diagonal*").c_str());
//         BOOST_WARN_GE(status, 0);
//         // write sparse grid data
//         std::string filenamePrefixToWrite = "dsgu_combine_" + std::to_string(testParams.sysNum);
//         std::string writeCompleteTokenFileName = filenamePrefixToWrite + "_complete.txt";
//         std::string filenamePrefixToRead =
//             "dsgu_combine_" + std::to_string((testParams.sysNum + 1) % 2);
//         std::string startReadingTokenFileName = filenamePrefixToRead + "_complete.txt";
//         manager.combineThirdLevelFileBasedWrite(filenamePrefixToWrite, writeCompleteTokenFileName);
// #ifdef DISCOTEC_USE_HIGHFIVE
//         // write interpolated values
//         std::vector<std::vector<real>> interpolationCoords = {
//             std::vector<real>(testParams.dim, 0.001), std::vector<real>(testParams.dim, 0.2),
//             std::vector<real>(testParams.dim, 0.5), std::vector<real>(testParams.dim, 0.7),
//             std::vector<real>(testParams.dim, 0.999)};
//         manager.writeInterpolatedValuesSingleFile(
//             interpolationCoords, "tl_group" + std::to_string(testParams.sysNum) + "_diagonal");
// #else   // def DISCOTEC_USE_HIGHFIVE
//         sleep(1);
// #endif  // def DISCOTEC_USE_HIGHFIVE
//         manager.combineThirdLevelFileBasedReadReduce(filenamePrefixToRead,
//                                                      startReadingTokenFileName);
//       }
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
      // l2ErrorSingle += std::pow(initialFunction(interpolationCoords[i], testParams.ncombi) -
      // values[i], 2);
      l2ErrorSingle += std::pow(initialFunction(interpolationCoords[i]) - values[i], 2);
    }

    std::cout << "Monte carlo errors are " << l2ErrorSingle << " on this system and "
              << l2ErrorTwoSystems << " in total. boundary: " << boundary << std::endl;
    BOOST_CHECK_LE(l2ErrorTwoSystems, l2ErrorSingle);

    std::string filename("thirdLevel_" + std::to_string(testParams.ncombi) + ".raw");
    Stats::startEvent("manager write solution");
    std::string dsgFileName = "thirdLevel_" + std::to_string(testParams.ncombi) + "_" +
                              std::to_string(testParams.sysNum) + "_dsgs";
    manager.writeDSGsToDisk(dsgFileName);
    manager.parallelEval(testParams.lmax, filename, 0);
    manager.waitForAllGroupsToWait();
    BOOST_TEST_CHECKPOINT("Waiting for file system...");
    while (!std::filesystem::exists(dsgFileName + "_0")) {
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    manager.readDSGsFromDisk(dsgFileName);
    Stats::stopEvent("manager write solution");

    manager.waitForAllGroupsToWait();
    manager.exit();

    // if output files are not needed, remove them right away
    remove((dsgFileName + "_0").c_str());
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
      BOOST_REQUIRE_NO_THROW(
          try { signal = pgroup.wait(); } catch (std::exception& e) {
            std::cout << "Exception in worker wait after " << std::to_string(signal) << ": "
                      << e.what() << std::endl;
            throw e;
          });
      if (signal == RUN_NEXT) {
        ++nrun;
      }
      // std::cout << "Worker with rank " << theMPISystem()->getLocalRank() << " processed signal "
      //           << signal << std::endl;
      if (signal == COMBINE_THIRD_LEVEL || signal == WAIT_FOR_TL_COMBI_RESULT ||
          signal == COMBINE_READ_DSGS_AND_REDUCE) {
        // after combination check workers' grids
        BOOST_CHECK(checkReducedFullGrid(pgroup, nrun));
      }
      if (signal == COMBINE_WRITE_DSGS) {
        // write partial stats (only one process group per system will get this signal)
        Stats::writePartial("third_level_partial" + std::to_string(testParams.sysNum) + ".json",
                            theMPISystem()->getLocalComm());
      }
      if (signal == REDUCE_SUBSPACE_SIZES_TL) {
        std::cout << "reduce ";
        for (auto& dsg : pgroup.getCombinedDSGVector()) {
          for (size_t size : dsg->getSubspaceDataSizes()) {
            std::cout << size << " ";
          }
        }
        std::cout << std::endl;
      }
      if (signal == REDUCE_SUBSPACE_SIZES_TL_AND_ALLOCATE_EXTRA_SG) {
        std::cout << "reduce extra ";
        for (auto& dsg : pgroup.getExtraDSGVector()) {
          for (size_t size : dsg->getSubspaceDataSizes()) {
            std::cout << size << " ";
          }
        }
        std::cout << std::endl;
      }
      if(signal == INIT_DSGUS){
        std::cout << "INIT DSGUS ";
        for (auto& dsg : pgroup.getCombinedDSGVector()) {
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
  BOOST_CHECK(thirdLevelExtraSparseGrid == false);
  BOOST_CHECK(testParams.comm != MPI_COMM_NULL);

  combigrid::Stats::initialize();
  theMPISystem()->initWorldReusable(testParams.comm, testParams.ngroup, testParams.nprocs);

  auto loadmodel = std::unique_ptr<LoadModel>(new LinearLoadModel());
  std::vector<BoundaryType> boundary(testParams.dim, testParams.boundary);

  std::vector<LevelVector> levels;
  std::vector<combigrid::real> coeffs;
  std::vector<size_t> taskNumbers;  // only used in case of static task assignment
  bool useStaticTaskAssignment = false;
  {
    // read in CT scheme
    std::unique_ptr<CombiMinMaxSchemeFromFile> scheme(new CombiMinMaxSchemeFromFile(
        testParams.dim, testParams.lmin, testParams.lmax, "test_scheme.json", MPI_COMM_SELF));
    const auto& pgNumbers = scheme->getProcessGroupNumbers();
    if (pgNumbers.size() > 0) {
      useStaticTaskAssignment = true;
      auto pgroupNumber = theMPISystem()->getProcessGroupNumber();
      size_t totalNumTasks =
          combigrid::getAssignedLevels(*scheme, pgroupNumber, levels, coeffs, taskNumbers);
      BOOST_CHECK_EQUAL(totalNumTasks, pgNumbers.size());
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
                                coeffs, taskIDs, testParams.ncombi, 1,
                                CombinationVariant::sparseGridReduce, parallelization,
                                LevelVector(testParams.dim, 0), LevelVector(testParams.dim, 1), 32,
                                false, testParams.host, testParams.port, 0);

    // create abstraction for Manager
    ProcessManager manager(pgroups, tasks, combiParams, std::move(loadmodel));

    if (useStaticTaskAssignment) {
      // read in CT scheme -- again!
      std::unique_ptr<CombiMinMaxSchemeFromFile> scheme(new CombiMinMaxSchemeFromFile(
          testParams.dim, testParams.lmin, testParams.lmax, "test_scheme.json", MPI_COMM_SELF));
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
        manager.unifySubspaceSizesThirdLevel(thirdLevelExtraSparseGrid);
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
          pgroup.initializeAllTasks<TaskConstParaboloid>(levels, coeffs, taskNumbers,
                                                         loadmodel.get());
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

void testCombineThirdLevelWithoutManagers(
    TestParams& testParams, CombinationVariant variant = CombinationVariant::sparseGridReduce) {
  bool thirdLevelExtraSparseGrid = true;  // TODO make parameter again
  BOOST_CHECK(testParams.comm != MPI_COMM_NULL);

  combigrid::Stats::initialize();
  theMPISystem()->initWorldReusable(testParams.comm, testParams.ngroup, testParams.nprocs, false,
                                    true);

  if (testParams.nprocs > 1 && thirdLevelExtraSparseGrid) {
    theMPISystem()->initOuputGroupComm(testParams.nprocs / 2);
  }
  auto loadmodel = std::unique_ptr<LoadModel>(new LinearLoadModel());
  std::vector<BoundaryType> boundary(testParams.dim, testParams.boundary);

  std::vector<LevelVector> levels;
  std::vector<combigrid::real> coeffs;
  std::vector<size_t> taskNumbers;  // only used in case of static task assignment
  {  // this is scoped so that anything created in here frees up memory for computation
    CombiMinMaxScheme combischeme(testParams.dim, testParams.lmin, testParams.lmax);
    combischeme.createAdaptiveCombischeme();
    // get full scheme first
    std::vector<LevelVector> fullLevels = combischeme.getCombiSpaces();
    std::vector<combigrid::real> fullCoeffs = combischeme.getCoeffs();
    MIDDLE_PROCESS_EXCLUSIVE_SECTION {
      BOOST_TEST_MESSAGE("total: " + std::to_string(fullLevels.size()) + " tasks");
    }

    // split scheme and assign each part to a system
    std::vector<LevelVector> systemLevels;
    std::vector<combigrid::real> systemCoeffs;
    std::vector<combigrid::real> fractions(testParams.numSystems,
                                           (1. / static_cast<double>(testParams.numSystems)));
    CombiThirdLevelScheme::createThirdLevelScheme(fullLevels, fullCoeffs, testParams.sysNum,
                                                  testParams.numSystems, systemLevels, systemCoeffs,
                                                  fractions);

    BOOST_REQUIRE_EQUAL(systemLevels.size(), systemCoeffs.size());

    for (size_t i = 0; i < systemLevels.size(); ++i) {
      // assign round-robin to process groups
      if (i % theMPISystem()->getNumGroups() == theMPISystem()->getProcessGroupNumber()) {
        // find index in full list
        auto position = std::find(systemLevels.begin(), systemLevels.end(), systemLevels[i]);
        BOOST_REQUIRE(position != systemLevels.end());
        taskNumbers.push_back(std::distance(systemLevels.begin(), position));
      }
    }
    for (size_t i = 0; i < systemLevels.size(); ++i) {
      if (std::find(taskNumbers.begin(), taskNumbers.end(), i) != taskNumbers.end()) {
        levels.push_back(systemLevels[i]);
        coeffs.push_back(systemCoeffs[i]);
      }
    }
    MASTER_EXCLUSIVE_SECTION {
      BOOST_TEST_MESSAGE("group " + std::to_string(theMPISystem()->getProcessGroupNumber()) +
                         " on system " + std::to_string(testParams.sysNum) + " will run " +
                         std::to_string(levels.size()) + " of " +
                         std::to_string(systemLevels.size()) + " tasks.");
    }
  }
  BOOST_REQUIRE_EQUAL(levels.size(), coeffs.size());

  WORLD_MANAGER_EXCLUSIVE_SECTION { BOOST_CHECK(false); }
  ProcessGroupWorker worker;
  // parallelization needs to be the same as in TaskConstParaboloid::init()
  std::vector<int> parallelization(testParams.dim, 1);
  parallelization[1] = static_cast<int>(testParams.nprocs);
  // create combiparameters
  CombiParameters combiParams(
      testParams.dim, testParams.lmin, testParams.lmax, boundary, testParams.ncombi, 1, variant,
      parallelization, LevelVector(testParams.dim, 0), LevelVector(testParams.dim, 1), 32, false);
  worker.setCombiParameters(std::move(combiParams));
  BOOST_CHECK_EQUAL(worker.getCombiParameters().getChunkSizeInMebibybtePerThread(), 32);

  // create Tasks
  worker.initializeAllTasks<TaskConstParaboloid>(levels, coeffs, taskNumbers, loadmodel.get());

  BOOST_TEST_CHECKPOINT("initialize sparse grid");
  worker.initCombinedDSGVector();
  std::string writeSubspaceSizeFile = "worker_subspace_sizes_" + std::to_string(testParams.sysNum);
  std::string writeSubspaceSizeFileToken = writeSubspaceSizeFile + "_token.txt";
  std::vector<std::string> readSubspaceSizeFiles;
  std::vector<std::string> readSubspaceSizeFileTokens;
  for (size_t i = 0; i < testParams.numSystems; ++i) {
    if (i != testParams.sysNum) {
      readSubspaceSizeFiles.push_back("worker_subspace_sizes_" + std::to_string(i));
      readSubspaceSizeFileTokens.push_back(readSubspaceSizeFiles.back() + "_token.txt");
    }
  }
  worker.reduceExtraSubspaceSizesFileBased(writeSubspaceSizeFile, writeSubspaceSizeFileToken,
                                           readSubspaceSizeFiles, readSubspaceSizeFileTokens);
  // remove subspace size files to avoid interference between multiple calls to this test function
  MPI_Barrier(testParams.comm);
  // // output sparse grid sizes -- for visual inspection if they are exactly the same between
  // systems
  // // now
  // OUTPUT_GROUP_EXCLUSIVE_SECTION {
  //   if (thirdLevelExtraSparseGrid) {
  //     for (auto& dsg : worker.getExtraDSGVector()) {
  //       for (size_t size : dsg->getSubspaceDataSizes()) {
  //         std::cout << size << " ";
  //       }
  //     }
  //   } else {
  //     for (auto& dsg : worker.getCombinedDSGVector()) {
  //       for (size_t size : dsg->getSubspaceDataSizes()) {
  //         std::cout << size << " ";
  //       }
  //     }
  //   }
  //   std::cout << std::endl;
  // }
  if (thirdLevelExtraSparseGrid) {
    OUTPUT_GROUP_EXCLUSIVE_SECTION {
      // make sure the first few subspaces are populated
      // (minus the first three, which may be 0 depending on boundary condition and decomposition)
      for (auto& dsg : worker.getExtraDSGVector()) {
        for (size_t i = 3; i < 5; ++i) {
          BOOST_CHECK_GT(dsg->getSubspaceDataSizes()[i], 0);
        }
      }
    }
    else {
      BOOST_CHECK_EQUAL(worker.getExtraDSGVector().size(), 0);
    }
  } else {
    // make sure all subspaces are populated on every rank in case of sparse grid reduce
    for (auto& dsg : worker.getCombinedDSGVector()) {
      if (worker.getCombiParameters().getCombinationVariant() ==
          CombinationVariant::sparseGridReduce) {
        for (size_t size : dsg->getSubspaceDataSizes()) {
          BOOST_CHECK_GT(size, 0);
        }
      }
    }
  }

  BOOST_CHECK_EQUAL(worker.getCurrentNumberOfCombinations(), 0);
  BOOST_TEST_CHECKPOINT("run first");
  auto start = std::chrono::high_resolution_clock::now();
  worker.runAllTasks();
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  MASTER_EXCLUSIVE_SECTION {
    BOOST_TEST_MESSAGE("worker run first solver step: " << duration.count() << " milliseconds");
  }

  for (size_t it = 0; it < testParams.ncombi - 1; ++it) {
    std::string writeSparseGridFile =
        "worker_sg_" + std::to_string(testParams.sysNum) + "_step_" + std::to_string(it);
    std::string writeSparseGridFileToken = writeSparseGridFile + "_token.txt";
    std::vector<std::string> readSparseGridFiles;
    std::vector<std::string> readSparseGridFileTokens;

    for (size_t i = 0; i < testParams.numSystems; ++i) {
      if (i != testParams.sysNum) {
        readSparseGridFiles.push_back("worker_sg_" + std::to_string(i) + "_step_" +
                                      std::to_string(it));
        readSparseGridFileTokens.push_back(readSparseGridFiles.back() + "_token.txt");
      }
    }
    BOOST_TEST_CHECKPOINT("combine system-wide");
    start = std::chrono::high_resolution_clock::now();
    worker.combineSystemWideAndWrite(writeSparseGridFile, writeSparseGridFileToken);

    // everyone writes partial stats (to show that we can continue to do stuff while waiting)
    Stats::writePartial("stats_thirdLevel_worker_" + std::to_string(testParams.sysNum) + ".json",
                        theMPISystem()->getWorldComm());
    BOOST_TEST_CHECKPOINT("combine read/reduce");
    std::cout << readSparseGridFiles << std::endl;
    worker.combineReadDistributeSystemWide(readSparseGridFiles, readSparseGridFileTokens, false,
                                           false);

    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    MASTER_EXCLUSIVE_SECTION {
      BOOST_TEST_MESSAGE("worker combine: " << duration.count() << " milliseconds");
    }
    BOOST_CHECK_EQUAL(worker.getCurrentNumberOfCombinations(), it + 1);
    // after combination check workers' grids
    BOOST_CHECK(checkReducedFullGrid(worker, worker.getCurrentNumberOfCombinations()));

    BOOST_TEST_CHECKPOINT("run next");
    worker.runAllTasks();
    MASTER_EXCLUSIVE_SECTION {
      BOOST_TEST_MESSAGE("worker run: " << Stats::getDuration("run") << " milliseconds");
    }
  }
  BOOST_TEST_CHECKPOINT("worker combine last time");
  std::string writeSparseGridFile = "worker_sg_" + std::to_string(testParams.sysNum) + "_final";
  std::string writeSparseGridFileToken = writeSparseGridFile + "_token.txt";
  std::string readSparseGridFile =
      "worker_sg_" + std::to_string((testParams.sysNum + 1) % 2) + "_final";
  std::string readSparseGridFileToken = readSparseGridFile + "_token.txt";
  // TODO combine

  OUTPUT_GROUP_EXCLUSIVE_SECTION {
    MASTER_EXCLUSIVE_SECTION {
      // remove some files we just read on this "system"
      int indexToDelete = (testParams.sysNum + 1) % testParams.numSystems;
      remove(("worker_subspace_sizes_" + std::to_string(indexToDelete)).c_str());
      remove(("worker_subspace_sizes_" + std::to_string(indexToDelete) + "_token.txt").c_str());
    }
  }

  // TODO monte-carlo interpolation? would need to read interpolated values from other system...
  BOOST_TEST_CHECKPOINT("worker exit");
  worker.exit();

  combigrid::Stats::finalize();
  combigrid::Stats::write("stats_thirdLevel_worker_" + std::to_string(testParams.sysNum) + ".json");
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
    CombiParameters combiParams(testParams.dim, testParams.lmin, testParams.lmax, boundary, levels,
                                coeffs, taskIDs, testParams.ncombi, 1,
                                CombinationVariant::sparseGridReduce, parallelization,
                                LevelVector(testParams.dim, 0), reduceCombinationDimsLmax, 32, true,
                                testParams.host, testParams.port, 0);

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

BOOST_FIXTURE_TEST_SUITE(thirdLevel, TestHelper::BarrierAtEnd, *boost::unit_test::timeout(2500))

BOOST_AUTO_TEST_CASE(test_0, *boost::unit_test::tolerance(TestHelper::tolerance) *
                                 boost::unit_test::disabled()) {
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
    assignProcsToSystems(ngroup * nprocs + 1, numSystems, sysNum, newcomm);

    if (newcomm != MPI_COMM_NULL) {  // remove unnecessary procs
      TestParams testParams(dim, lmin, lmax, boundary, ngroup, nprocs, ncombi, sysNum, newcomm);
      startInfrastructure();
      testCombineThirdLevel(testParams, false);
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }
}

BOOST_AUTO_TEST_CASE(test_2, *boost::unit_test::tolerance(TestHelper::tolerance) *
                                 boost::unit_test::disabled()) {
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
    assignProcsToSystems(ngroup * nprocs + 1, numSystems, sysNum, newcomm);

    if (newcomm != MPI_COMM_NULL) {  // remove unnecessary procs
      TestParams testParams(dim, lmin, lmax, boundary, ngroup, nprocs, ncombi, sysNum, newcomm);
      startInfrastructure();
      testCombineThirdLevel(testParams, false);
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }
}

BOOST_AUTO_TEST_CASE(test_3, *boost::unit_test::tolerance(TestHelper::tolerance) *
                                 boost::unit_test::disabled()) {
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
    for (bool extraSparseGrid : {false}) { //TODO add true again
      assignProcsToSystems(ngroup * nprocs + 1, numSystems, sysNum, newcomm);

      if (newcomm != MPI_COMM_NULL) {  // remove unnecessary procs
        TestParams testParams(dim, lmin, lmax, boundary, ngroup, nprocs, ncombi, sysNum, newcomm);
        startInfrastructure();
        testCombineThirdLevel(testParams, extraSparseGrid);
      }

      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
}

BOOST_AUTO_TEST_CASE(test_4, *boost::unit_test::tolerance(TestHelper::tolerance) *
                                 boost::unit_test::disabled()) {
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
    assignProcsToSystems(ngroup * nprocs + 1, numSystems, sysNum, newcomm);

    if (newcomm != MPI_COMM_NULL) {  // remove unnecessary procs
      TestParams testParams(dim, lmin, lmax, boundary, ngroup, nprocs, ncombi, sysNum, newcomm);
      startInfrastructure();
      testCombineThirdLevel(testParams, false);
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }
}

BOOST_AUTO_TEST_CASE(test_5, *boost::unit_test::tolerance(TestHelper::tolerance) *
                                 boost::unit_test::disabled()) {
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
    assignProcsToSystems(ngroup * nprocs + 1, numSystems, sysNum, newcomm);

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
BOOST_AUTO_TEST_CASE(test_6, *boost::unit_test::tolerance(TestHelper::tolerance) *
                                 boost::unit_test::disabled()) {
  unsigned int numSystems = 2;
  unsigned int ngroup = 3;
  unsigned int nprocs = 1;
  unsigned int ncombi = 10;
  DimType dim = 2;
  LevelVector lmin = {3, 6};
  LevelVector lmax = {7, 10};

  unsigned int sysNum;
  CommunicatorType newcomm;

  for (auto boundary : std::vector<BoundaryType>({0, 1, 2})) {
    assignProcsToSystems(ngroup * nprocs + 1, numSystems, sysNum, newcomm);
    BOOST_TEST_CHECKPOINT("static group assignment. sysNum: " + std::to_string(sysNum));
    if (newcomm != MPI_COMM_NULL) {  // remove unnecessary procs
      BOOST_TEST_CHECKPOINT("static sysNum: " + std::to_string(sysNum));
      for (bool extraSparseGrid : {false}) { //TODO add true again
        TestParams testParams(dim, lmin, lmax, boundary, ngroup, nprocs, ncombi, sysNum, newcomm);
        startInfrastructure();
        testCombineThirdLevelStaticTaskAssignment(testParams, extraSparseGrid);
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }
}

// only with workers (which needs static task assignment)
BOOST_AUTO_TEST_CASE(test_workers_only, *boost::unit_test::tolerance(TestHelper::tolerance)) {
  unsigned int numSystems = 2;
  unsigned int ncombi = 10;
  DimType dim = 2;
  LevelVector lmin = {3, 6};
  LevelVector lmax = {7, 10};

  unsigned int sysNum;
  CommunicatorType newcomm;
  for (auto ngroup : std::vector<unsigned int>({1, 2, 3})) {
    for (auto nprocs : std::vector<unsigned int>({1, 2, 4})) {
      if (ngroup * nprocs < 4) {
        for (auto boundary : std::vector<BoundaryType>({1, 2})) {  // TODO no-boundary case
          assignProcsToSystems(ngroup * nprocs, numSystems, sysNum, newcomm);
          if (newcomm != MPI_COMM_NULL) {  // remove unnecessary procs
            BOOST_TEST_CHECKPOINT("worker sysNum: " + std::to_string(sysNum));
            TestParams testParams(dim, lmin, lmax, boundary, ngroup, nprocs, ncombi, sysNum,
                                  newcomm);
            if (getCommRank(MPI_COMM_WORLD) == 0) {
              std::cout << "third level worker only " << boundary << " " << ngroup << "*" << nprocs
                        << std::endl;
            }
            CHECK_NO_THROW_LOUD(testCombineThirdLevelWithoutManagers(testParams));
            MPI_Barrier(newcomm);
            BOOST_TEST_CHECKPOINT("sysNum " + std::to_string(sysNum) + " finished");
          }
          MPI_Barrier(MPI_COMM_WORLD);
        }
      }
    }
  }
}

// like test_5, but send only dummy data between managers
BOOST_AUTO_TEST_CASE(test_7, *boost::unit_test::tolerance(TestHelper::tolerance) *
                                 boost::unit_test::disabled()) {
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
      assignProcsToSystems(ngroup * nprocs + 1, numSystems, sysNum, newcomm);

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
BOOST_AUTO_TEST_CASE(test_8, *boost::unit_test::tolerance(TestHelper::tolerance) *
                                 boost::unit_test::disabled()) {
  unsigned int numSystems = 2;
  unsigned int nprocs = 1;
  unsigned int ncombi = 10;
  DimType dim = 6;
  BoundaryType boundary = 1;
  // LevelVector lmin(dim, 4);
  // LevelVector lmax(dim, 7);
  LevelVector lmin(dim, 2);
  LevelVector lmax(dim, 5);

  for (unsigned int ngroup : std::vector<unsigned int>({2, 3})) {
    unsigned int sysNum;
    CommunicatorType newcomm = MPI_COMM_NULL;
    assignProcsToSystems(ngroup * nprocs + 1, numSystems, sysNum, newcomm);

    if (newcomm != MPI_COMM_NULL) {  // remove unnecessary procs
      TestParams testParams(dim, lmin, lmax, boundary, ngroup, nprocs, ncombi, sysNum, newcomm);
      startInfrastructure();
      BOOST_CHECK_NO_THROW(testCombineThirdLevel(testParams, false)); //TODO make true again
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }
}

BOOST_AUTO_TEST_CASE(test_workers_2d, *boost::unit_test::tolerance(TestHelper::tolerance) *
                                          boost::unit_test::timeout(100)) {
  unsigned int numSystems = 2;
  unsigned int nprocs = 1;
  unsigned int ncombi = 5;
  DimType dim = 2;
  BoundaryType boundary = 1;
  LevelVector lmin(dim, 2);
  LevelVector lmax(dim, 8);

  for (unsigned int ngroup : std::vector<unsigned int>({1, 2, 3, 4})) {
    unsigned int sysNum;
    CommunicatorType newcomm = MPI_COMM_NULL;
    assignProcsToSystems(ngroup * nprocs, numSystems, sysNum, newcomm);
    for (CombinationVariant variant :
         {CombinationVariant::sparseGridReduce, CombinationVariant::outgroupSparseGridReduce,
          CombinationVariant::chunkedOutgroupSparseGridReduce}) {
      if (newcomm != MPI_COMM_NULL) {  // remove unnecessary procs
        TestParams testParams(dim, lmin, lmax, boundary, ngroup, nprocs, ncombi, sysNum, newcomm);
        BOOST_TEST_MESSAGE("test_workers_small: " + std::to_string(variant));
        CHECK_NO_THROW_LOUD(testCombineThirdLevelWithoutManagers(testParams, variant));
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
}

// same as test_8 but only with workers
BOOST_AUTO_TEST_CASE(test_8_workers, *boost::unit_test::tolerance(TestHelper::tolerance) *
                                         boost::unit_test::timeout(1250)) {
  unsigned int numSystems = 2;
  unsigned int nprocs = 1;
  unsigned int ncombi = 3;
  DimType dim = 6;
  BoundaryType boundary = 1;
  LevelVector lmin(dim, 2);
#ifdef NDEBUG
  LevelVector lmax(dim, 9);
#else
  LevelVector lmax(dim, 5);
#endif

  for (unsigned int ngroup : std::vector<unsigned int>({2, 3, 4})) {
    BOOST_TEST_MESSAGE("num groups: " + std::to_string(ngroup));
    unsigned int sysNum;
    CommunicatorType newcomm = MPI_COMM_NULL;
    assignProcsToSystems(ngroup * nprocs, numSystems, sysNum, newcomm);
    for (CombinationVariant variant :
         {CombinationVariant::sparseGridReduce, CombinationVariant::outgroupSparseGridReduce,
          CombinationVariant::chunkedOutgroupSparseGridReduce}) {
      if (newcomm != MPI_COMM_NULL) {  // remove unnecessary procs
        TestParams testParams(dim, lmin, lmax, boundary, ngroup, nprocs, ncombi, sysNum, newcomm);
        BOOST_TEST_MESSAGE("test_8_workers: " + std::to_string(variant));
        CHECK_NO_THROW_LOUD(testCombineThirdLevelWithoutManagers(testParams, variant));
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
}

BOOST_AUTO_TEST_CASE(test_workers_three_systems_2d,
                     *boost::unit_test::tolerance(TestHelper::tolerance) *
                         boost::unit_test::timeout(250)) {
  unsigned int numSystems = 3;
  unsigned int ncombi = 3;
  DimType dim = 2;
  BoundaryType boundary = 1;
  LevelVector lmin(dim, 2);
  LevelVector lmax(dim, 9);

  for (unsigned int ngroup : std::vector<unsigned int>({1, 2})) {
    unsigned int nprocs = 2 / ngroup;
    unsigned int sysNum;
    CommunicatorType newcomm = MPI_COMM_NULL;
    assignProcsToSystems(ngroup * nprocs, numSystems, sysNum, newcomm);
    auto variant = CombinationVariant::chunkedOutgroupSparseGridReduce;
    if (newcomm != MPI_COMM_NULL) {  // remove unnecessary procs
      TestParams testParams(dim, lmin, lmax, boundary, ngroup, nprocs, ncombi, sysNum, newcomm,
                            numSystems);
      BOOST_TEST_MESSAGE("test_workers_three_systems: " + std::to_string(variant));
      CHECK_NO_THROW_LOUD(testCombineThirdLevelWithoutManagers(testParams, variant));
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

BOOST_AUTO_TEST_CASE(test_workers_three_systems_6d,
                     *boost::unit_test::tolerance(TestHelper::tolerance) *
                         boost::unit_test::timeout(550)) {
  unsigned int numSystems = 3;
  unsigned int ncombi = 3;
  DimType dim = 6;
  BoundaryType boundary = 1;
  LevelVector lmin(dim, 2);
  LevelVector lmax(dim, 5);

  for (unsigned int ngroup : std::vector<unsigned int>({1, 2})) {
    unsigned int nprocs = 2 / ngroup;
    unsigned int sysNum;
    CommunicatorType newcomm = MPI_COMM_NULL;
    assignProcsToSystems(ngroup * nprocs, numSystems, sysNum, newcomm);
    auto variant = CombinationVariant::chunkedOutgroupSparseGridReduce;
    if (newcomm != MPI_COMM_NULL) {  // remove unnecessary procs
      TestParams testParams(dim, lmin, lmax, boundary, ngroup, nprocs, ncombi, sysNum, newcomm,
                            numSystems);
      BOOST_TEST_MESSAGE("test_workers_three_systems: " + std::to_string(variant));
      CHECK_NO_THROW_LOUD(testCombineThirdLevelWithoutManagers(testParams, variant));
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

BOOST_AUTO_TEST_SUITE_END()
