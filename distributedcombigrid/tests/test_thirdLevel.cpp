#define BOOST_TEST_DYN_LINK
#include <mpi.h>
#include <boost/serialization/export.hpp>
#include <boost/test/unit_test.hpp>
#include "TaskConstParaboloid.hpp"
#include "TaskConst.hpp"
#include "test_helper.hpp"
#include "stdlib.h"

#include "sgpp/distributedcombigrid/sparsegrid/DistributedSparseGridUniform.hpp"

#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"
#include "sgpp/distributedcombigrid/combischeme/CombiThirdLevelScheme.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LearningLoadModel.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LinearLoadModel.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupManager.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupWorker.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessManager.hpp"
#include "sgpp/distributedcombigrid/task/Task.hpp"
#include "sgpp/distributedcombigrid/utils/Config.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"

BOOST_CLASS_EXPORT(TaskConstParaboloid)

class TestParams{
  public:
    DimType dim = 2;
    LevelVector lmin;
    LevelVector lmax;
    unsigned int ngroup = 1;
    unsigned int nprocs = 1;
    unsigned int ncombi = 1;
    unsigned int sysNum = 0;
    const CommunicatorType& comm;
    std::string host = "localhost";
    unsigned short port = 9999;

    TestParams(DimType dim, LevelVector& lmin, LevelVector& lmax,
               unsigned int ngroup, unsigned int nprocs, unsigned int ncombi,
               unsigned int sysNum, const CommunicatorType& comm,
               const std::string& host = "localhost",
               unsigned short dataPort = 9999)
               : dim(dim), lmin(lmin), lmax(lmax), ngroup(ngroup),
                 nprocs(nprocs), ncombi(ncombi), sysNum(sysNum), comm(comm),
                 host(host), port(dataPort)
  {
  }
};

/**
* Checks if combination was successful.
* Since the tasks don't evolve over time the expected result should match the
* initial function values.
*/
bool checkReducedFullGrid(ProcessGroupWorker& worker) {
  TaskContainer& tasks = worker.getTasks();
  int numGrids = (int) worker.getCombiParameters().getNumGrids();

  for (Task* t : tasks) {
    for (int g = 0; g < numGrids; g++) {
      DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid(g);
      ParaboloidFn<CombiDataType> initialFunction;
      for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
        std::vector<double> coords(dfg.getDimension());
        dfg.getCoordsLocal(li, coords);
        CombiDataType expected = initialFunction(coords);
        CombiDataType occuring = dfg.getData()[li];
        double diff = abs(occuring - expected);
        if (diff > TestHelper::tolerance) {
          std::cout << "Found value that does not match expected at index "
                    << li  << ": "  << occuring << " != "
                    << expected << ", diff = " << diff << std::endl;
          return false;
        }
      }
    }
  }
  return true;
}

void assignProcsToSystems(unsigned int ngroup, unsigned int nprocs,
                          unsigned int numSystems, unsigned int& sysNum,
                          CommunicatorType& newcomm) {
  int rank, size, color, key;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  unsigned int procsPerSys = ngroup * nprocs + 1;

  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(int(numSystems * procsPerSys)));

  // assign procs to systems
  sysNum = unsigned(rank) / procsPerSys;
  color = int(sysNum);
  key = rank % (int) procsPerSys;

  MPI_Comm_split(MPI_COMM_WORLD, color, key, &newcomm);
}

/** Runs the third level manager in the background as a forked child process */
void runThirdLevelManager() {
  std::cout << "starting thirdLevelManager..." << std::endl;
  std::string command = "../../distributedcombigrid/third_level_manager/run.sh &";
  system(command.c_str());
}

/** Runs the tl manager*/
void startInfrastructure() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    runThirdLevelManager();
  }
  // give infrastructure some time to set up
  sleep(2);
}

void testCombineThirdLevel(TestParams& testParams) {

  BOOST_CHECK(testParams.comm != MPI_COMM_NULL);

  size_t procsPerSys = testParams.ngroup * testParams.nprocs + 1;

  combigrid::Stats::initialize();

  theMPISystem()->initWorldReusable(testParams.comm, testParams.ngroup,
                                    testParams.nprocs);

  WORLD_MANAGER_EXCLUSIVE_SECTION {

    ProcessGroupManagerContainer pgroups;
    for (size_t i = 0; i < testParams.ngroup; ++i) {
      int pgroupRootID((int) i);
      pgroups.emplace_back(std::make_shared<ProcessGroupManager>(pgroupRootID));
    }

    auto loadmodel = std::unique_ptr<LoadModel>(new LinearLoadModel());
    std::vector<bool> boundary(testParams.dim, false);

    // create third level specific scheme
    CombiMinMaxScheme combischeme(testParams.dim, testParams.lmin, testParams.lmax);
    combischeme.createClassicalCombischeme();
    //combischeme.createAdaptiveCombischeme();
    // get full scheme first
    std::vector<LevelVector> fullLevels = combischeme.getCombiSpaces();
    std::vector<combigrid::real> fullCoeffs = combischeme.getCoeffs();
    // split scheme and assign each half to a system
    std::vector<LevelVector> levels;
    std::vector<combigrid::real> coeffs;
    CombiThirdLevelScheme::createThirdLevelScheme(fullLevels, fullCoeffs,
                                                  boundary, testParams.sysNum,
                                                  2, levels, coeffs);

    //std::cout << "Combischeme " << systemName << ":" << std::endl;
    //for (const auto& l : levels)
    //  std::cout << toString(l) << std::endl;


    BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable((int)procsPerSys));

    // create Tasks
    TaskContainer tasks;
    std::vector<int> taskIDs;
    for (size_t i = 0; i < levels.size(); i++) {
      Task* t = new TaskConstParaboloid(levels[i], boundary, coeffs[i], loadmodel.get());

      tasks.push_back(t);
      taskIDs.push_back(t->getID());
    }

    // create combiparameters
    IndexVector parallelization = {static_cast<long>(testParams.nprocs), 1};
    CombiParameters combiParams(testParams.dim, testParams.lmin, testParams.lmax,
                                boundary, levels, coeffs, taskIDs, testParams.ncombi,
                                1, parallelization, std::vector<IndexType>(0),
                                std::vector<IndexType>(0), testParams.host,
                                testParams.port, 0);

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

        Stats::startEvent("manager init dsgus");
        manager.initDsgus();
        Stats::stopEvent("manager init dsgus");

        // exchange subspace sizes to unify the dsgs with the remote system
        Stats::startEvent("manager unify subspace sizes with remote");
        manager.unifySubspaceSizesThirdLevel(),
        Stats::startEvent("manager unify subspace sizes with remote");
      } else {
        Stats::startEvent("manager run");
        manager.runnext();
        Stats::stopEvent("manager run");
      }
      //combine grids
      Stats::startEvent("manager combine third level");
      manager.combineThirdLevel();
      Stats::stopEvent("manager combine third level");
    }
    manager.exit();
  }
  else {
    ProcessGroupWorker pgroup;
    SignalType signal = -1;
    signal = pgroup.wait();
    while (signal != EXIT) {
      signal = pgroup.wait();
      std::cout << "Worker with rank " << theMPISystem()->getLocalRank() << " processed signal " << signal << std::endl; 
    }

    // after combination check workers grids
    BOOST_CHECK(checkReducedFullGrid(pgroup));
  }

  combigrid::Stats::finalize();
  MPI_Barrier(testParams.comm);
}

BOOST_AUTO_TEST_SUITE(thirdLevel)

/*
BOOST_AUTO_TEST_CASE(test_1, *boost::unit_test::tolerance(TestHelper::tolerance)) {
  unsigned int numSystems = 2;
  unsigned int ngroup     = 1;
  unsigned int nprocs     = 1;
  unsigned int ncombi     = 10;
  DimType dim   = 2;
  LevelVector lmin(dim, 4);
  LevelVector lmax(dim, 7);

  unsigned int sysNum;
  CommunicatorType newcomm;
  assignProcsToSystems(ngroup, nprocs, numSystems, sysNum, newcomm);
  TestParams testParams(dim, lmin, lmax, ngroup, nprocs, ncombi, sysNum, newcomm);

  startInfrastructure();

  testCombineThirdLevel(testParams, newcomm, sysNum);

  MPI_Barrier(MPI_COMM_WORLD);
}

BOOST_AUTO_TEST_CASE(test_2, *boost::unit_test::tolerance(TestHelper::tolerance)) {
  unsigned int numSystems = 2;
  unsigned int ngroup     = 1;
  unsigned int nprocs     = 4;
  unsigned int ncombi     = 10;
  LevelVector lmin(dim, 4);
  LevelVector lmax(dim, 7);
  DimType dim   = 2;

  unsigned int sysNum;
  CommunicatorType newcomm;
  assignProcsToSystems(ngroup, nprocs, numSystems, sysNum, newcomm);
  TestParams testParams(dim, lmin, lmax, ngroup, nprocs, ncombi, sysNum, newcomm);

  startInfrastructure();

  testCombineThirdLevel(testParams, newcomm, sysNum);

  MPI_Barrier(MPI_COMM_WORLD);
}*/

BOOST_AUTO_TEST_CASE(test_3, *boost::unit_test::tolerance(TestHelper::tolerance)) {
  unsigned int numSystems = 2;
  unsigned int ngroup     = 2;
  unsigned int nprocs     = 2;
  unsigned int ncombi     = 10;
  DimType dim             = 2;
  LevelVector lmin(dim, 4);
  LevelVector lmax(dim, 7);

  unsigned int sysNum;
  CommunicatorType newcomm;
  assignProcsToSystems(ngroup, nprocs, numSystems, sysNum, newcomm);

  if (sysNum < numSystems) { // remove unnecessary procs
    TestParams testParams(dim, lmin, lmax, ngroup, nprocs, ncombi, sysNum, newcomm);

    startInfrastructure();

    testCombineThirdLevel(testParams);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}


BOOST_AUTO_TEST_SUITE_END()
