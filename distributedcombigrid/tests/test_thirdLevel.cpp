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
    size_t ngroup = 1;
    size_t nprocs = 1;
    size_t ncombi = 1;
    std::string systemName = "";
    std::string host = "localhost";
    unsigned short dataPort = 9999;

    TestParams(DimType dim, LevelVector& lmin, LevelVector& lmax, size_t ngroup,
               size_t nprocs, size_t ncombi, const std::string& systemName,
               const std::string& host = "localhost", unsigned short dataPort = 9999)
               : dim(dim), lmin(lmin), lmax(lmax), ngroup(ngroup), nprocs(nprocs),
                 ncombi(ncombi), host(host), systemName(systemName),
                 dataPort(dataPort) {}
};

/**
* Checks if combination was successful.
* Since the task doesn't evolve over time the expected result should match the
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

void assignProcToSystem(size_t ngroup, size_t nprocs, CommunicatorType& newcomm,
                        std::string& systemName) {
  int rank, size, color;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  size_t procsPerSys = ngroup * nprocs + 1;

  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable( 2*procsPerSys ));

  // assign procs to systems
  if (size_t(rank) < procsPerSys) {
    color = 0;
    systemName = "system0";
  } else if (size_t(rank) < 2 * procsPerSys) {
    color = 1;
    systemName = "system1";
  } else {
    color = 2;
  }

  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &newcomm);

  // remove unnecessary procs
  if (color == 2)
    return;
}

/** Runs the RabbitMQ server */
void runRabbitMQServer() {
  std::cout << "starting RabbitMQ server..." << std::endl;
  std::string command = "rabbitmq-server &";
  system(command.c_str());

  // give rabbitmq some time to set up
  sleep(10);
}

/** Runs the third level manager in the background as a forked child process */
void runThirdLevelManager() {
  std::cout << "starting thirdLevelManager..." << std::endl;
  std::string command = "../../distributedcombigrid/third_level_manager/run.sh &";
  system(command.c_str());

  // give the thirdLevelManger some time to set up
  sleep(1);
}

/** Kills the tl manager and rabbit mq server */
void killInfrastructure() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    //system("killall rabbitmq-server");
    //system("killall thirdLevelManager");
  }
}

/** Runs the tl manager and RabbitMQ server*/
void startInfrastructure() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    //runRabbitMQServer();
    //runThirdLevelManager();
  }
}

void testCombineThirdLevel(TestParams& testParams, CommunicatorType& comm) {

  BOOST_CHECK(comm != MPI_COMM_NULL);

  size_t size = testParams.ngroup * testParams.nprocs + 1;

  combigrid::Stats::initialize();

  theMPISystem()->initWorldReusable(comm, testParams.ngroup, testParams.nprocs);

  WORLD_MANAGER_EXCLUSIVE_SECTION {

    ProcessGroupManagerContainer pgroups;
    for (int i = 0; i < testParams.ngroup; ++i) {
      int pgroupRootID(i);
      pgroups.emplace_back(std::make_shared<ProcessGroupManager>(pgroupRootID));
    }

    auto loadmodel = std::unique_ptr<LoadModel>(new LinearLoadModel());
    std::vector<bool> boundary(testParams.dim, false);

    // create third level specific scheme
    CombiMinMaxScheme combischeme(testParams.dim, testParams.lmin, testParams.lmax);
    combischeme.createClassicalCombischeme();
    //combischeme.createAdaptiveCombischeme();


    std::vector<LevelVector> levels = combischeme.getCombiSpaces();
    std::vector<combigrid::real> coeffs = combischeme.getCoeffs();

    std::vector<LevelVector> commonSubspaces;

    unsigned int sysNum;
    if (testParams.systemName == "system0") {
      sysNum = 0;
      CombiThirdLevelScheme::createThirdLevelScheme(levels, coeffs, commonSubspaces, boundary, sysNum, 2);
      std::cout << "Common Subspace" << std::endl;
      for (const auto& ss : commonSubspaces)
        std::cout << toString(ss) << std::endl;
    } else {
      sysNum = 1;
      CombiThirdLevelScheme::createThirdLevelScheme(levels, coeffs, commonSubspaces, boundary, sysNum, 2);
    }

    //std::cout << "Combischeme " << systemName << ":" << std::endl;
    //for (const auto& l : levels)
    //  std::cout << toString(l) << std::endl;


    BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(size));

    // create Tasks
    TaskContainer tasks;
    std::vector<int> taskIDs;
    for (size_t i = 0; i < levels.size(); i++) {
      Task* t = new TaskConstParaboloid(levels[i], boundary, coeffs[i], loadmodel.get());
      BOOST_CHECK(true);

      tasks.push_back(t);
      taskIDs.push_back(t->getID());
    }

    BOOST_CHECK(true);

    IndexVector parallelization = {static_cast<long>(testParams.nprocs), 1};
    // create combiparameters
    CombiParameters combiParams(testParams.dim, testParams.lmin, testParams.lmax,
                                boundary, levels, coeffs, taskIDs, testParams.ncombi,
                                1, parallelization, std::vector<IndexType>(0),
                                std::vector<IndexType>(0), testParams.host,
                                testParams.dataPort, testParams.systemName,
                                commonSubspaces, 0);

    // create abstraction for Manager
    ProcessManager manager(pgroups, tasks, combiParams, std::move(loadmodel));

    // the combiparameters are sent to all process groups before the
    // computations start
    manager.updateCombiParameters();

    for (int i = 0; i < testParams.ncombi; i++) {
      if (i == 0) {
        Stats::startEvent("manager run");
        manager.runfirst();
        Stats::stopEvent("manager run");
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
  MPI_Barrier(comm);
}

BOOST_AUTO_TEST_SUITE(thirdLevel)

/*
BOOST_AUTO_TEST_CASE(test_1, *boost::unit_test::tolerance(TestHelper::tolerance)) {
  DimType dim   = 2;
  size_t ngroup = 1;
  size_t nprocs = 1;
  size_t ncombi  = 10;
  LevelVector lmin(dim, 4);
  LevelVector lmax(dim, 7);
  std::string systemName = "";

  CommunicatorType newcomm;
  assignProcToSystem(ngroup, nprocs, newcomm, systemName);
  TestParams testParams(dim, lmin, lmax, ngroup, nprocs, ncombi, systemName);

  startInfrastructure();

  if (systemName != "") {
    testCombineThirdLevel(testParams, newcomm);
  }

  MPI_Barrier(MPI_COMM_WORLD);
}

BOOST_AUTO_TEST_CASE(test_2, *boost::unit_test::tolerance(TestHelper::tolerance)) {
  DimType dim   = 2;
  size_t ngroup = 1;
  size_t nprocs = 4;
  size_t ncombi  = 10;
  LevelVector lmin(dim, 4);
  LevelVector lmax(dim, 7);
  std::string systemName = "";

  CommunicatorType newcomm;
  assignProcToSystem(ngroup, nprocs, newcomm, systemName);
  TestParams testParams(dim, lmin, lmax, ngroup, nprocs, ncombi, systemName);

  startInfrastructure();

  if (systemName != "")
    testCombineThirdLevel(testParams, newcomm);

  MPI_Barrier(MPI_COMM_WORLD);
}*/

BOOST_AUTO_TEST_CASE(test_3, *boost::unit_test::tolerance(TestHelper::tolerance)) {
  DimType dim   = 2;
  size_t ngroup = 2;
  size_t nprocs = 2;
  size_t ncombi  = 10;
  LevelVector lmin(dim, 4);
  LevelVector lmax(dim, 7);
  std::string systemName = "";

  CommunicatorType newcomm;
  assignProcToSystem(ngroup, nprocs, newcomm, systemName);
  TestParams testParams(dim, lmin, lmax, ngroup, nprocs, ncombi, systemName);

  startInfrastructure();

  if (systemName != "")
    testCombineThirdLevel(testParams, newcomm);

  MPI_Barrier(MPI_COMM_WORLD);

  killInfrastructure();
}


BOOST_AUTO_TEST_SUITE_END()
