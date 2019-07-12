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
  std::string command = "../examples/gene_distributed_third_level/third_level_manager/run.sh &";
  system(command.c_str());

  // give the thirdLevelManger some time to set up
  sleep(1);
}

/** Kills the tl manager and rabbit mq server */
void killInfrastructure() {
  //system("killall rabbitmq-server");
  //system("killall thirdLevelManager");
}

/** Runs the tl manager and RabbitMQ server*/
void startInfrastructure() {
  //runRabbitMQServer();
  //runThirdLevelManager();
}

void testCombineThirdLevel(size_t ngroup = 1, size_t nprocs = 1, const std::string& sysName = "system1", CommunicatorType comm = MPI_COMM_NULL) {
  size_t size = ngroup * nprocs + 1;

  if (comm == MPI_COMM_NULL) {
    return;
  }

  combigrid::Stats::initialize();

  theMPISystem()->initWorldReusable(comm, ngroup, nprocs);

  WORLD_MANAGER_EXCLUSIVE_SECTION {
    int argc = boost::unit_test::framework::master_test_suite().argc;
    char** argv = boost::unit_test::framework::master_test_suite().argv;

    ProcessGroupManagerContainer pgroups;
    for (int i = 0; i < ngroup; ++i) {
      int pgroupRootID(i);
      pgroups.emplace_back(std::make_shared<ProcessGroupManager>(pgroupRootID));
    }

    auto loadmodel = std::unique_ptr<LoadModel>(new LinearLoadModel());

    DimType dim = 2;
    LevelVector lmin(dim, 4);
    LevelVector lmax(dim, 7);

    size_t ncombi = 1;
    std::vector<bool> boundary(dim, false);

    // create third level specific scheme
    CombiMinMaxScheme combischeme(dim, lmin, lmax);
    combischeme.createClassicalCombischeme();
    //combischeme.createAdaptiveCombischeme();


    std::vector<LevelVector> levels = combischeme.getCombiSpaces();
    std::vector<combigrid::real> coeffs = combischeme.getCoeffs();

    std::vector<LevelVector> commonSubspaces;

    unsigned int sysNum;
    if (sysName == "system1") {
      sysNum = 0;
      CombiThirdLevelScheme::createThirdLevelScheme(levels, coeffs, commonSubspaces, boundary, sysNum, 2);
      std::cout << "Common Subspace" << std::endl;
      for (const auto& ss : commonSubspaces)
        std::cout << toString(ss) << std::endl;
    } else {
      sysNum = 1;
      CombiThirdLevelScheme::createThirdLevelScheme(levels, coeffs, commonSubspaces, boundary, sysNum, 2);
    }

    //std::cout << "Combischeme " << sysName << ":" << std::endl;
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

    IndexVector parallelization = {static_cast<long>(nprocs), 1};
    // create combiparameters
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs, ncombi, 1,
                           parallelization, std::vector<IndexType>(0),
                           std::vector<IndexType>(0), "localhost", 9999, sysName, commonSubspaces, 0);

    // create abstraction for Manager
    ProcessManager manager(pgroups, tasks, params, std::move(loadmodel));

    // the combiparameters are sent to all process groups before the
    // computations start
    manager.updateCombiParameters();

    std::cout << "running first";
    manager.runfirst();
    std::cout << "running combineThirdLevel";
    manager.combineThirdLevel();
    manager.exit();

    if (argc == 1) {
      killInfrastructure();
    }
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

BOOST_AUTO_TEST_CASE(test_1, *boost::unit_test::tolerance(TestHelper::tolerance)) {

  int argc = boost::unit_test::framework::master_test_suite().argc;
  char** const & argv =  boost::unit_test::framework::master_test_suite().argv;

  if (argc != 3) {
    std::cout << "Please specify number of groups and procs per group" << std::endl;
    return;
  }

  size_t ngroup = static_cast<size_t>(std::stoi(argv[1]));
  size_t nprocs = static_cast<size_t>(std::stoi(argv[2]));

  int rank, size, color;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  size_t procsPerSys = ngroup * nprocs + 1;

  if (size_t(size) < 2*procsPerSys) {
    std::cout << "Failed running test: Number of Processes given: " << size << " required: " 
              << 2*procsPerSys << std::endl;
    return;
  }

  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable( 2*procsPerSys ));

  if (rank == 0)
    startInfrastructure();

  // assign procs to systems
  std::string systemName;
  if (size_t(rank) < procsPerSys) {
    color = 0;
    systemName = "system1";
  } else if (size_t(rank) < 2 * procsPerSys) {
    color = 1;
    systemName = "system2";
  } else {
    color = 2;
  }

  MPI_Comm newcomm;
  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &newcomm);

  // remove unnecessary procs
  if (color == 2)
    return;

  testCombineThirdLevel(ngroup, nprocs, systemName, newcomm);

  if (rank == 0)
    killInfrastructure();
}

BOOST_AUTO_TEST_SUITE_END()
