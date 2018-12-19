#define BOOST_TEST_DYN_LINK
#include <mpi.h>
#include <boost/serialization/export.hpp>
#include <boost/test/unit_test.hpp>
#include "TaskConst.hpp"
#include "test_helper.hpp"
// BOOST_CLASS_EXPORT_IMPLEMENT(TaskConst)

#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LearningLoadModel.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LinearLoadModel.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupManager.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupWorker.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessManager.hpp"
#include "sgpp/distributedcombigrid/task/Task.hpp"
#include "sgpp/distributedcombigrid/utils/Config.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"

void checkManagerRecv() {
  // Recv algorithm:
  // get handle to process group to write to
  // set up receive channel to middleman
  // after every first combine from manager{
  //  listen & receive
  //  write into extra process group's recv-memory, one after the other process
  //  //  notify manager if completed
  //  notify process of process group of completion, upon which they call the second
  //  combine/allreduce
  // }

  // technology considerations If receiver as extra thread: MPI+threads or MPI+MPI(with RMA)?
  // cf wgropp.cs.illinois.edu/courses/cs598-s16/lectures/lecture36.pdf
  // => rather advises towards RMA
}

void checkManagerSend() {
  // Send algorithm:
  // get handle to process group to read from
  // set up send channel to middleman
  // after every first combine from manager{
  //  gather grid
  //  stream it to middleman/other machines
  //  notify process of process group of completion, upon which they call the second
  //  combine/allreduce
  // }
}

void checkGatherSparseGridFromProcessGroup(size_t ngroup = 1, size_t nprocs = 1) {
  size_t size = ngroup * nprocs + 1;
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(size));

  CommunicatorType comm = TestHelper::getComm(size);
  if (comm == MPI_COMM_NULL) {
    return;
  }

  combigrid::Stats::initialize();

  theMPISystem()->initWorldReusable(comm, ngroup, nprocs);
  // theMPISystem()->init(ngroup, nprocs);

  // set up a constant valued distributed fullgrid in the process groups
  WORLD_MANAGER_EXCLUSIVE_SECTION {
    ProcessGroupManagerContainer pgroups;
    for (int i = 0; i < ngroup; ++i) {
      int pgroupRootID(i);
      pgroups.emplace_back(std::make_shared<ProcessGroupManager>(pgroupRootID));
    }

    auto loadmodel = std::unique_ptr<LoadModel>(new LinearLoadModel());

    DimType dim = 2;
    LevelVector lmin(dim, 2);
    LevelVector lmax(dim, 4), leval(dim, 4);

    size_t ncombi = 2;
    std::vector<bool> boundary(dim, false);

    CombiMinMaxScheme combischeme(dim, lmin, lmax);
    // combischeme.createClassicalCombischeme();
    combischeme.createAdaptiveCombischeme();
    std::vector<LevelVector> levels = combischeme.getCombiSpaces();
    std::vector<combigrid::real> coeffs = combischeme.getCoeffs();

    BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(size));

    // create Tasks
    TaskContainer tasks;
    std::vector<int> taskIDs;
    for (size_t i = 0; i < levels.size(); i++) {
      Task* t = new TaskConst(levels[i], boundary, coeffs[i], loadmodel.get());
      tasks.push_back(t);
      taskIDs.push_back(t->getID());
    }

    IndexVector parallelization = {nprocs, 1};
    // create combiparameters
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs, ncombi, 1,
                           parallelization);

    // create abstraction for Manager
    ProcessManager manager(pgroups, tasks, params, std::move(loadmodel));

    manager.runfirst();
    manager.combine();
    std::cerr << "manager get DSG" << std::endl;
    manager.getDSGFromProcessGroup();
    manager.exit();

    int numGrids = params.getNumGrids();
    auto& combinedUniDSGVector = manager.getOutboundUniDSGVector();
    // iterate first process group
    bool found = false;
    for (size_t i = 0; i < nprocs; ++i) {
      for (int g = 0; g < numGrids; ++g) {
        const auto& dsgu = combinedUniDSGVector[g];
        BOOST_CHECK(dsgu->getDim() > 0);
        BOOST_CHECK(!dsgu->getBoundaryVector().empty());
        BOOST_CHECK(dsgu->getNMax()[0] >= 0);
        BOOST_CHECK(dsgu->getNumSubspaces() > 0);

        LevelVector smallestSubspace = {1, 1};
        std::cerr << combigrid::toString(dsgu->getDataVector(smallestSubspace)) << std::endl;
        auto it = std::find_if(dsgu->getDataVector(smallestSubspace).begin(),
                               dsgu->getDataVector(smallestSubspace).end(), [](CombiDataType& dve) {
                                 return (abs(dve) - 1.333333333333333) < TestHelper::tolerance;
                               });
        if (it != dsgu->getDataVector(smallestSubspace).end()) {
          found = true;
        }
      }
    }
    BOOST_TEST(found);
  }
  else {
    ProcessGroupWorker pgroup;
    SignalType signal = -1;
    signal = pgroup.wait();
    while (signal != EXIT) {
      signal = pgroup.wait();
      // for introspection:
      // if (signal == RUN_FIRST || signal == RUN_NEXT) {
      //   pgroup.initCombinedUniDSGVector();
      //   pgroup.hierarchizeFullGrids();
      //   pgroup.addFullGridsToUniformSG();
      //   pgroup.reduceUniformSG();
      // for (auto& dsg : pgroup.combinedUniDSGVector_) {
      //   for (size_t j = 0; j < dsg->getNumSubspaces(); ++j) {
      //     std::vector<CombiDataType>& subspaceData = dsg->getDataVector(j);
      //     std::cerr << combigrid::toString(subspaceData) << std::endl;
      //   }
      // }
      //   pgroup.dehierarchizeFullGrids();
      // }
    }
  }

  combigrid::Stats::finalize();
  MPI_Barrier(comm);
}

void checkAddSparseGridToProcessGroup() {
  // put subspace data into buffer for allreduce //TODO
  // iterate process group
}

BOOST_AUTO_TEST_SUITE(managerSendRecv)

BOOST_AUTO_TEST_CASE(test_1, *boost::unit_test::tolerance(TestHelper::tolerance) *
                                 boost::unit_test::timeout(40)) {
  checkGatherSparseGridFromProcessGroup(1, 1);
}

BOOST_AUTO_TEST_CASE(test_2, *boost::unit_test::tolerance(TestHelper::tolerance) *
                                 boost::unit_test::timeout(40)) {
  checkGatherSparseGridFromProcessGroup(1, 2);
}

BOOST_AUTO_TEST_CASE(test_3, *boost::unit_test::tolerance(TestHelper::tolerance) *
                                 boost::unit_test::timeout(40)) {
  checkGatherSparseGridFromProcessGroup(2, 2);
}

BOOST_AUTO_TEST_SUITE_END()