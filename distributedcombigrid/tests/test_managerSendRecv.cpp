#define BOOST_TEST_DYN_LINK
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/serialization/export.hpp>
#include "test_helper.hpp"
#include "TaskConst.hpp"
// BOOST_CLASS_EXPORT_IMPLEMENT(TaskConst)

#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LinearLoadModel.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LearningLoadModel.hpp"
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
  //  notify process of process group of completion, upon which they call the second combine/allreduce
  // }


  //technology considerations If receiver as extra thread: MPI+threads or MPI+MPI(with RMA)?
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
  //  notify process of process group of completion, upon which they call the second combine/allreduce
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
    std::vector<bool> boundary(dim, true);

    CombiMinMaxScheme combischeme(dim, lmin, lmax);
    // combischeme.createClassicalCombischeme();
    combischeme.createAdaptiveCombischeme();
    std::vector<LevelVector> levels = combischeme.getCombiSpaces();
    std::vector<combigrid::real> coeffs = combischeme.getCoeffs();
    std::cout << "levels " << std::endl;
    for (const auto& level : levels) {
      std::cout << toString(level) << std::endl;
    }

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
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs, ncombi, 1, parallelization);

    // create abstraction for Manager
    ProcessManager manager(pgroups, tasks, params, std::move(loadmodel));

    manager.runfirst();
    manager.combine();

    manager.exit();

    std::cout << "manager exit" << std::endl;
    //TODO introduce signal for grid gathering
    int numGrids = params.getNumGrids();  
    std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>> combinedUniDSGVector(numGrids);
    // iterate first process group
    for(size_t i = 0; i < nprocs; ++i){
      for (int g = 0; g < numGrids; ++g) {
        combinedUniDSGVector[g].reset(recvDSGUniform<CombiDataType>(i, theMPISystem()->getWorldComm()));
        const auto & dsgu = combinedUniDSGVector[g];
        BOOST_CHECK(dsgu->getDim() > 0);
        BOOST_CHECK(!dsgu->getBoundaryVector().empty());
        BOOST_CHECK(dsgu->getNMax()[0] >= 0);
        BOOST_CHECK(dsgu->getNumSubspaces() > 0);
        BOOST_CHECK(!dsgu->getDataVector(0).empty());
        BOOST_CHECK(dsgu->getDataVector(dsgu->getNumSubspaces() - 1).size() >= 0);

        bool found = false;
        for (size_t j = 0; j < dsgu->getNumSubspaces(); ++j){
          std::cerr << combigrid::toString(dsgu->getDataVector(j)) << std::endl;
          auto it = std::find (dsgu->getDataVector(j).begin(), dsgu->getDataVector(j).end(), CombiDataType(4.0/3.0,0) );
          if ( it != dsgu->getDataVector(j).end() ) {
            found = true;
          }
        }
        BOOST_TEST(found);
      }
    }

  }
  else{
    ProcessGroupWorker pgroup;
    SignalType signal = -1;
    while (signal != EXIT) signal = pgroup.wait();

    std::cout << "worker exit" << std::endl;
    pgroup.initCombinedUniDSGVector();
    pgroup.hierarchizeUniformSG();
    pgroup.reduceUniformSG();
    pgroup.sendSparseGridToManager();
  }

  combigrid::Stats::finalize();
  MPI_Barrier(comm);
}

void checkAddSparseGridToProcessGroup(){
  // put subspace data into buffer for allreduce //TODO
  // iterate process group
}

BOOST_AUTO_TEST_SUITE(managerSendRecv)

BOOST_AUTO_TEST_CASE(test_1, * boost::unit_test::tolerance(TestHelper::tolerance) * boost::unit_test::timeout(40)) {
  // use recombination
  // checkManager(true, false, 1.54369, 11.28857);
  checkGatherSparseGridFromProcessGroup();
}

BOOST_AUTO_TEST_CASE(test_2, * boost::unit_test::tolerance(TestHelper::tolerance) * boost::unit_test::timeout(40)) {
  // use recombination
  // checkManager(true, false, 1.54369, 11.28857);
  checkGatherSparseGridFromProcessGroup(2,2);
}

BOOST_AUTO_TEST_SUITE_END()