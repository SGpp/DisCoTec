#define BOOST_TEST_DYN_LINK
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <complex>
#include <cstdarg>
#include <iostream>
#include <vector>

#include <boost/serialization/export.hpp>
#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"
#include "sgpp/distributedcombigrid/fault_tolerance/FaultCriterion.hpp"
#include "sgpp/distributedcombigrid/fault_tolerance/StaticFaults.hpp"
#include "sgpp/distributedcombigrid/fault_tolerance/WeibullFaults.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LinearLoadModel.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupManager.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupWorker.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessManager.hpp"
#include "sgpp/distributedcombigrid/task/Task.hpp"
#include "sgpp/distributedcombigrid/utils/Config.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "test_helper.hpp"
#include "TaskConst.hpp"

using namespace combigrid;

BOOST_CLASS_EXPORT(TaskConst)

void checkCombine(size_t ngroup = 1, size_t nprocs = 1) {
  size_t size = ngroup * nprocs + 1;
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(static_cast<int>(size)));

  CommunicatorType comm = TestHelper::getComm(size);
  if (comm == MPI_COMM_NULL) {
    return;
  }

  combigrid::Stats::initialize();

  theMPISystem()->initWorldReusable(comm, ngroup, nprocs);
  // theMPISystem()->init(ngroup, nprocs);

  WORLD_MANAGER_EXCLUSIVE_SECTION {
    ProcessGroupManagerContainer pgroups;
    for (int i = 0; i < static_cast<int>(ngroup); ++i) {
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
    // std::cout << "levels " << std::endl;
    // for (const auto& level : levels) {
    //   std::cout << toString(level) << std::endl;
    // }

    std::vector<combigrid::real> coeffs = combischeme.getCoeffs();
    // std::cout << "coeffs " << std::endl;
    // for (const auto& coeff : coeffs) {
    //   std::cout << coeff << std::endl;
    // }

    BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(static_cast<int>(size)));

    // create Tasks
    TaskContainer tasks;
    std::vector<int> taskIDs;
    for (size_t i = 0; i < levels.size(); i++) {
      Task* t = new TaskConst(levels[i], boundary, coeffs[i], loadmodel.get());
      tasks.push_back(t);
      taskIDs.push_back(t->getID());
    }

    IndexVector parallelization = {static_cast<IndexType>(nprocs), 1};
    // create combiparameters
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs, ncombi);
    params.setParallelization({static_cast<long int>(nprocs), 1}); //TODO why??

    // create abstraction for Manager
    ProcessManager manager(pgroups, tasks, params, std::move(loadmodel));

    // the combiparameters are sent to all process groups before the
    // computations start
    manager.updateCombiParameters(); //TODO move to manager constructor or runfirst?

    /* distribute task according to load model and start computation for
     * the first time */
    std::cout << "run first " << std::endl;
    manager.runfirst();

    for (size_t it = 0; it < ncombi; ++it) {
      std::cout << "combine " << std::endl;
      manager.combine();
    }

    // evaluate solution
    FullGrid<CombiDataType> fg_eval(dim, leval, boundary);
    manager.gridEval(fg_eval);

    // compare with known results:
    // point in the middle
    CombiDataType midResult = fg_eval.getData()[fg_eval.getNrElements() / 2];
    std::cout << "midResult " << fabs(midResult) << std::endl;
    BOOST_TEST(fabs(midResult) == 1.333333333);

    manager.exit();
  }
  else {
    ProcessGroupWorker pgroup;
    SignalType signal = -1;
    while (signal != EXIT) signal = pgroup.wait();
  }

  combigrid::Stats::finalize();
  MPI_Barrier(comm);
  TestHelper::testStrayMessages(comm);
}

BOOST_AUTO_TEST_SUITE(reduce)

BOOST_AUTO_TEST_CASE(test_1, *boost::unit_test::tolerance(TestHelper::higherTolerance) *
                                 boost::unit_test::timeout(60)) {
  std::cout << "reduce/test_1"<< std::endl;
  checkCombine(1,1);
}

BOOST_AUTO_TEST_CASE(test_2, *boost::unit_test::tolerance(TestHelper::higherTolerance) *
                                 boost::unit_test::timeout(60)) {
  std::cout << "reduce/test_2"<< std::endl;
  checkCombine(1,2);
}

BOOST_AUTO_TEST_CASE(test_3, *boost::unit_test::tolerance(TestHelper::higherTolerance) *
                                 boost::unit_test::timeout(60)) {
  std::cout << "reduce/test_3"<< std::endl;
  checkCombine(2,2);
}

BOOST_AUTO_TEST_CASE(test_4, *boost::unit_test::tolerance(TestHelper::higherTolerance) *
                                 boost::unit_test::timeout(60)) {
  std::cout << "reduce/test_4"<< std::endl;
  checkCombine(2,4);
}

BOOST_AUTO_TEST_SUITE_END()
