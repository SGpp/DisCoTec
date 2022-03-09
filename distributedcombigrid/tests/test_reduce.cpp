#define BOOST_TEST_DYN_LINK
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

using namespace combigrid;

/* simple task class to set all values on the grid to $levelVector_1 / levelVector_2$
 */
class TaskConst : public combigrid::Task {
 public:
  TaskConst(LevelVector& l, std::vector<bool>& boundary, real coeff, LoadModel* loadModel)
      : Task(2, l, boundary, coeff, loadModel){}

  void init(CommunicatorType lcomm, std::vector<IndexVector> decomposition) {
    // parallelization
    // assert(dfg_ == nullptr);
    long nprocs = getCommSize(lcomm);
    // IndexVector p = {nprocs, 1};
    IndexVector p = {nprocs,1};

    // decomposition = std::vector<IndexVector>(2);
    // size_t l1 = getLevelVector()[1];
    // size_t npoint_x1 = pow(2, l1) + 1;

    // decomposition[0].push_back(0);
    // decomposition[1].push_back(0);

    // // std::cout << "decomposition" << std::endl;
    // for (int r = 1; r < nprocs_; ++r) {
    //   decomposition[1].push_back(r * (npoint_x1 / nprocs_));

    //   // std::cout << decomposition[1].back() << std::endl;
    // }

    dfg_ = new DistributedFullGrid<CombiDataType>(getDim(), getLevelVector(), lcomm, getBoundary(),
                                                  p, false, decomposition);

    std::vector<CombiDataType>& elements = dfg_->getElementVector();
    for (auto& element : elements) {
      element = 10;
    }
  }

  void run(CommunicatorType lcomm) {

    std::cout << "run " << getCommRank(lcomm) << std::endl;    
    
    std::vector<CombiDataType>& elements = dfg_->getElementVector();
    for (auto& element : elements) {
      // BOOST_CHECK(abs(dfg_->getData()[li]));
      element = static_cast<double>(getLevelVector()[0]) / getLevelVector()[1];
    }
    BOOST_CHECK(dfg_);

    setFinished(true);
    
    MPI_Barrier(lcomm);
    std::cerr << "barrier" << std::endl;
  }

  void getFullGrid(FullGrid<CombiDataType>& fg, RankType r, CommunicatorType lcomm, int n = 0) {
    dfg_->gatherFullGrid(fg, r);
  }

  DistributedFullGrid<CombiDataType>& getDistributedFullGrid(int n = 0) { return *dfg_; }

  void setZero() {}

  ~TaskConst() {
    if (dfg_ != NULL) delete dfg_;
  }

 protected:
  TaskConst() : dfg_(NULL) {}

 private:
  friend class boost::serialization::access;

  DistributedFullGrid<CombiDataType>* dfg_;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& boost::serialization::base_object<Task>(*this);
    // ar& nprocs_;
  }
};

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
    std::vector<size_t> taskIDs;
    for (size_t i = 0; i < levels.size(); i++) {
      Task* t = new TaskConst(levels[i], boundary, coeffs[i], loadmodel.get());
      tasks.push_back(t);
      taskIDs.push_back(t->getID());
    }

    IndexVector parallelization = {static_cast<IndexType>(nprocs), 1};
    // create combiparameters
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs, ncombi);
    params.setParallelization(parallelization); //TODO why??

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
      // manager.combine();

      // std::cout << "run next " << std::endl;
      // manager.runnext();
      std::cout << "combine " << std::endl;
      manager.combine();
    }

    // evaluate solution
    FullGrid<CombiDataType> fg_eval(dim, leval, boundary);
    manager.gridEval(fg_eval);
    // std::cout << "fg_eval "<< fg_eval << std::endl;

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
    while (signal != EXIT) {
      signal = pgroup.wait();
      BOOST_TEST_CHECKPOINT("Last Successful Worker Signal " + std::to_string(signal));
    }
  }

  combigrid::Stats::finalize();
  MPI_Barrier(comm);
  TestHelper::testStrayMessages(comm);
}

BOOST_FIXTURE_TEST_SUITE(reduce, TestHelper::BarrierAtEnd, *boost::unit_test::timeout(60))

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

BOOST_AUTO_TEST_SUITE_END()
