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

using namespace combigrid;

class TestingTaskRescheduler : public TaskRescheduler {
 public:
  std::vector<std::pair<LevelVector, int>> eval(
      const std::map<LevelVector, int>& levelVectorToProcessGroupIndex,
      const std::map<LevelVector, unsigned long>& levelVectorToTaskDuration,
      LoadModel *loadModel) override {
    // Find arbitrary tasks to reschedule! (but at least 1 must be left per process group)

    // Find groups with more than 1 task
    // * Get Task ids
    // Move them to the next group

    std::vector<std::pair<LevelVector, int>> result{};
    std::map<int, bool> processGroupIndexToHasAtLeastOneTask{};
    std::set<int> activeProcessGroupIndices;

    // initialize for every process group index the value to false
    for (const auto& p : levelVectorToProcessGroupIndex) {
      processGroupIndexToHasAtLeastOneTask.insert({p.second, false});
      activeProcessGroupIndices.insert(p.second);
    }

    for (const auto& p : levelVectorToProcessGroupIndex) {
      if (processGroupIndexToHasAtLeastOneTask[p.second]) {
        // reschedule task into random different process group:
        auto r = rand() % activeProcessGroupIndices.size();
        auto iterator = activeProcessGroupIndices.begin();
        std::advance(iterator, r);
        auto randomProcessGroup = *iterator;
        result.push_back({p.first, randomProcessGroup});
      } else {
        processGroupIndexToHasAtLeastOneTask[p.second] = true;
      }
    }

    
    return result;
  }
};

/* simple task class to set all values on the grid to $levelVector_1 / levelVector_2$
 */
class TestingTask : public combigrid::Task {
 public:
  TestingTask(LevelVector& l, std::vector<bool>& boundary, real coeff, LoadModel* loadModel)
      : Task(2, l, boundary, coeff, loadModel){}

  void init(CommunicatorType lcomm, std::vector<IndexVector> decomposition) override {
    // parallelization
    // assert(dfg_ == nullptr);
    long nprocs = getCommSize(lcomm);
    IndexVector p = {nprocs,1};

    dfg_ = new DistributedFullGrid<CombiDataType>(getDim(), getLevelVector(), lcomm, getBoundary(),
                                                  p, false, decomposition);

    std::vector<CombiDataType>& elements = dfg_->getElementVector();
    for (auto& element : elements) {
      element = 0; // default state is 0
    }
  }

  void run(CommunicatorType lcomm) override {
    
    ++valueToPersist_;

    std::cout << "run " << getCommRank(lcomm) << std::endl;    
    
    std::vector<CombiDataType>& elements = dfg_->getElementVector();
    for (auto& element : elements) {
      element = 10; // after run was executed the state is 10
    }
    BOOST_CHECK(dfg_);

    setFinished(true);
    
    MPI_Barrier(lcomm);
  }

  void getFullGrid(FullGrid<CombiDataType>& fg, RankType r, CommunicatorType lcomm, int n = 0) override {
    dfg_->gatherFullGrid(fg, r);
  }

  DistributedFullGrid<CombiDataType>& getDistributedFullGrid(int n = 0) override { return *dfg_; }

  void setZero() override {}

  ~TestingTask() override {
    if (dfg_ != nullptr) delete dfg_;
  }

  int valueToPersist_{0};

 protected:
  TestingTask() {}

 private:
  friend class boost::serialization::access;

  DistributedFullGrid<CombiDataType>* dfg_{nullptr};


  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar & boost::serialization::base_object<Task>(*this);
    ar & valueToPersist_;
  }
};

BOOST_CLASS_EXPORT(TestingTask)

bool tasksContainSameValue(const TaskContainer& tasks) {
  std::cout << "task size" << tasks.size() << "\n";
  if (tasks.size() <= 1) {
    return true;
  }

  auto firstValue = dynamic_cast<TestingTask *>(tasks[0])->valueToPersist_;

  for (auto i = tasks.cbegin() + 1; i < tasks.cend(); ++i) {
    if (firstValue != dynamic_cast<TestingTask *>(*i)->valueToPersist_) {
      return false;
    }
  }
  return true;
}

void checkRescheduling(size_t ngroup = 1, size_t nprocs = 1) {
  size_t size = ngroup * nprocs + 1;
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(size));

  CommunicatorType comm = TestHelper::getComm(size);
  if (comm == MPI_COMM_NULL) {
    return;
  }

  combigrid::Stats::initialize();

  theMPISystem()->initWorldReusable(comm, ngroup, nprocs);
  // theMPISystem()->init(ngroup, nprocs);

  WORLD_MANAGER_EXCLUSIVE_SECTION {
    ProcessGroupManagerContainer pgroups;
    for (int i = 0; i < ngroup; ++i) {
      int pgroupRootID(i);
      pgroups.emplace_back(std::make_shared<ProcessGroupManager>(pgroupRootID));
    }

    auto loadmodel = std::unique_ptr<LoadModel>(new LinearLoadModel());
    auto rescheduler = std::unique_ptr<TaskRescheduler>(new TestingTaskRescheduler());

    DimType dim = 2;
    LevelVector lmin(dim, 2);
    LevelVector lmax(dim, 4);

    size_t ncombi = 2;
    std::vector<bool> boundary(dim, true);

    CombiMinMaxScheme combischeme(dim, lmin, lmax);
    combischeme.createAdaptiveCombischeme();

    std::vector<LevelVector> levels = combischeme.getCombiSpaces();
    std::vector<combigrid::real> coeffs = combischeme.getCoeffs();

    // create Tasks
    TaskContainer tasks;
    std::vector<int> taskIDs;
    for (size_t i = 0; i < levels.size(); i++) {
      Task* t = new TestingTask(levels[i], boundary, coeffs[i], loadmodel.get());
      tasks.push_back(t);
      taskIDs.push_back(t->getID());
    }

    // create combiparameters
    //
    // Reduce combination dims lmin and lmax are 0!!
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs, ncombi);
    params.setParallelization({static_cast<IndexType>(nprocs), 1});


    // create abstraction for Manager
    ProcessManager manager{pgroups, tasks, params, std::move(loadmodel), std::move(rescheduler)};

    // the combiparameters are sent to all process groups before the
    // computations start
    manager.updateCombiParameters(); //TODO move to manager constructor or runfirst?

    /* distribute task according to load model and start computation for
     * the first time */
    BOOST_TEST_CHECKPOINT("run first");
    manager.runfirst();


    for (size_t it = 0; it < ncombi - 1; ++it) {
      BOOST_TEST_CHECKPOINT("combine");
      manager.combine();

      BOOST_TEST_CHECKPOINT("reschedule");
      manager.reschedule();

      BOOST_TEST_CHECKPOINT("run next");
      manager.runnext();
    }
    manager.combine();

    manager.exit();
  }
  else {
    BOOST_TEST_CHECKPOINT("Worker startet");
    ProcessGroupWorker pgroup;
    SignalType signal = -1;
    while (signal != EXIT){ 
      signal = pgroup.wait();

      // test conditions:
      // ================

      BOOST_REQUIRE(tasksContainSameValue(pgroup.getTasks()));

      for (auto& t : pgroup.getTasks()) {
        std::vector<CombiDataType>& elements = t->getDistributedFullGrid().getElementVector();
        for (auto e : elements) {
          // Elements need to be always equal to 10 because
          // * first run: run is executed
          // * next run: run is executed
          // * combination: element values do not change (every element is 
          //                equal to 10)
          // * Task was added: element values are restored
          // * Task was removed: all remaining elements were already checked in 
          //                     a previous iteration
          BOOST_REQUIRE(e == 10);
        }
      }
    }
  }

  combigrid::Stats::finalize();
  MPI_Barrier(comm);
  TestHelper::testStrayMessages(comm);
}

BOOST_AUTO_TEST_SUITE(rescheduling)

BOOST_AUTO_TEST_CASE(test_1, *boost::unit_test::tolerance(TestHelper::higherTolerance) *
                                 boost::unit_test::timeout(60)) {
  std::cout << "rescheduling/test_1"<< std::endl;
  checkRescheduling(3,1);
}

BOOST_AUTO_TEST_CASE(test_2, *boost::unit_test::tolerance(TestHelper::higherTolerance) *
                                 boost::unit_test::timeout(60)) {
  std::cout << "rescheduling/test_2"<< std::endl;
  checkRescheduling(3,2);
}


BOOST_AUTO_TEST_SUITE_END()
