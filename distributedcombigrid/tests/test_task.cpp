#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <mpi.h>
#include <iostream>
#include <complex>
#include <cstdarg>
#include <vector>

#include <boost/serialization/export.hpp>
#include "sgpp/distributedcombigrid/utils/Config.hpp"
#include "sgpp/distributedcombigrid/task/Task.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LinearLoadModel.hpp"

#include "test_helper.hpp"

class TaskTest : public combigrid::Task {
public:
  int test;

  TaskTest(DimType dim, LevelVector& l, std::vector<bool>& boundary,
              real coeff, LoadModel* loadModel, int t) :
    Task(dim, l, boundary, coeff, loadModel), test(t) {
  }

  void init(CommunicatorType lcomm, std::vector<IndexVector> decomposition = std::vector<IndexVector>()) {
    // create dummy dfg
    IndexVector p(getDim(), 1);
    dfg_ = new DistributedFullGrid<CombiDataType>(getDim(), getLevelVector(),
                                                  lcomm, getBoundary(), p);
  }

  void run(CommunicatorType lcomm) {
  }

  void getFullGrid(FullGrid<CombiDataType>& fg, RankType r,
                   CommunicatorType lcomm, int n = 0) {
    dfg_->gatherFullGrid(fg, r);
  }

  DistributedFullGrid<CombiDataType>& getDistributedFullGrid(int n = 0) {
    return *dfg_;
  }

  void setZero() {

  }

protected:
  TaskTest() :
    dfg_(NULL) {
  }

  ~TaskTest() {
    if (dfg_ != NULL)
      delete dfg_;
  }

private:
  friend class boost::serialization::access;

  DistributedFullGrid<CombiDataType>* dfg_;

  template<class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& boost::serialization::base_object<Task>(*this);
    ar& test;
  }
};

BOOST_CLASS_EXPORT(TaskTest)

BOOST_AUTO_TEST_SUITE(task)

BOOST_AUTO_TEST_CASE(test) {
  int size = 8;

  BOOST_REQUIRE(TestHelper::checkNumProcs(size));
  CommunicatorType comm = TestHelper::getComm(size);
  if (comm == MPI_COMM_NULL) return;

  DimType dim = 3;
  LevelVector l(dim, 2);
  std::vector<bool> boundary(dim, true);
  LoadModel* loadmodel = new combigrid::LinearLoadModel();

  Task* t;
  if (TestHelper::getRank(comm) == 0) {
    t = new TaskTest(dim, l, boundary, 1, loadmodel, 42);
  }

  // test broadcast
  Task::broadcast(&t, 0, comm);
  BOOST_CHECK(static_cast<TaskTest*>(t)->test == 42);

  // test send / receive
  static_cast<TaskTest*>(t)->test += TestHelper::getRank(comm);

  if (TestHelper::getRank(comm) != 0) {
    Task::send(&t, 0, comm);
  } else {
    for (int i = 1; i < size; ++i) {
      Task::receive(&t, i, comm);
      BOOST_CHECK(static_cast<TaskTest*>(t)->test == 42 + i);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
