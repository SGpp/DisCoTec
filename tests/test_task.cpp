#define BOOST_TEST_DYN_LINK
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <complex>
#include <cstdarg>
#include <iostream>
#include <vector>

#include <boost/serialization/export.hpp>
#include "../include/discotec/loadmodel/LinearLoadModel.hpp"
#include "../include/discotec/Task.hpp"
#include "../include/discotec/utils/Config.hpp"

#include "test_helper.hpp"

using namespace combigrid;

class TaskTest : public combigrid::Task {
 public:
  int test;

  TaskTest(DimType dim, const LevelVector& l, const std::vector<BoundaryType>& boundary, real coeff,
           LoadModel* loadModel, int t)
      : Task(l, boundary, coeff, loadModel), test(t) {}

  void init(CommunicatorType lcomm,
            const std::vector<IndexVector>& decomposition = std::vector<IndexVector>()) override {
    // create dummy dfg
    std::vector<int> p(getDim(), 1);
    dfg_ = new OwningDistributedFullGrid<CombiDataType>(getDim(), getLevelVector(), lcomm,
                                                        getBoundary(), p);
  }

  void run(CommunicatorType lcomm) override {}

  void getFullGrid(FullGrid<CombiDataType>& fg, RankType r, CommunicatorType lcomm, int n = 0) override {
    dfg_->gatherFullGrid(fg, r);
  }

  DistributedFullGrid<CombiDataType>& getDistributedFullGrid(size_t n = 0) override { return *dfg_; }

  void setZero() override {}

  ~TaskTest() {
    if (dfg_ != NULL) delete dfg_;
  }

 protected:
  TaskTest() : dfg_(NULL) {}

 private:
  friend class boost::serialization::access;

  OwningDistributedFullGrid<CombiDataType>* dfg_;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& boost::serialization::base_object<Task>(*this);
    ar& test;
  }
};

BOOST_CLASS_EXPORT(TaskTest)
BOOST_FIXTURE_TEST_SUITE(task, TestHelper::BarrierAtEnd, *boost::unit_test::timeout(60))

BOOST_AUTO_TEST_CASE(test) {
  int size = 8;

  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(size));
  CommunicatorType comm = TestHelper::getComm(size);
  if (comm == MPI_COMM_NULL) return;

  DimType dim = 3;
  LevelVector l(dim, 2);
  std::vector<BoundaryType> boundary(dim, 2);

  std::unique_ptr<LoadModel> loadmodel = std::unique_ptr<LinearLoadModel>(new LinearLoadModel());

  Task* t;
  if (TestHelper::getRank(comm) == 0) {
    t = new TaskTest(dim, l, boundary, 1, loadmodel.get(), 42);
  }

  // test broadcast
  Task::broadcast(&t, 0, comm);
  BOOST_CHECK(static_cast<TaskTest*>(t)->test == 42);

  // test send / receive
  static_cast<TaskTest*>(t)->test += TestHelper::getRank(comm);

  if (TestHelper::getRank(comm) != 0) {
    Task::send(t, 0, comm);
  } else {
    for (int i = 1; i < size; ++i) {
      Task::receive(&t, i, comm);
      BOOST_CHECK(static_cast<TaskTest*>(t)->test == 42 + i);
    }
  }

  assert(loadmodel);
  LevelVector test_l = t->getLevelVector();
  loadmodel->eval(test_l);
}

BOOST_AUTO_TEST_SUITE_END()
