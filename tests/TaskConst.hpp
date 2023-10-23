#ifndef TASKCONST_HPP_
#define TASKCONST_HPP_

#define BOOST_TEST_DYN_LINK

#include <boost/serialization/export.hpp>
#include "task/Task.hpp"

using namespace combigrid;

/* simple task class to set all values on the grid to $levelVector_1 / levelVector_2$
 */
class TaskConst : public combigrid::Task {
 public:
  TaskConst(LevelVector& l, std::vector<bool>& boundary, real coeff, LoadModel* loadModel)
      : Task(l, boundary, coeff, loadModel) {
    assert(l.size() == 2);
  }

  void init(CommunicatorType lcomm, std::vector<IndexVector> decomposition) {
    // parallelization
    // assert(dfg_ == nullptr);
    auto nprocs = getCommSize(lcomm);
    std::vector<int> p = {nprocs,1};

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

    dfg_ = new OwningDistributedFullGrid<CombiDataType>(getDim(), getLevelVector(), lcomm,
                                                        getBoundary(), p, false, decomposition);
    auto elements = dfg_->getData();
    for (size_t i = 0; i < dfg_->getNrLocalElements(); ++i) {
      elements[i] = 10;
    }
    BOOST_CHECK(true);
  }

  void run(CommunicatorType lcomm) {
    auto elements = dfg_->getData();
    for (size_t i = 0; i < dfg_->getNrLocalElements(); ++i) {
      elements[i] = getLevelVector()[0] / (double)getLevelVector()[1];
    }

    BOOST_CHECK(dfg_);

    setFinished(true);
    
    MPI_Barrier(lcomm);
    BOOST_CHECK(true);
  }

  void getFullGrid(FullGrid<CombiDataType>& fg, RankType r, CommunicatorType lcomm, int n = 0) {
    BOOST_CHECK(true);
    dfg_->gatherFullGrid(fg, r);
  }

  DistributedFullGrid<CombiDataType>& getDistributedFullGrid(int n = 0) { return *dfg_; }

  void setZero() { BOOST_CHECK(true); }

  ~TaskConst() {
    if (dfg_ != NULL) delete dfg_;
  }

 protected:
  TaskConst() : dfg_(NULL) {}

 private:
  friend class boost::serialization::access;

  OwningDistributedFullGrid<CombiDataType>* dfg_;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& boost::serialization::base_object<Task>(*this);
    // ar& nprocs_;
  }
};

#endif // def TASKCONST_HPP_    BOOST_CHECK(true);
