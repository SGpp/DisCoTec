#pragma once

#define BOOST_TEST_DYN_LINK

#include <boost/serialization/export.hpp>
#include "sgpp/distributedcombigrid/task/Task.hpp"

using namespace combigrid;


/**
 * functor for test function $f(x) = \sum_{i=0}^d x_i * (i+1)$
 * which maps to points on a hyperplane
 */
template <typename FG_ELEMENT>
class TestFnCount {
 public:
  // function value
  FG_ELEMENT operator()(const std::vector<double>& coords, size_t nrun = 1) {
    FG_ELEMENT result = 0.;
    for (DimType d = 0; d < coords.size(); ++d) {
      result += coords[d] * std::pow(10,d);
    }
    result *= static_cast<double>(nrun);
    return result;
  }
};

/* simple task class to set all values on the grid to n_timestep * sum( (d_i + 1) * x_i)
 */
class TaskCount : public combigrid::Task {
 public:
  TaskCount(DimType dim, LevelVector& l, std::vector<bool>& boundary, real coeff, LoadModel* loadModel)
      : Task(dim, l, boundary, coeff, loadModel), dfg_(nullptr), nrun_(0) {
        BOOST_TEST_CHECKPOINT("TaskCount constructor");
      }

  void init(CommunicatorType lcomm, std::vector<IndexVector> decomposition) {

    long nprocs = getCommSize(lcomm);
    std::vector<int> p;
    if (decomposition.size() == 0) {
      p = {nprocs,1};
    } else {
      for (const auto& d : decomposition) {
        p.push_back(d.size());
      }
    }
    dfg_ = new DistributedFullGrid<CombiDataType>(getDim(), getLevelVector(), lcomm, getBoundary(), p, true, decomposition);

    std::vector<CombiDataType>& elements = dfg_->getElementVector();
    for (auto& element : elements) {
      element = -0.;
    }
    nrun_ = 0;
    BOOST_TEST_CHECKPOINT("TaskCount init");
  }

  void run(CommunicatorType lcomm) {

    // std::cout << "run " << getCommRank(lcomm) << std::endl;

    // increase each value by sum( (d_i+1) * x_i)
    TestFnCount<CombiDataType> f;
    for (IndexType li = 0; li < dfg_->getNrLocalElements(); ++li) {
      std::vector<double> coords(getDim());
      dfg_->getCoordsLocal(li, coords);
      dfg_->getData()[li] += f(coords);
    }

    ++nrun_;
    BOOST_CHECK(dfg_);

    time_ += 1.;

    setFinished(true);

    MPI_Barrier(lcomm);
    BOOST_TEST_CHECKPOINT("TaskCount run");
  }

  void getFullGrid(FullGrid<CombiDataType>& fg, RankType r, CommunicatorType lcomm, int n = 0) {
    BOOST_TEST_CHECKPOINT("TaskCount getFullGrid");
    dfg_->gatherFullGrid(fg, r);
  }

  DistributedFullGrid<CombiDataType>& getDistributedFullGrid(int n = 0) { return *dfg_; }

  void setZero() {
    BOOST_TEST_CHECKPOINT("TaskCount setZero");
  }

  ~TaskCount() {
    BOOST_TEST_CHECKPOINT("TaskCount destructor begin");
    if (dfg_ != NULL) delete dfg_;
    BOOST_TEST_CHECKPOINT("TaskCount destructor");
  }

  CombiDataType analyticalSolution(const std::vector<real>& coords, int n = 0) const override {
    TestFnCount<CombiDataType> f;
    return f(coords, nrun_);
  }

  real getCurrentTime() const override {
    return time_;
  }

 protected:
  TaskCount() : dfg_(NULL) {}

 private:
  friend class boost::serialization::access;

  DistributedFullGrid<CombiDataType>* dfg_;

  size_t nrun_;
  combigrid::real time_ = 0.;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& boost::serialization::base_object<Task>(*this);
    // ar& nprocs_;
  }
};
