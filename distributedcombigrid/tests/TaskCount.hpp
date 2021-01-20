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
  FG_ELEMENT operator()(std::vector<double>& coords, size_t nrun = 1) {
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
      : Task(dim, l, boundary, coeff, loadModel){}

  void init(CommunicatorType lcomm, std::vector<IndexVector> decomposition) {

    long nprocs = getCommSize(lcomm);
    IndexVector p = {nprocs,1};
    dfg_ = new DistributedFullGrid<CombiDataType>(getDim(), getLevelVector(), lcomm, getBoundary(),
                                                  p, false, decomposition);

    std::vector<CombiDataType>& elements = dfg_->getElementVector();
    for (auto& element : elements) {
      element = -0.;
    }
    BOOST_CHECK(true);
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

  ~TaskCount() {
    if (dfg_ != NULL) delete dfg_;
  }

 protected:
  TaskCount() : dfg_(NULL) {}

 private:
  friend class boost::serialization::access;

  DistributedFullGrid<CombiDataType>* dfg_;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& boost::serialization::base_object<Task>(*this);
    // ar& nprocs_;
  }
};
