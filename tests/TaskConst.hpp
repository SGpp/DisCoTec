#ifndef TASKCONST_HPP_
#define TASKCONST_HPP_

#define BOOST_TEST_DYN_LINK

#include <boost/serialization/export.hpp>
#include <optional>

#include "discotec/task/Task.hpp"

using namespace combigrid;

/* simple task class to set all values on the grid to $levelVector_1 / levelVector_2$
 */
class TaskConst : public combigrid::Task<> {
  static constexpr DimType DIM = 2;

 public:
  TaskConst(LevelVector& l, std::vector<bool>& boundary, real coeff, LoadModel* loadModel)
      : Task<>(l, boundary, coeff, loadModel) {
    assert(l.size() == DIM);
  }

  void init(CommunicatorType lcomm, const std::vector<IndexVector>& decomposition) {
    // parallelization
    auto nprocs = getCommSize(lcomm);
    std::vector<int> p = {nprocs, 1};

    dfg_.emplace(DIM, getLevelVector(), lcomm, getBoundary(), p, false, decomposition);
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

    BOOST_CHECK(dfg_.has_value());

    setFinished(true);

    MPI_Barrier(lcomm);
    BOOST_CHECK(true);
  }

  void getFullGrid(FullGrid<CombiDataType>& fg, RankType r, CommunicatorType lcomm, int n = 0) {
    BOOST_CHECK(true);
    dfg_->gatherFullGrid(fg, r);
  }

  DistributedFullGridRef<CombiDataType> getDistributedFullGrid(size_t n = 0) override {
    return std::ref(static_cast<DistributedFullGrid<CombiDataType, DIM>&>(*dfg_));
  }

  ConstDistributedFullGridRef<CombiDataType> getDistributedFullGrid(size_t n = 0) const override {
    return std::cref(static_cast<const DistributedFullGrid<CombiDataType, DIM>&>(*dfg_));
  }

  void setZero() { BOOST_CHECK(true); }

  ~TaskConst() {}

 protected:
  TaskConst() {}

 private:
  friend class boost::serialization::access;

  std::optional<OwningDistributedFullGrid<CombiDataType, DIM>> dfg_;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& boost::serialization::base_object<Task<>>(*this);
  }
};

#endif  // def TASKCONST_HPP_
