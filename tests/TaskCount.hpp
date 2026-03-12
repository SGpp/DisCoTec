#pragma once

#define BOOST_TEST_DYN_LINK

#include <boost/serialization/export.hpp>
#include <optional>

#include "task/Task.hpp"

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
      result += coords[d] * std::pow(10, d);
    }
    result *= static_cast<double>(nrun);
    return result;
  }
};

/* simple task class to set all values on the grid to n_timestep * sum( (d_i + 1) * x_i)
 */
class TaskCount : public combigrid::Task<> {
 public:
  TaskCount(const LevelVector& l, const std::vector<BoundaryType>& boundary, real coeff,
            LoadModel* loadModel)
      : Task<>(l, boundary, coeff, loadModel), nrun_(0) {
    BOOST_TEST_CHECKPOINT("TaskCount constructor");
  }

  void init(CommunicatorType lcomm, const std::vector<IndexVector>& decomposition) override {
    auto nprocs = getCommSize(lcomm);
    std::vector<int> p;
    if (decomposition.size() == 0) {
      p = {nprocs, 1};
    } else {
      for (const auto& d : decomposition) {
        p.push_back(static_cast<int>(d.size()));
      }
    }
    dfg_.emplace(makeOwningDistributedFullGrid<CombiDataType>(
        getDim(), getLevelVector(), lcomm, getBoundary(), p, false, decomposition));

    std::visit(
        [](auto& dfg) {
          auto element = dfg.getData();
          for (size_t i = 0; i < dfg.getNrLocalElements(); ++i) {
            element[i] = -0.;
          }
        },
        *dfg_);
    nrun_ = 0;
    BOOST_TEST_CHECKPOINT("TaskCount init");
  }

  void run(CommunicatorType lcomm) override {
    TestFnCount<CombiDataType> f;
    std::visit(
        [&](auto& dfg) {
          for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
            std::vector<double> coords(getDim());
            dfg.getCoordsLocal(li, coords);
            dfg.getData()[li] += f(coords);
          }
        },
        *dfg_);

    ++nrun_;
    BOOST_CHECK(dfg_.has_value());

    time_ += 1.;

    setFinished(true);

    MPI_Barrier(lcomm);
    BOOST_TEST_CHECKPOINT("TaskCount run");
  }

  void getFullGrid(FullGrid<CombiDataType>& fg, RankType r, CommunicatorType lcomm,
                   int n = 0) override {
    BOOST_TEST_CHECKPOINT("TaskCount getFullGrid");
    std::visit([&](auto& dfg) { dfg.gatherFullGrid(fg, r); }, *dfg_);
  }

  DistributedFullGridRef<CombiDataType> getDistributedFullGrid(size_t n = 0) override {
    return toRef<CombiDataType>(*dfg_);
  }

  void setZero() override { BOOST_TEST_CHECKPOINT("TaskCount setZero"); }

  ~TaskCount() {
    BOOST_TEST_CHECKPOINT("TaskCount destructor begin");
    BOOST_TEST_CHECKPOINT("TaskCount destructor");
  }

  CombiDataType analyticalSolution(const std::vector<real>& coords, int n = 0) const override {
    TestFnCount<CombiDataType> f;
    return f(coords, nrun_);
  }

  real getCurrentTime() const override { return time_; }

 protected:
  TaskCount() = default;

 private:
  friend class boost::serialization::access;

  std::optional<OwningDistributedFullGridVariant<CombiDataType>> dfg_;

  size_t nrun_;
  combigrid::real time_ = 0.;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& boost::serialization::base_object<Task<>>(*this);
  }
};
