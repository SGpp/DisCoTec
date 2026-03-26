#ifndef TASKCONSTPARABOLOID_HPP_
#define TASKCONSTPARABOLOID_HPP_

#define BOOST_TEST_DYN_LINK

#include <boost/serialization/export.hpp>
#include <optional>

#include "task/Task.hpp"
#include "utils/PowerOfTwo.hpp"

using namespace combigrid;

template <typename FG_ELEMENT>
class ParaboloidFn {
 public:
  ParaboloidFn() = default;

  FG_ELEMENT operator()(std::vector<double>& coords) {
    auto dim = static_cast<DimType>(coords.size());
    FG_ELEMENT sign;
    (dim % 2) ? sign = 1. : sign = -1.;
    FG_ELEMENT result(sign);
    for (DimType d = 0; d < dim; ++d) {
      result *= coords[d] * (coords[d] - 1.);
    }
    return result;
  }
};

/* simple task class to set all values on the grid to $levelVector_1 / levelVector_2$
 */
class TaskConstParaboloid : public combigrid::Task<> {
 public:
  TaskConstParaboloid(const LevelVector& l, const std::vector<BoundaryType>& boundary, real coeff,
                      LoadModel* loadModel)
      : Task<>(l, boundary, coeff, loadModel) {
    BOOST_TEST_CHECKPOINT("TaskConstParaboloid constructor");
  }

  void init(CommunicatorType lcomm, const std::vector<IndexVector>& decomposition) override {
    // parallelization
    auto nprocs = getCommSize(lcomm);
    std::vector<int> p(getDim(), 1);
    p[1] = nprocs;

    dfg_.emplace(makeOwningDistributedFullGrid<CombiDataType>(
        getDim(), getLevelVector(), lcomm, getBoundary(), p, false, decomposition));

    // set paraboloid function values
    ParaboloidFn<CombiDataType> f;
    auto ref = toRef(*dfg_);
    visitDFG(
        [&](auto& dfg) {
          for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
            std::vector<double> coords(getDim());
            dfg.getCoordsLocal(li, coords);
            dfg.getData()[li] = f(coords);
          }
        },
        ref);
    BOOST_TEST_CHECKPOINT("TaskConstParaboloid init");
  }

  void run(CommunicatorType lcomm) override {
    // constant run method
    ++nsteps_;
    setFinished(true);
    MPI_Barrier(lcomm);
    BOOST_TEST_CHECKPOINT("TaskConstParaboloid run");
  }

  void getFullGrid(FullGrid<CombiDataType>& fg, RankType r, CommunicatorType lcomm,
                   int n = 0) override {
    BOOST_TEST_CHECKPOINT("TaskConstParaboloid getFullGrid");
    auto ref = toRef(*dfg_);
    visitDFG([&](auto& dfg) { dfg.gatherFullGrid(fg, r); }, ref);
  }

  DistributedFullGridRef<CombiDataType> getDistributedFullGrid(size_t n = 0) override {
    BOOST_TEST_CHECKPOINT("TaskConstParaboloid getDFG");
    return toRef(*dfg_);
  }

  ConstDistributedFullGridRef<CombiDataType> getDistributedFullGrid(size_t n = 0) const override {
    BOOST_TEST_CHECKPOINT("TaskConstParaboloid getDFG const");
    return toConstRef(*dfg_);
  }

  real getCurrentTime() const override { return static_cast<real>(nsteps_); }

  void setZero() override { BOOST_CHECK(true); }

  ~TaskConstParaboloid() { BOOST_TEST_CHECKPOINT("TaskConstParaboloid destructor"); }

 protected:
  TaskConstParaboloid() {}

 private:
  friend class boost::serialization::access;

  std::optional<OwningDistributedFullGridVariant<CombiDataType>> dfg_;
  size_t nsteps_ = 0;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& boost::serialization::base_object<Task<>>(*this);
    ar & nsteps_;
  }
};

#endif  // def TASKCONSTPARABOLOID_HPP
