#pragma once

#define BOOST_TEST_DYN_LINK

#include <numeric>

#include <boost/serialization/export.hpp>
#include <boost/math/constants/constants.hpp>

#include "sgpp/distributedcombigrid/task/Task.hpp"

using namespace combigrid;

static constexpr double pi = boost::math::constants::pi<double>();

/**
 * functor for a sine-wave that travels through the domain
 * with periodic boundary conditions
 */
template <typename FG_ELEMENT>
class WaveFunction {
 public:
  // Notation like in Atanasov / Schnetter equation (47)
  FG_ELEMENT operator()(const std::vector<double>& coords, const std::vector<short>& wavenumber,
                        double phase = 0., double amplitude = 1.) {
    BOOST_CHECK(coords.size() > 0);
    BOOST_CHECK_EQUAL(coords.size(), wavenumber.size());
    FG_ELEMENT result = amplitude * std::cos(2. * pi *
                                                 std::inner_product(coords.begin(), coords.end(),
                                                                    wavenumber.begin(), 0.) +
                                             phase);
    BOOST_TEST_CHECKPOINT("calculated wave value");
    return result;
  }
};

/* simple task class to interpolate values of wave function
 */
class TaskInterpolateWave : public combigrid::Task {
 public:
  TaskInterpolateWave(DimType dim, LevelVector& l, std::vector<bool>& boundary, real coeff,
                      LoadModel* loadModel, double timeStep)
      : Task(dim, l, boundary, coeff, loadModel), dfg_(nullptr), timeStep_(timeStep) {
    // initialize wavenumber with increasing numbers according to dimensionality of problem
    waveNumber_ = std::vector<short>(dim);
    std::iota(waveNumber_.begin(), waveNumber_.end(), 1);
    BOOST_TEST_CHECKPOINT("TaskInterpolateWave constructor");
  }

  void init(CommunicatorType lcomm, std::vector<IndexVector> decomposition) {
    BOOST_TEST_CHECKPOINT("TaskInterpolateWave init");
    // set the cartesian decomposition of the grid that this task represents
    IndexVector p;
    if (decomposition.size() == 0) {
      auto nprocs = getCommSize(lcomm);
      p = {nprocs, 1};
    } else {
      for (const auto& d : decomposition) {
        p.push_back(d.size());
      }
    }
    BOOST_TEST_CHECKPOINT("Creating grid");
    dfg_ = new DistributedFullGrid<CombiDataType>(getDim(), getLevelVector(), lcomm, getBoundary(),
                                                  p, true, decomposition);
    BOOST_CHECK(dfg_);
    BOOST_TEST_CHECKPOINT("Created grid");

    // initialize element vector
    WaveFunction<CombiDataType> f;
    BOOST_TEST_CHECKPOINT("Created functor");
    BOOST_CHECK(&waveNumber_);
    for (IndexType li = 0; li < dfg_->getNrLocalElements(); ++li) {
      std::vector<double> coords(getDim());
      BOOST_TEST_CHECKPOINT("queried dimensions");
      dfg_->getCoordsLocal(li, coords);
      BOOST_TEST_CHECKPOINT("converted coords");
      dfg_->getData()[li] = f(coords, waveNumber_);
    }
    BOOST_TEST_CHECKPOINT("TaskInterpolateWave initialized");
  }

  void run(CommunicatorType lcomm) {
    // std::cout << "run " << getCommRank(lcomm) << std::endl;

    time_ += timeStep_;

    // update each value with new time step, assume everything travels through the domain once over each time step
    WaveFunction<CombiDataType> f;
    for (IndexType li = 0; li < dfg_->getNrLocalElements(); ++li) {
      std::vector<double> coords(getDim());
      dfg_->getCoordsLocal(li, coords);
      dfg_->getData()[li] += f(coords, waveNumber_, time_*2*pi);
    }
    BOOST_CHECK(dfg_);

    setFinished(true);

    // MPI_Barrier(lcomm);
    BOOST_TEST_CHECKPOINT("TaskInterpolateWave run");
  }

  void getFullGrid(FullGrid<CombiDataType>& fg, RankType r, CommunicatorType lcomm, int n = 0) {
    BOOST_TEST_CHECKPOINT("TaskInterpolateWave getFullGrid");
    dfg_->gatherFullGrid(fg, r);
  }

  DistributedFullGrid<CombiDataType>& getDistributedFullGrid(int n = 0) { return *dfg_; }

  void setZero() { BOOST_TEST_CHECKPOINT("TaskInterpolateWave setZero"); }

  ~TaskInterpolateWave() {
    BOOST_TEST_CHECKPOINT("TaskInterpolateWave destructor begin");
    if (dfg_ != NULL) delete dfg_;
    BOOST_TEST_CHECKPOINT("TaskInterpolateWave destructor");
  }

  real getCurrentTime() const override { return time_; }

 protected:
  TaskInterpolateWave() : dfg_(NULL) {}

 private:
  friend class boost::serialization::access;

  DistributedFullGrid<CombiDataType>* dfg_;

  combigrid::real time_ = 0.;
  combigrid::real timeStep_;
  std::vector<short> waveNumber_;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& boost::serialization::base_object<Task>(*this);
    ar& waveNumber_;
    ar& time_;
    ar& timeStep_;
  }
};
