/*
 * TaskExample.hpp
 *
 *  Created on: Sep 25, 2015
 *      Author: heenemo
 */

#ifndef TASKEXAMPLE_HPP_
#define TASKEXAMPLE_HPP_

#include <optional>

#include "fault_tolerance/FTUtils.hpp"
#include "fullgrid/DistributedFullGrid.hpp"
#include "task/Task.hpp"

namespace combigrid {

class TaskExample : public Task<> {
 public:
  TaskExample(DimType dim, const LevelVector& l, const std::vector<BoundaryType>& boundary,
              real coeff, LoadModel* loadModel, real dt, size_t nsteps,
              const std::vector<int>& p = std::vector<int>(0),
              FaultCriterion* faultCrit = (new StaticFaults({0, IndexVector(0), IndexVector(0)})))
      : Task(l, boundary, coeff, loadModel, faultCrit),
        dt_(dt),
        nsteps_(nsteps),
        stepsTotal_(0),
        p_(p),
        initialized_(false),
        combiStep_(0) {}

  void init(CommunicatorType lcomm,
            const std::vector<IndexVector>& decomposition = std::vector<IndexVector>()) {
    assert(!initialized_);

    int lrank = theMPISystem()->getLocalRank();

    int np;
    MPI_Comm_size(lcomm, &np);

    if (!((np > 0) && ((np & (~np + 1)) == np)))
      assert(false && "number of processes not power of two");

    DimType dim = this->getDim();
    std::vector<int> p(dim, 1);
    const LevelVector& l = this->getLevelVector();

    if (p_.size() == 0) {
      IndexType prod_p(1);

      while (prod_p != static_cast<IndexType>(np)) {
        DimType dimMaxRatio = 0;
        real maxRatio = 0.0;

        for (DimType k = 0; k < dim; ++k) {
          real ratio = std::pow(2.0, l[k]) / p[k];

          if (ratio > maxRatio) {
            maxRatio = ratio;
            dimMaxRatio = k;
          }
        }

        p[dimMaxRatio] *= 2;
        prod_p = 1;

        for (DimType k = 0; k < dim; ++k) prod_p *= p[k];
      }
    } else {
      p = p_;
    }

    if (lrank == 0) {
      std::cout << "group " << theMPISystem()->getGlobalRank() << " "
                << "computing task " << this->getID() << " with l = " << this->getLevelVector()
                << " and p = " << p << std::endl;
    }

    dfg_.emplace(
        makeOwningDistributedFullGrid<CombiDataType>(dim, l, lcomm, this->getBoundary(), p));

    std::visit(
        [&](auto& dfg) {
          auto elements = dfg.getData();
          for (size_t i = 0; i < dfg.getNrLocalElements(); ++i) {
            IndexType globalLinearIndex = dfg.getGlobalLinearIndex(i);
            std::vector<real> globalCoords(dim);
            dfg.getCoordsGlobal(globalLinearIndex, globalCoords);
            elements[i] = TaskExample::myfunction(globalCoords, 0.0);
          }
        },
        *dfg_);

    initialized_ = true;
  }

  void run(CommunicatorType lcomm) {
    assert(initialized_);

    std::visit(
        [&](auto& dfg) {
          auto elements = dfg.getData();

          for (size_t step = stepsTotal_; step < stepsTotal_ + nsteps_; ++step) {
            real time = (step + 1) * dt_;

            for (size_t i = 0; i < dfg.getNrLocalElements(); ++i) {
              IndexType globalLinearIndex = dfg.getGlobalLinearIndex(i);
              std::vector<real> globalCoords(this->getDim());
              dfg.getCoordsGlobal(globalLinearIndex, globalCoords);
              elements[i] = TaskExample::myfunction(globalCoords, time);
            }
          }
        },
        *dfg_);

    stepsTotal_ += nsteps_;
    this->setFinished(true);
    decideToKill();
    ++combiStep_;
  }

  void getFullGrid(FullGrid<CombiDataType>& fg, RankType r, CommunicatorType lcomm, int n = 0) {
    std::visit([&](auto& dfg) { dfg.gatherFullGrid(fg, r); }, *dfg_);
  }

  DistributedFullGridRef<CombiDataType> getDistributedFullGrid(size_t n) override {
    return toRef<CombiDataType>(*dfg_);
  }

  ConstDistributedFullGridRef<CombiDataType> getDistributedFullGrid(size_t n) const override {
    return toConstRef<CombiDataType>(*dfg_);
  }

  static real myfunction(std::vector<real>& coords, real t) {
    real u = std::cos(M_PI * t);

    for (size_t d = 0; d < coords.size(); ++d) u *= std::cos(2.0 * M_PI * coords[d]);

    return u;
  }

  inline void setStepsTotal(size_t stepsTotal);

  inline void setZero() {}

  ~TaskExample() = default;

 protected:
  TaskExample() : initialized_(false), combiStep_(0) {}

 private:
  friend class boost::serialization::access;

  real dt_;
  size_t nsteps_ = 0;
  size_t stepsTotal_;
  std::vector<int> p_;

  bool initialized_;
  size_t combiStep_;

  std::optional<OwningDistributedFullGridVariant<CombiDataType>> dfg_;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& boost::serialization::base_object<Task<>>(*this);

    ar & dt_;
    ar & nsteps_;
    ar & stepsTotal_;
    ar & p_;
  }

  void decideToKill() {
    using namespace std::chrono;

    int globalRank = theMPISystem()->getGlobalRank();

    if (combiStep_ != 0 && faultCriterion_->failNow(combiStep_, -1.0, globalRank)) {
      std::cout << "Rank " << globalRank << " failed at iteration " << combiStep_ << std::endl;
      StatusType status = PROCESS_GROUP_FAIL;
      MASTER_EXCLUSIVE_SECTION {
        simft::Sim_FT_MPI_Send(&status, 1, MPI_INT, theMPISystem()->getManagerRank(),
                               TRANSFER_STATUS_TAG, theMPISystem()->getGlobalCommFT());
      }
      theMPISystem()->sendFailedSignal();
      simft::Sim_FT_kill_me();
    }
  }
};

inline void TaskExample::setStepsTotal(size_t stepsTotal) { stepsTotal_ = stepsTotal; }

}  // namespace combigrid

#endif /* TASKEXAMPLE_HPP_ */
