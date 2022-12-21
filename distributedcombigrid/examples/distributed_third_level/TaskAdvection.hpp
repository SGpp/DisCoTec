#pragma once

#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGrid.hpp"
#include "sgpp/distributedcombigrid/task/Task.hpp"
#include "sgpp/distributedcombigrid/mpi/MPIMemory.hpp"

namespace combigrid {

// exact solution
class TestFn {
 public:
  // function value
  double operator()(std::vector<double>& coords, double t) const {
    double exponent = 0;
    const double sigma = 1./3.;
    const double sigmaSquared = sigma * sigma;
    const double sigmaSquaredInv = 1. / sigmaSquared;
    for (DimType d = 0; d < coords.size(); ++d) {
      coords[d] = std::fmod(1.0 + std::fmod(coords[d] - t, 1.0), 1.0);
      assert(coords[d] >= 0.);
      assert(coords[d] <= 1.);
      exponent -= std::pow(coords[d] - 0.5, 2);
    }
    exponent *= sigmaSquaredInv;
    // leave out normalization, such that maximum is always 1
    return std::exp(exponent);  // / std::sqrt( std::pow(2*pi*sigmaSquared, coords.size()));
  }
};

class TaskAdvection : public Task {
 public:
  /* if the constructor of the base task class is not sufficient we can provide an
   * own implementation. here, we add dt, nsteps, and p as a new parameters.
   */
  TaskAdvection(DimType dim, LevelVector& l, std::vector<BoundaryType>& boundary, real coeff,
                LoadModel* loadModel, real dt, size_t nsteps,
                std::vector<int> p = std::vector<int>(0),
                FaultCriterion* faultCrit = (new StaticFaults({0, IndexVector(0), IndexVector(0)})))
      : Task(dim, l, boundary, coeff, loadModel, faultCrit),
        dt_(dt),
        nsteps_(nsteps),
        p_(p),
        initialized_(false),
        stepsTotal_(0),
        dfg_(nullptr) {
          for (const auto& b : boundary) {
            assert(b == 1);
          }
        }

  void init(CommunicatorType lcomm,
            std::vector<IndexVector> decomposition = std::vector<IndexVector>()) {
    assert(!initialized_);
    assert(dfg_ == NULL);

    auto lrank = theMPISystem()->getLocalRank();
    DimType dim = this->getDim();
    const LevelVector& l = this->getLevelVector();

    // create local subgrid on each process
    dfg_ = new DistributedFullGrid<CombiDataType>(dim, l, lcomm, this->getBoundary(), p_, false,
                                                  decomposition);
    if (phi_ == nullptr) {
      phi_ = new std::vector<CombiDataType>(dfg_->getNrLocalElements());
    }

    std::vector<double> h = dfg_->getGridSpacing();
    auto sumOneOverH = 0.;
    for (const auto& h_x : h) {
      sumOneOverH += 1. / h_x;
    }
    if (!(dt_ * sumOneOverH < 1.)) {
      throw std::runtime_error("CFL condition not satisfied!");
    }

    TestFn f;
    for (IndexType li = 0; li < dfg_->getNrLocalElements(); ++li) {
      static std::vector<double> coords(this->getDim());
      dfg_->getCoordsLocal(li, coords);
      dfg_->getData()[li] = f(coords, 0.);
    }

    initialized_ = true;
  }

  /* this is were the application code kicks in and all the magic happens.
   * do whatever you have to do, but make sure that your application uses
   * only lcomm or a subset of it as communicator.
   * important: don't forget to set the isFinished flag at the end of the computation.
   */
  void run(CommunicatorType lcomm) {
    assert(initialized_);
    // dfg_->print(std::cout);

    std::vector<CombiDataType> velocity(this->getDim(), 1);

    std::vector<double> h = dfg_->getGridSpacing();
    const auto& fullOffsets = dfg_->getLocalOffsets();

    phi_->resize(dfg_->getNrLocalElements());
    std::fill(phi_->begin(), phi_->end(), 0.);

    for (size_t i = 0; i < nsteps_; ++i) {
      // compute the gradient in the original dfg_, then update into phi_ and
      // swap at the end of each time step
      auto& u_dot_dphi = *phi_;

      for (unsigned int d = 0; d < this->getDim(); ++d) {
        // to update the values in the "lowest" layer, we need the ghost values from the lower
        // neighbor
        IndexVector subarrayExtents;
        std::vector<CombiDataType> phi_ghost{};
        dfg_->exchangeGhostLayerUpward(d, subarrayExtents, phi_ghost);
        int subarrayOffset = 1;
        IndexVector subarrayOffsets(this->getDim());
        for (DimType d_j = 0; d_j < dim_; ++d_j) {
          subarrayOffsets[d_j] = subarrayOffset;
          subarrayOffset *= subarrayExtents[d_j];
        }

        for (IndexType li = 0; li < dfg_->getNrLocalElements(); ++li) {
          // calculate local axis index of backward neighbor
          static IndexVector locAxisIndex(this->getDim());
          dfg_->getLocalVectorIndex(li, locAxisIndex);
          // TODO can be unrolled into ghost and other part, avoiding if-statement
          CombiDataType phi_neighbor = 0.;
          if (locAxisIndex[d] == 0){
            assert(phi_ghost.size() > 0);
            // then use values from boundary exchange
            IndexType gli = 0;
            for (DimType d_j = 0; d_j < this->getDim(); ++d_j) {
              gli = gli + subarrayOffsets[d_j] * locAxisIndex[d_j];
            }
            assert(gli > -1);
            assert(gli < phi_ghost.size());
            phi_neighbor = phi_ghost[gli];
          } else{
            // --locAxisIndex[d];
            // IndexType lni = dfg_->getLocalLinearIndex(locAxisIndex);
            // assert(lni == (li - fullOffsets[d]));
            phi_neighbor = dfg_->getElementVector()[li - fullOffsets[d]];
          }
          // calculate gradient of phi with backward differential quotient
          auto dphi = (dfg_->getElementVector()[li] - phi_neighbor) / h[d];

          u_dot_dphi[li] += velocity[d] * dphi;
        }
      }
      for (IndexType li = 0; li < dfg_->getNrLocalElements(); ++li) {
        (*phi_)[li] = dfg_->getElementVector()[li] - u_dot_dphi[li] * dt_;
      }
      phi_->swap(dfg_->getElementVector());
    }
    stepsTotal_ += nsteps_;

    this->setFinished(true);
  }

  /* this function evaluates the combination solution on a given full grid.
   * here, a full grid representation of your task's solution has to be created
   * on the process of lcomm with the rank r.
   * typically this would require gathering your (in whatever way) distributed
   * solution on one process and then converting it to the full grid representation.
   * the DistributedFullGrid class offers a convenient function to do this.
   */
  void getFullGrid(FullGrid<CombiDataType>& fg, RankType r, CommunicatorType lcomm, int n = 0) {
    assert(fg.getLevels() == dfg_->getLevels());

    dfg_->gatherFullGrid(fg, r);
  }

  DistributedFullGrid<CombiDataType>& getDistributedFullGrid(int n = 0) { return *dfg_; }

  void setZero() {
    dfg_->setZero();
    std::fill(phi_->begin(), phi_->end(), 0.);
  }

  ~TaskAdvection() {
    if (dfg_ != nullptr) delete dfg_;
    if (phi_ != nullptr) delete phi_;
  }

  CombiDataType analyticalSolution(const std::vector<real>& coords, int n = 0) const override {
    assert(n == 0);
    auto coordsCopy = coords;
    TestFn f;
    return f(coordsCopy, stepsTotal_*dt_);
  }

 protected:
  /* if there are local variables that have to be initialized at construction
   * you have to do it here. the worker processes will create the task using
   * this constructor before overwriting the variables that are set by the
   * manager. here we need to set the initialized variable to make sure it is
   * set to false. */
  TaskAdvection() : initialized_(false), stepsTotal_(0), dfg_(nullptr) {}

 private:
  friend class boost::serialization::access;

  // new variables that are set by manager. need to be added to serialize
  real dt_;
  size_t nsteps_;
  std::vector<int> p_;

  // pure local variables that exist only on the worker processes
  bool initialized_;
  size_t stepsTotal_;
  DistributedFullGrid<CombiDataType>* dfg_;
  static std::vector<CombiDataType>* phi_;

  /**
   * The serialize function has to be extended by the new member variables.
   * However this concerns only member variables that need to be exchanged
   * between manager and workers. We do not need to add "local" member variables
   * that are only needed on either manager or worker processes.
   * For serialization of the parent class members, the class must be
   * registered with the BOOST_CLASS_EXPORT macro.
   */
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    // handles serialization of base class
    ar& boost::serialization::base_object<Task>(*this);

    // add our new variables
    ar& dt_;
    ar& nsteps_;
    ar& p_;
  }
};

std::vector<CombiDataType>* TaskAdvection::phi_ = nullptr;

}  // namespace combigrid
