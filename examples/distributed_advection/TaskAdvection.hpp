#pragma once

#include "fullgrid/DistributedFullGrid.hpp"
#include "task/Task.hpp"

namespace combigrid {

// exact solution
class TestFn {
 public:
  // function value
  double operator()(const std::vector<double>& coords, double t) const {
    double exponent = 0;
    for (DimType d = 0; d < coords.size(); ++d) {
      auto shifted_coords = std::fmod(1.0 + std::fmod(coords[d] - t, 1.0), 1.0);
      assert(shifted_coords >= 0.);
      assert(shifted_coords <= 1.);
      exponent -= std::pow(shifted_coords - 0.5, 2);
    }
    exponent *= sigmaSquaredInv_;
    // leave out normalization, such that maximum is always 1
    return std::exp(exponent);  // / std::sqrt( std::pow(2*pi*sigmaSquared, coords.size()));
  }

 private:
  static constexpr double sigmaSquaredInv_ = 1. / ((1. / 3.) * (1. / 3.));
};

class TaskAdvection : public Task {
 public:
  /* if the constructor of the base task class is not sufficient we can provide an
   * own implementation. here, we add dt, nsteps, and p as a new parameters.
   */
  TaskAdvection(const LevelVector& l, const std::vector<BoundaryType>& boundary, real coeff,
                LoadModel* loadModel, real dt, size_t nsteps,
                std::vector<int> p = std::vector<int>(0),
                FaultCriterion* faultCrit = (new StaticFaults({0, IndexVector(0), IndexVector(0)})))
      : Task(l, boundary, coeff, loadModel, faultCrit),
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
            std::vector<IndexVector> decomposition = std::vector<IndexVector>()) override {
    assert(!initialized_);
    assert(dfg_ == NULL);

    DimType dim = this->getDim();
    const LevelVector& l = this->getLevelVector();

    // create local subgrid on each process
    dfg_ = new OwningDistributedFullGrid<CombiDataType>(dim, l, lcomm, this->getBoundary(), p_,
                                                        false, decomposition);
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
  void run(CommunicatorType lcomm) override {
    assert(initialized_);
    // dfg_->print(std::cout);
    const auto numLocalElements = dfg_->getNrLocalElements();

    const std::vector<CombiDataType> velocity(this->getDim(), 1);

    const std::vector<double> oneOverH = dfg_->getInverseGridSpacing();
    const auto& fullOffsets = dfg_->getLocalOffsets();

    phi_->resize(numLocalElements);

    for (size_t i = 0; i < nsteps_; ++i) {
      std::memset(phi_->data(), 0, phi_->size() * sizeof(CombiDataType));
      // compute the gradient in the original dfg_, then update into phi_ and
      // swap at the end of each time step
      auto& u_dot_dphi = *phi_;
      auto const ElementVector = dfg_->getData();
      for (unsigned int d = 0; d < this->getDim(); ++d) {
        static std::vector<int> subarrayExtents;
        std::vector<CombiDataType> phi_ghost{};
        dfg_->exchangeGhostLayerUpward(d, subarrayExtents, phi_ghost);

        // update all values; this will also (wrongly) update the lowest layer's values
        for (IndexType li = 0; li < numLocalElements; ++li) {
#ifndef NDEBUG
          IndexVector locAxisIndex(this->getDim());
          dfg_->getLocalVectorIndex(li, locAxisIndex);
          if (locAxisIndex[d] > 0) {
            --locAxisIndex[d];
            IndexType lni = dfg_->getLocalLinearIndex(locAxisIndex);
            assert(lni == (li - fullOffsets[d]));
          }
#endif  // NDEBUG
          // ensure modulo is positive, cf. https://stackoverflow.com/a/12277233
          const auto neighborLinearIndex = (li - fullOffsets[d] + numLocalElements) % numLocalElements;
          assert(neighborLinearIndex >= 0);
          assert(neighborLinearIndex < numLocalElements);
          CombiDataType phi_neighbor = ElementVector[neighborLinearIndex];
          auto dphi = (ElementVector[li] - phi_neighbor) * oneOverH[d];
          u_dot_dphi[li] += velocity[d] * dphi;
        }
        // iterate the lowest layer and update the values, compensating for the wrong update
        // before
        assert(dfg_->getNrLocalElements() / dfg_->getLocalSizes()[d] == phi_ghost.size());
        const auto& stride = dfg_->getLocalOffsets()[d];
        const IndexType jump = stride * dfg_->getLocalSizes()[d];
        const IndexType numberOfPolesHigherDimensions = dfg_->getNrLocalElements() / jump;
        IndexType dfgLowestLayerIteratedIndex;
        IndexType ghostIndex = 0;
        for (IndexType nHigher = 0; nHigher < numberOfPolesHigherDimensions; ++nHigher) {
          dfgLowestLayerIteratedIndex = nHigher * jump;  // local linear index
          for (IndexType nLower = 0; nLower < dfg_->getLocalOffsets()[d];
               ++nLower && ++ghostIndex && ++dfgLowestLayerIteratedIndex) {
#ifndef NDEBUG
            assert(dfgLowestLayerIteratedIndex < numLocalElements);
            IndexVector locAxisIndex(this->getDim());
            dfg_->getLocalVectorIndex(dfgLowestLayerIteratedIndex, locAxisIndex);
            if (locAxisIndex[d] != 0) {
              std::cout << "lowest index = " << dfgLowestLayerIteratedIndex << std::endl;
              std::cout << "locAxisIndex[" << d << "] = " << locAxisIndex << std::endl;
              std::cout << "ghostIndex = " << ghostIndex << " of " << phi_ghost.size() << std::endl;
              std::cout << "offset = " << fullOffsets << " index " << d << std::endl;
            }
            assert(locAxisIndex[d] == 0);
#endif  // NDEBUG
        // compute wrong term to "subtract" again
            const auto wrongNeighborLinearIndex =
                (dfgLowestLayerIteratedIndex - fullOffsets[d] + numLocalElements) %
                numLocalElements;
            assert(wrongNeighborLinearIndex >= 0);
            assert(wrongNeighborLinearIndex < numLocalElements);
            const CombiDataType wrongPhiNeighbor = ElementVector[wrongNeighborLinearIndex];
            // auto wrongDPhi = (dfg_->getDataVector()[ghostIndex] - wrongPhiNeigbor) / h[d];
            const auto dphi = (wrongPhiNeighbor - phi_ghost[ghostIndex]) * oneOverH[d];
            u_dot_dphi[dfgLowestLayerIteratedIndex] += velocity[d] * dphi;
          }
        }
      }
      for (IndexType li = 0; li < dfg_->getNrLocalElements(); ++li) {
        (*phi_)[li] = ElementVector[li] - u_dot_dphi[li] * dt_;
      }
      dfg_->swapDataVector(*phi_);
    }
    stepsTotal_ += nsteps_;

    this->setFinished(true);
  }

  void getFullGrid(FullGrid<CombiDataType>& fg, RankType r, CommunicatorType lcomm,
                   int n = 0) override {
    assert(fg.getLevels() == dfg_->getLevels());

    dfg_->gatherFullGrid(fg, r);
  }

  DistributedFullGrid<CombiDataType>& getDistributedFullGrid(int n = 0) override { return *dfg_; }

  void setZero() override {
    dfg_->setZero();
    std::memset(phi_->data(), 0, phi_->size() * sizeof(CombiDataType));
  }

  ~TaskAdvection() {
    if (dfg_ != nullptr) delete dfg_;
    dfg_ = nullptr;
    if (phi_ != nullptr) delete phi_;
    phi_ = nullptr;
  }

  real getCurrentTime() const override { return stepsTotal_ * dt_; }

  CombiDataType analyticalSolution(const std::vector<real>& coords, int n = 0) const override {
    assert(n == 0);
    auto coordsCopy = coords;
    TestFn f;
    return f(coordsCopy, stepsTotal_ * dt_);
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
  OwningDistributedFullGrid<CombiDataType>* dfg_{};
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
