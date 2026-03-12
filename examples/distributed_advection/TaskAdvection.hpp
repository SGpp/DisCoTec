#pragma once

#include <optional>

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
    return std::exp(exponent);
  }

 private:
  static constexpr double sigmaSquaredInv_ = 1. / ((1. / 3.) * (1. / 3.));
};

class TaskAdvection : public Task<> {
 public:
  TaskAdvection(const LevelVector& l, const std::vector<BoundaryType>& boundary, real coeff,
                LoadModel* loadModel, real dt, size_t nsteps,
                const std::vector<int>& p = std::vector<int>(0),
                FaultCriterion* faultCrit = (new StaticFaults({0, IndexVector(0), IndexVector(0)})))
      : Task(l, boundary, coeff, loadModel, faultCrit),
        dt_(dt),
        nsteps_(nsteps),
        p_(std::move(p)),
        initialized_(false),
        stepsTotal_(0) {
    for ([[maybe_unused]] const auto& b : boundary) {
      assert(b == 1);
    }
  }

  void init(CommunicatorType lcomm,
            const std::vector<IndexVector>& decomposition = std::vector<IndexVector>()) override {
    assert(!initialized_);

    DimType dim = this->getDim();
    const LevelVector& l = this->getLevelVector();

    dfg_.emplace(makeOwningDistributedFullGrid<CombiDataType>(dim, l, lcomm, this->getBoundary(),
                                                              p_, false, decomposition));

    std::visit(
        [&](auto& dfg) {
          if (phi_ == nullptr) {
            phi_ = new std::vector<CombiDataType>(dfg.getNrLocalElements());
          }

          const auto& h = dfg.getGridSpacing();
          auto sumOneOverH = 0.;
          for (const auto& h_x : h) {
            sumOneOverH += 1. / h_x;
          }
          if (!(dt_ * sumOneOverH < 1.)) {
            throw std::runtime_error("CFL condition not satisfied!");
          }

          TestFn f;
          static thread_local std::vector<double> coords(this->getDim());
#pragma omp parallel for schedule(static) default(none) shared(f, dfg)
          for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
            dfg.getCoordsLocal(li, coords);
            dfg.getData()[li] = f(coords, 0.);
          }
        },
        *dfg_);

    initialized_ = true;
  }

  void run(CommunicatorType lcomm) override {
    assert(initialized_);

    std::visit(
        [&](auto& dfg) {
          const auto numLocalElements = dfg.getNrLocalElements();

          const std::vector<CombiDataType> velocity(this->getDim(), 1);

          const std::vector<double> oneOverH = dfg.getInverseGridSpacing();
          const auto& fullOffsets = dfg.getLocalOffsets();

          phi_->resize(numLocalElements);

          for (size_t i = 0; i < nsteps_; ++i) {
            std::memset(phi_->data(), 0, phi_->size() * sizeof(CombiDataType));
            auto& u_dot_dphi = *phi_;
            auto const ElementVector = dfg.getData();
            for (unsigned int d = 0; d < this->getDim(); ++d) {
              static std::vector<int> subarrayExtents;
              std::vector<CombiDataType> phi_ghost{};
              MPI_Request recvRequest;
              dfg.exchangeGhostLayerUpward(d, subarrayExtents, phi_ghost, &recvRequest);
              auto fullOffsetsInThisDimension = fullOffsets[d];
              {
#pragma omp parallel for schedule(static) default(none)           \
    firstprivate(d, numLocalElements, fullOffsetsInThisDimension) \
    shared(u_dot_dphi, ElementVector, oneOverH, velocity)
                for (IndexType li = 0; li < fullOffsetsInThisDimension; ++li) {
#ifndef NDEBUG
                  IndexVector locAxisIndex(this->getDim());
                  dfg.getLocalVectorIndex(li, locAxisIndex);
                  if (locAxisIndex[d] > 0) {
                    --locAxisIndex[d];
                    IndexType lni = dfg.getLocalLinearIndex(locAxisIndex);
                    assert(lni == (li - fullOffsetsInThisDimension));
                  }
#endif
                  const auto neighborLinearIndex =
                      (li - fullOffsetsInThisDimension + numLocalElements);
                  assert(neighborLinearIndex >= 0);
                  assert(neighborLinearIndex < numLocalElements);
                  CombiDataType phi_neighbor = ElementVector[neighborLinearIndex];
                  auto dphi = (ElementVector[li] - phi_neighbor) * oneOverH[d];
                  u_dot_dphi[li] += velocity[d] * dphi;
                }
#pragma omp parallel for schedule(static) default(none)           \
    firstprivate(d, numLocalElements, fullOffsetsInThisDimension) \
    shared(u_dot_dphi, ElementVector, oneOverH, velocity)
                for (IndexType li = fullOffsetsInThisDimension; li < numLocalElements; ++li) {
#ifndef NDEBUG
                  IndexVector locAxisIndex(this->getDim());
                  dfg.getLocalVectorIndex(li, locAxisIndex);
                  if (locAxisIndex[d] > 0) {
                    --locAxisIndex[d];
                    IndexType lni = dfg.getLocalLinearIndex(locAxisIndex);
                    assert(lni == (li - fullOffsetsInThisDimension));
                  }
#endif
                  const auto neighborLinearIndex = (li - fullOffsetsInThisDimension);
                  assert(neighborLinearIndex >= 0);
                  assert(neighborLinearIndex < numLocalElements);
                  CombiDataType phi_neighbor = ElementVector[neighborLinearIndex];
                  auto dphi = (ElementVector[li] - phi_neighbor) * oneOverH[d];
                  u_dot_dphi[li] += velocity[d] * dphi;
                }
              }
              assert(static_cast<size_t>(numLocalElements) / dfg.getLocalSizes()[d] ==
                     phi_ghost.size());
              const auto& stride = fullOffsetsInThisDimension;
              const IndexType jump = stride * dfg.getLocalSizes()[d];
              const IndexType numberOfPolesHigherDimensions = numLocalElements / jump;
              MPI_Wait(&recvRequest, MPI_STATUS_IGNORE);
              {
#pragma omp parallel for schedule(static) default(none)                            \
    firstprivate(d, numLocalElements, stride, jump, numberOfPolesHigherDimensions, \
                     fullOffsetsInThisDimension)                                   \
    shared(u_dot_dphi, ElementVector, oneOverH, phi_ghost, velocity, std::cout)
                for (IndexType nLower = 0; nLower < stride; ++nLower) {
                  IndexType dfgLowestLayerIteratedIndex = nLower;
                  IndexType ghostIndex = nLower;
#ifndef NDEBUG
                  assert(dfgLowestLayerIteratedIndex < numLocalElements);
                  IndexVector locAxisIndex(this->getDim());
                  dfg.getLocalVectorIndex(dfgLowestLayerIteratedIndex, locAxisIndex);
                  assert(locAxisIndex[d] == 0);
#endif
                  const auto wrongNeighborLinearIndex =
                      (dfgLowestLayerIteratedIndex - fullOffsetsInThisDimension + numLocalElements);
                  assert(wrongNeighborLinearIndex >= 0);
                  assert(wrongNeighborLinearIndex < numLocalElements);
                  const CombiDataType wrongPhiNeighbor = ElementVector[wrongNeighborLinearIndex];
                  const auto dphi = (wrongPhiNeighbor - phi_ghost[ghostIndex]) * oneOverH[d];
                  u_dot_dphi[dfgLowestLayerIteratedIndex] += velocity[d] * dphi;
                }

#pragma omp parallel for collapse(2) schedule(static) default(none)                \
    firstprivate(d, numLocalElements, stride, jump, numberOfPolesHigherDimensions, \
                     fullOffsetsInThisDimension)                                   \
    shared(u_dot_dphi, ElementVector, oneOverH, phi_ghost, velocity, std::cout)
                for (IndexType nHigher = 1; nHigher < numberOfPolesHigherDimensions; ++nHigher) {
                  for (IndexType nLower = 0; nLower < stride; ++nLower) {
                    IndexType dfgLowestLayerIteratedIndex = nHigher * jump + nLower;
                    IndexType ghostIndex = nLower + nHigher * stride;
#ifndef NDEBUG
                    assert(dfgLowestLayerIteratedIndex < numLocalElements);
                    IndexVector locAxisIndex(this->getDim());
                    dfg.getLocalVectorIndex(dfgLowestLayerIteratedIndex, locAxisIndex);
                    assert(locAxisIndex[d] == 0);
#endif
                    const auto wrongNeighborLinearIndex =
                        (dfgLowestLayerIteratedIndex - fullOffsetsInThisDimension);
                    assert(wrongNeighborLinearIndex >= 0);
                    assert(wrongNeighborLinearIndex < numLocalElements);
                    const CombiDataType wrongPhiNeighbor = ElementVector[wrongNeighborLinearIndex];
                    const auto dphi = (wrongPhiNeighbor - phi_ghost[ghostIndex]) * oneOverH[d];
                    u_dot_dphi[dfgLowestLayerIteratedIndex] += velocity[d] * dphi;
                  }
                }
              }
            }
#pragma omp parallel for simd schedule(simd : static) default(none) \
    shared(ElementVector, u_dot_dphi) firstprivate(dt_, numLocalElements)
            for (IndexType li = 0; li < numLocalElements; ++li) {
              (*phi_)[li] = ElementVector[li] - u_dot_dphi[li] * dt_;
            }
            dfg.swapDataVector(*phi_);
          }
        },
        *dfg_);
    stepsTotal_ += nsteps_;

    this->setFinished(true);
  }

  void getFullGrid(FullGrid<CombiDataType>& fg, RankType r, CommunicatorType lcomm,
                   int n = 0) override {
    std::visit([&](auto& dfg) { dfg.gatherFullGrid(fg, r); }, *dfg_);
  }

  DistributedFullGridRef<CombiDataType> getDistributedFullGrid(size_t n = 0) override {
    return toRef<CombiDataType>(*dfg_);
  }

  void setZero() override {
    std::visit([&](auto& dfg) { dfg.setZero(); }, *dfg_);
    std::memset(phi_->data(), 0, phi_->size() * sizeof(CombiDataType));
  }

  ~TaskAdvection() {
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
  TaskAdvection() : initialized_(false), stepsTotal_(0) {}

 private:
  friend class boost::serialization::access;

  real dt_;
  size_t nsteps_;
  std::vector<int> p_;

  bool initialized_;
  size_t stepsTotal_;
  std::optional<OwningDistributedFullGridVariant<CombiDataType>> dfg_;
  static std::vector<CombiDataType>* phi_;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& boost::serialization::base_object<Task<>>(*this);

    ar & dt_;
    ar & nsteps_;
    ar & p_;
  }
};

std::vector<CombiDataType>* TaskAdvection::phi_ = nullptr;

}  // namespace combigrid
