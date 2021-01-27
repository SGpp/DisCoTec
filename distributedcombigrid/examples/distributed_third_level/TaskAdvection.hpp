#pragma once

#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGrid.hpp"
#include "sgpp/distributedcombigrid/task/Task.hpp"

namespace combigrid {

// exact solution
class TestFn {
 public:
  // function value
  double operator()(std::vector<double>& coords, double t) {
    double exponent = 0;
    for (DimType d = 0; d < coords.size(); ++d) {
      coords[d] = std::fmod(1.0 + std::fmod(coords[d] - t, 1.0), 1.0);
      exponent -= std::pow(coords[d] - 0.5, 2);
    }
    return std::exp(exponent*100.0) * 2;
  }
};

class TaskAdvection : public Task {
 public:
  /* if the constructor of the base task class is not sufficient we can provide an
   * own implementation. here, we add dt, nsteps, and p as a new parameters.
   */
  TaskAdvection(DimType dim, LevelVector& l, std::vector<bool>& boundary, real coeff,
                LoadModel* loadModel, real dt, size_t nsteps, IndexVector p = IndexVector(0),
                FaultCriterion* faultCrit = (new StaticFaults({0, IndexVector(0), IndexVector(0)})))
      : Task(dim, l, boundary, coeff, loadModel, faultCrit),
        dt_(dt),
        nsteps_(nsteps),
        p_(p),
        initialized_(false),
        stepsTotal_(0),
        dfg_(NULL) {}

  void init(CommunicatorType lcomm,
            std::vector<IndexVector> decomposition = std::vector<IndexVector>()) {
    assert(!initialized_);
    assert(dfg_ == NULL);

    auto lrank = theMPISystem()->getLocalRank();
    int np;
    MPI_Comm_size(lcomm, &np);

    /* create distributed full grid. we try to find a balanced ratio between
     * the number of grid points and the number of processes per dimension
     * by this very simple algorithm. to keep things simple we require powers
     * of two for the number of processes here. */

    // check if power of two
    if (!((np > 0) && ((np & (~np + 1)) == np)))
      assert(false && "number of processes not power of two");

    DimType dim = this->getDim();
    IndexVector p(dim, 1);
    const LevelVector& l = this->getLevelVector();

    if (p_.size() == 0) {
      // compute domain decomposition
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
      std::cout << "init task " << this->getID() << " with l = " << this->getLevelVector()
                << " and p = " << p << std::endl;
    }

    // create local subgrid on each process
    dfg_ = new DistributedFullGrid<CombiDataType>(dim, l, lcomm, this->getBoundary(), p);
    phi_ = new DistributedFullGrid<CombiDataType>(dim, l, lcomm, this->getBoundary(), p);

    for (IndexType li = 0; li < dfg_->getNrLocalElements(); ++li) {
      std::vector<double> coords(this->getDim());
      dfg_->getCoordsLocal(li, coords);

      double exponent = 0;
      for (DimType d = 0; d < this->getDim(); ++d) {
        exponent -= std::pow(coords[d] - 0.5, 2);
      }
      dfg_->getData()[li] = std::exp(exponent * 100.0) * 2;
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

    std::vector<CombiDataType> u(this->getDim(), 1);

    // gradient of phi
    std::vector<CombiDataType> dphi(this->getDim());

    std::vector<double> h(this->getDim());

    auto l = dfg_->getLocalSizes();
    for (unsigned int i = 0; i < this->getDim(); i++) {
      h[i] = 1.0 / (double)l[i];
    }

    for (size_t i = 0; i < nsteps_; ++i) {
      phi_->getElementVector().swap(dfg_->getElementVector());

      auto u_dot_dphi = std::vector<CombiDataType> (dfg_->getNrLocalElements(), 0.);

      for (unsigned int d = 0; d < this->getDim(); ++d) {
        // to update the values in the "lowest" layer, we need the ghost values from the lower neighbor
        IndexVector subarrayExtents;
        auto phi_ghost = phi_->exchangeGhostLayerUpward(d, subarrayExtents);
        // std::cout << "phi_ghost " << phi_ghost << std::endl;
        int offset = 1;
        IndexVector offsets (this->getDim());
        for (DimType d_j = 0; d_j < dim_; ++d_j) {
          offsets[d_j] = offset;
          offset *= subarrayExtents[d_j];
        }

        for (IndexType li = 0; li < dfg_->getNrLocalElements(); ++li) {
          // calculate local axis index of backward neighbor
          IndexVector locAxisIndex(this->getDim());
          dfg_->getLocalVectorIndex(li, locAxisIndex);
          //TODO can be unrolled into ghost and other part, avoiding if-statement
          CombiDataType phi_neighbor;
          if (locAxisIndex[d] == 0){
            // if we are in the lowest layer in d,
            // make sure we are not on the lowest global layer
            IndexVector globAxisIndex(this->getDim());
            dfg_->getGlobalVectorIndex(locAxisIndex, globAxisIndex);
            if (globAxisIndex[d] == 0){
              assert(phi_ghost.size()==0);
              continue;
              // remove for periodic BC
            }
            // then use values from boundary exchange
            IndexType gli = 0;
            for (DimType d_j = 0; d_j < this->getDim(); ++d_j) {
              gli = gli + offsets[d_j] * locAxisIndex[d_j];
            }
            assert(gli > -1);
            assert(gli < phi_ghost.size());
            phi_neighbor = phi_ghost[gli];
          } else{
            --locAxisIndex[d];
            IndexType lni = dfg_->getLocalLinearIndex(locAxisIndex);
            phi_neighbor = phi_->getElementVector()[lni];
          }
          // calculate gradient of phi with backward differential quotient
          auto dphi = (phi_->getElementVector()[li] - phi_neighbor) / h[d];

          u_dot_dphi[li] += u[d] * dphi;
        }
      }
      for (IndexType li = 0; li < dfg_->getNrLocalElements(); ++li) {
        dfg_->getData()[li] = phi_->getElementVector()[li] - u_dot_dphi[li] * dt_;
      }
      MPI_Barrier(lcomm);
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
    phi_->setZero();
  }

  ~TaskAdvection() {
    if (dfg_ != NULL) delete dfg_;
    if (phi_ != NULL) delete phi_;
  }

 protected:
  /* if there are local variables that have to be initialized at construction
   * you have to do it here. the worker processes will create the task using
   * this constructor before overwriting the variables that are set by the
   * manager. here we need to set the initialized variable to make sure it is
   * set to false. */
  TaskAdvection() : initialized_(false), stepsTotal_(1), dfg_(nullptr), phi_(nullptr) {}

 private:
  friend class boost::serialization::access;

  // new variables that are set by manager. need to be added to serialize
  real dt_;
  size_t nsteps_;
  IndexVector p_;

  // pure local variables that exist only on the worker processes
  bool initialized_;
  size_t stepsTotal_;
  DistributedFullGrid<CombiDataType>* dfg_;
  DistributedFullGrid<CombiDataType>* phi_;

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

}  // namespace combigrid
