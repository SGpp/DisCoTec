/*
 * TaskExample.hpp
 *
 *  Created on: Sep 25, 2015
 *      Author: heenemo
 */

#ifndef TASKEXAMPLE_HPP_
#define TASKEXAMPLE_HPP_

#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGrid.hpp"
#include "sgpp/distributedcombigrid/task/Task.hpp"

namespace combigrid {

class TaskExample: public Task {

 public:
  /* if the constructor of the base task class is not sufficient we can provide an
   * own implementation. here, we add dt, nsteps, and p as a new parameters.
   */
  TaskExample(DimType dim, LevelVector& l, std::vector<bool>& boundary,
              real coeff, LoadModel* loadModel, real dt,
              size_t nsteps, IndexVector p = IndexVector(0),FaultCriterion *faultCrit = (new StaticFaults({0,IndexVector(0),IndexVector(0)})) ) :
    Task(dim, l, boundary, coeff, loadModel, faultCrit), dt_(dt), nsteps_(
      nsteps), p_(p), initialized_(false), stepsTotal_(0), dfg_(NULL) {
  }

    void init(CommunicatorType lcomm, std::vector<IndexVector> decomposition = std::vector<IndexVector>()){
    assert(!initialized_);
    assert(dfg_ == NULL);

    int lrank;
    MPI_Comm_rank(lcomm, &lrank);

    /* create distributed full grid. we try to find a balanced ratio between
     * the number of grid points and the number of processes per dimension
     * by this very simple algorithm. to keep things simple we require powers
     * of two for the number of processes here. */
    int np;
    MPI_Comm_size(lcomm, &np);

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

        for (DimType k = 0; k < dim; ++k)
          prod_p *= p[k];
      }
    } else {
      p = p_;
    }

    if (lrank == 0) {
      std::cout << "init task " << this->getID() << " with l = "
                << this->getLevelVector() << " and p = " << p << std::endl;
    }

    // create local subgrid on each process
    dfg_ = new DistributedFullGrid<CombiDataType>(dim, l, lcomm,
        this->getBoundary(), p);

    /* loop over local subgrid and set initial values */
    std::vector<CombiDataType>& elements = dfg_->getElementVector();

    phi_.resize(elements.size());
    assert(phi_.size() == elements.size());

    for (IndexType li = 0; li < elements.size(); ++li) {
      std::vector<double> coords(this->getDim());
      dfg_->getCoordsGlobal(li, coords);

      double exponent = 0;
      for (DimType d = 0; d < this->getDim(); ++d) {
        exponent -= std::pow(coords.at(d) - 0.5, 2);
      }
      dfg_->getElementVector().at(li) = std::exp(exponent*100.0) * 2;
    }

    initialized_ = true;
  }


  /* this is were the application code kicks in and all the magic happens.
   * do whatever you have to do, but make sure that your application uses
   * only lcomm or a subset of it as communicator.
   * important: don't forget to set the isFinished flag at the end of the computation.
   */
  void run(CommunicatorType lcomm) {
    std::cout << "run my task\n";
    assert(initialized_);

    int lrank;
    MPI_Comm_rank(lcomm, &lrank);

    /* pseudo timestepping to demonstrate the behaviour of your typical
     * time-dependent simulation problem. */
     std::vector<CombiDataType> u(this->getDim(), 0.0000001);

     // gradient of phi
     std::vector<CombiDataType> dphi(this->getDim());

     std::vector<IndexType> l(this->getDim());
     std::vector<double> h(this->getDim());

     for (int i = 0; i < this->getDim(); i++){
       l[i] = dfg_->length(i);
       h[i] = 1.0 / (double)l[i];
     }

     for (size_t i = 0; i < nsteps_; ++i) {
       phi_.swap(dfg_->getElementVector());

       for (IndexType li = 0; li < dfg_->getElementVector().size(); ++li) {
         IndexVector ai(this->getDim());
         dfg_->getGlobalVectorIndex(li, ai);

         //neighbour
         std::vector<IndexVector> ni(this->getDim(), ai);
         std::vector<IndexType> lni(this->getDim());

         CombiDataType u_dot_dphi = 0;

         for(int j = 0; j < this->getDim(); j++){
           ni[j][j] = (l[j] + ni[j][j] - 1) % l[j];
           lni[j] = dfg_->getGlobalLinearIndex(ni[j]);
         }

         for(int j = 0; j < this->getDim(); j++){
           //calculate gradient of phi with backward differential quotient
           dphi.at(j) = (phi_.at(li) - phi_.at(lni.at(j))) / h.at(j);

           u_dot_dphi += u[j] * dphi[j];
         }

         dfg_->getData()[li] = phi_[li] - u_dot_dphi * dt_;
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
  void getFullGrid(FullGrid<CombiDataType>& fg, RankType r,
                   CommunicatorType lcomm, int n = 0) {
    assert(fg.getLevels() == dfg_->getLevels());

    dfg_->gatherFullGrid(fg, r);
  }

  DistributedFullGrid<CombiDataType>& getDistributedFullGrid(int n = 0) {
    return *dfg_;
  }


  void setZero(){

  }

 protected:
  /* if there are local variables that have to be initialized at construction
   * you have to do it here. the worker processes will create the task using
   * this constructor before overwriting the variables that are set by the
   * manager. here we need to set the initialized variable to make sure it is
   * set to false. */
  TaskExample() :
    initialized_(false), stepsTotal_(1), dfg_(NULL) {
  }

  ~TaskExample() {
    if (dfg_ != NULL)
      delete dfg_;
  }

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
  std::vector<CombiDataType> phi_;

  /**
   * The serialize function has to be extended by the new member variables.
   * However this concerns only member variables that need to be exchanged
   * between manager and workers. We do not need to add "local" member variables
   * that are only needed on either manager or worker processes.
   * For serialization of the parent class members, the class must be
   * registered with the BOOST_CLASS_EXPORT macro.
   */
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    // handles serialization of base class
    ar& boost::serialization::base_object<Task>(*this);

    // add our new variables
    ar& dt_;
    ar& nsteps_;
    ar& p_;
  }
};

} // namespace combigrid

#endif /* TASKEXAMPLE_HPP_ */
