/*
 * TaskExample.hpp
 *
 *  Created on: Sep 25, 2015
 *      Author: heenemo
 */

#ifndef TASKEXAMPLE_HPP_
#define TASKEXAMPLE_HPP_

#include "fullgrid/DistributedFullGrid.hpp"
#include "task/Task.hpp"

namespace combigrid {

class TaskExample : public Task {
 public:
  /* if the constructor of the base task class is not sufficient we can provide an
   * own implementation. here, we add dt, nsteps, and p as a new parameters.
   */
  TaskExample(const LevelVector& l, const std::vector<BoundaryType>& boundary, real coeff,
              LoadModel* loadModel, real dt, size_t nsteps,
              const std::vector<int>& p = std::vector<int>(0),
              FaultCriterion* faultCrit = (new StaticFaults({0, IndexVector(0), IndexVector(0)})))
      : Task(l, boundary, coeff, loadModel, faultCrit),
        dt_(dt),
        nsteps_(nsteps),
        p_(p),
        initialized_(false),
        stepsTotal_(0),
        dfg_(NULL) {}

  void init(CommunicatorType lcomm,
            const std::vector<IndexVector>& decomposition = std::vector<IndexVector>()) {
    assert(!initialized_);
    assert(dfg_ == NULL);

    int lrank = theMPISystem()->getLocalRank();
    int pgroupNumber = theMPISystem()->getProcessGroupNumber();

    /* create distributed full grid.*/

    DimType dim = this->getDim();
    const LevelVector& l = this->getLevelVector();

    if (lrank == 0) {
      std::cout << "init task " << this->getID() << " with l = " << this->getLevelVector()
                << " and p_ = " << p_ << " on group " << pgroupNumber << std::endl;
    }

    // create local subgrid on each process
    dfg_ = new OwningDistributedFullGrid<CombiDataType>(dim, l, lcomm, this->getBoundary(), p_,
                                                        false, decomposition);

    /* loop over local subgrid and set initial values */
    auto elements = dfg_->getData();

    for (size_t i = 0; i < dfg_->getNrLocalElements(); ++i) {
      IndexType globalLinearIndex = dfg_->getGlobalLinearIndex(i);
      std::vector<real> globalCoords(dim);
      dfg_->getCoordsGlobal(globalLinearIndex, globalCoords);
      elements[i] = TaskExample::myfunction(globalCoords, 0.0);
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

    auto elements = dfg_->getData();
    // TODO if your Example uses another data structure, you need to copy
    // the data from elements to that data structure

    /* pseudo timestepping to demonstrate the behaviour of your typical
     * time-dependent simulation problem. */
    // TODO replace by your time stepping algorithm

    for (size_t step = stepsTotal_; step < stepsTotal_ + nsteps_; ++step) {
      real time = step * dt_;

      for (size_t i = 0; i < dfg_->getNrLocalElements(); ++i) {
        IndexType globalLinearIndex = dfg_->getGlobalLinearIndex(i);
        std::vector<real> globalCoords(this->getDim());
        dfg_->getCoordsGlobal(globalLinearIndex, globalCoords);
        elements[i] = TaskExample::myfunction(globalCoords, time);
      }

      MPI_Barrier(lcomm);
    }

    stepsTotal_ += nsteps_;

    // TODO if your Example uses another data structure, you need to copy
    // the data from that data structure to elements/dfg_

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

  DistributedFullGrid<CombiDataType>& getDistributedFullGrid(size_t n = 0) override { return *dfg_; }

  const DistributedFullGrid<CombiDataType>& getDistributedFullGrid(size_t n = 0) const override { return *dfg_; }

  static real myfunction(std::vector<real>& coords, real t) {
    real u = std::cos(M_PI * t);

    for (size_t d = 0; d < coords.size(); ++d) u *= std::cos(2.0 * M_PI * coords[d]);

    return u;
  }

  void setZero() {}

  ~TaskExample() {
    if (dfg_ != NULL) delete dfg_;
  }

 protected:
  /* if there are local variables that have to be initialized at construction
   * you have to do it here. the worker processes will create the task using
   * this constructor before overwriting the variables that are set by the
   * manager. here we need to set the initialized variable to make sure it is
   * set to false. */
  TaskExample() : initialized_(false), stepsTotal_(1), dfg_(NULL) {}

 private:
  friend class boost::serialization::access;

  // new variables that are set by manager. need to be added to serialize
  real dt_;
  size_t nsteps_ = 0;
  std::vector<int> p_;

  // pure local variables that exist only on the worker processes
  bool initialized_;
  size_t stepsTotal_;
  OwningDistributedFullGrid<CombiDataType>* dfg_;

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

#endif /* TASKEXAMPLE_HPP_ */
