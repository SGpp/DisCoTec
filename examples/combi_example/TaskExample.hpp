/*
 * TaskExample.hpp
 *
 *  Created on: Sep 25, 2015
 *      Author: heenemo
 */

#ifndef TASKEXAMPLE_HPP_
#define TASKEXAMPLE_HPP_

#include <optional>

#include "fullgrid/DistributedFullGrid.hpp"
#include "task/Task.hpp"

namespace combigrid {

class TaskExample : public Task<> {
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
        stepsTotal_(0) {}

  void init(CommunicatorType lcomm,
            const std::vector<IndexVector>& decomposition = std::vector<IndexVector>()) {
    assert(!initialized_);

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
    dfg_.emplace(makeOwningDistributedFullGrid<CombiDataType>(dim, l, lcomm, this->getBoundary(),
                                                              p_, false, decomposition));

    /* loop over local subgrid and set initial values */
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

  /* this is were the application code kicks in and all the magic happens.
   * do whatever you have to do, but make sure that your application uses
   * only lcomm or a subset of it as communicator.
   * important: don't forget to set the isFinished flag at the end of the computation.
   */
  void run(CommunicatorType lcomm) {
    assert(initialized_);

    std::visit(
        [&](auto& dfg) {
          auto elements = dfg.getData();

          /* pseudo timestepping to demonstrate the behaviour of your typical
           * time-dependent simulation problem. */

          for (size_t step = stepsTotal_; step < stepsTotal_ + nsteps_; ++step) {
            real time = step * dt_;

            for (size_t i = 0; i < dfg.getNrLocalElements(); ++i) {
              IndexType globalLinearIndex = dfg.getGlobalLinearIndex(i);
              std::vector<real> globalCoords(this->getDim());
              dfg.getCoordsGlobal(globalLinearIndex, globalCoords);
              elements[i] = TaskExample::myfunction(globalCoords, time);
            }

            MPI_Barrier(lcomm);
          }
        },
        *dfg_);

    stepsTotal_ += nsteps_;

    this->setFinished(true);
  }

  void getFullGrid(FullGrid<CombiDataType>& fg, RankType r, CommunicatorType lcomm, int n = 0) {
    std::visit([&](auto& dfg) { dfg.gatherFullGrid(fg, r); }, *dfg_);
  }

  DistributedFullGridRef<CombiDataType> getDistributedFullGrid(size_t n = 0) override {
    return toRef<CombiDataType>(*dfg_);
  }

  ConstDistributedFullGridRef<CombiDataType> getDistributedFullGrid(size_t n = 0) const override {
    return toConstRef<CombiDataType>(*dfg_);
  }

  static real myfunction(std::vector<real>& coords, real t) {
    real u = std::cos(M_PI * t);

    for (size_t d = 0; d < coords.size(); ++d) u *= std::cos(2.0 * M_PI * coords[d]);

    return u;
  }

  void setZero() {}

  ~TaskExample() = default;

 protected:
  TaskExample() : initialized_(false), stepsTotal_(1) {}

 private:
  friend class boost::serialization::access;

  // new variables that are set by manager. need to be added to serialize
  real dt_;
  size_t nsteps_ = 0;
  std::vector<int> p_;

  // pure local variables that exist only on the worker processes
  bool initialized_;
  size_t stepsTotal_;
  std::optional<OwningDistributedFullGridVariant<CombiDataType>> dfg_;

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
    ar& boost::serialization::base_object<Task<>>(*this);

    // add our new variables
    ar & dt_;
    ar & nsteps_;
    ar & p_;
  }
};

}  // namespace combigrid

#endif /* TASKEXAMPLE_HPP_ */
