/*
 * SelalibTask.hpp
 *
 *  Created on: Nov 19, 2020
 *      Author: obersteiner
 */

#ifndef SELALIBTASK_HPP_
#define SELALIBTASK_HPP_

#include <stddef.h>

#include <cassert>
#include <functional>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "fault_tolerance/FTUtils.hpp"
#include "fullgrid/DistributedFullGrid.hpp"
#include "hierarchization/DistributedHierarchization.hpp"
#include "task/Task.hpp"
#include "utils/IndexVector.hpp"
#include "utils/LevelVector.hpp"
#include "utils/Types.hpp"

extern "C" {
void sll_s_allocate_collective();
void sll_s_set_communicator_collective(MPI_Fint* mpi_comm);
void sll_s_halt_collective();

void sim_bsl_vp_3d3v_cart_dd_slim_init(void** sim, const char* filename);
void sim_bsl_vp_3d3v_cart_dd_slim_run(void** sim);
void sim_bsl_vp_3d3v_cart_dd_slim_delete(void** sim);
void sim_bsl_vp_3d3v_cart_dd_slim_get_distribution(void** sim, double*& cPtr);
void sim_bsl_vp_3d3v_cart_dd_slim_set_distribution(void** sim, double*& cPtr);
void sim_bsl_vp_3d3v_cart_dd_slim_get_local_size(void** sim, int32_t* cPtr);
void sim_bsl_vp_3d3v_cart_dd_slim_advect_v(void** sim, double* delta_t);
void sim_bsl_vp_3d3v_cart_dd_slim_advect_x(void** sim);
void sim_bsl_vp_3d3v_cart_dd_slim_print_etas(void** sim);
void sim_bsl_vp_3d3v_cart_dd_slim_write_diagnostics_init(void** sim);
void sim_bsl_vp_3d3v_cart_dd_slim_write_diagnostics(void** sim, int32_t* timeStepNumber);
}

namespace combigrid {

class SelalibTask;
inline std::ostream& operator<<(std::ostream& os, const SelalibTask& t);

static constexpr DimType sixdimensions = 6;

class SelalibTask : public combigrid::Task {
 public:
  SelalibTask(LevelVector& l, const std::vector<BoundaryType>& boundary, real coeff,
              LoadModel* loadModel, std::string& path, real dt, size_t nsteps,
              const std::vector<int>& p = std::vector<int>(0),
              FaultCriterion* faultCrit = (new StaticFaults({0, IndexVector(0), IndexVector(0)})))
      : Task(sixdimensions, l, boundary, coeff, loadModel, faultCrit),
        path_(path),
        p_(p),
        localSize_(),
        localDistribution_(nullptr),
        dfg_(nullptr),
        selalibSimPointer_(nullptr),
        initialized_(false),
        diagnosticsInitialized_(false),
        currentNumTimeStepsRun_(0),
        dt_(dt),
        nsteps_(nsteps) {}

  SelalibTask()
      : path_(),
        p_(),
        localSize_(),
        localDistribution_(nullptr),
        dfg_(nullptr),
        initialized_(false) {}

  virtual ~SelalibTask() {
    changeDir(MPI_COMM_SELF, false);
    if (dfg_ != nullptr) {
      delete dfg_;
    }
    sim_bsl_vp_3d3v_cart_dd_slim_delete(simPtrPtr_);
    changeDir(MPI_COMM_SELF, true);
  }

  /**
   * lcomm is the local communicator of the process group.
   */
  void run(CommunicatorType lcomm) {
    currentNumTimeStepsRun_ += nsteps_;
    // only run anything if the coefficient is not 0., i.e. it is not the diagnostics task
    if (coeff_ != 0.) {
      changeDir(lcomm);
      // MASTER_EXCLUSIVE_SECTION{
      //   std::cout << "run " << *this << std::endl;
      // }
      // bool haveResolution = coeff_ == std::numeric_limits<combigrid::real>::max();
      Stats::startEvent("BSL run");
      if (selalibSimPointer_ == nullptr) {
        throw std::runtime_error("selalibSimPointer_ is null");
      }
      sim_bsl_vp_3d3v_cart_dd_slim_run(simPtrPtr_);
      Stats::stopEvent("BSL run");

      int32_t* iPtr = &currentNumTimeStepsRun_;
      sim_bsl_vp_3d3v_cart_dd_slim_write_diagnostics(simPtrPtr_, iPtr);
      changeDir(lcomm, true);
    }
    setFinished(true);
  }

  /**
   * This method changes the folder to the folder of the task
   * lcomm is the local communicator of the process group.
   */
  void changeDir(CommunicatorType lcomm, bool backToBaseFolder = false) {
    if (!backToBaseFolder && chdir(path_.c_str())) {
      printf("could not change to directory %s \n", path_.c_str());
      MPI_Abort(MPI_COMM_WORLD, 1);
    } else if (backToBaseFolder && chdir("..")) {
      printf("could not change to base directory \n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    // it is safer to wait here until all procs are in the right directory
    MPI_Barrier(lcomm);
  }

  /**
   * This method initializes the task
   * lcomm is the local communicator of the process group.
   * decomposition is the spatial decomposition of the component grid
   */
  void init(CommunicatorType lcomm,
            std::vector<IndexVector> decomposition = std::vector<IndexVector>()) {
    static_assert(!reverseOrderingDFGPartitions,
                  "BSL needs this flag to be false, "
                  "change only if the partitioning does not match between dfg and distribution");
    assert(lcomm != MPI_COMM_NULL);

    auto f_lComm = MPI_Comm_c2f(lcomm);
    sll_s_set_communicator_collective(&f_lComm);

    changeDir(lcomm);
    auto paramFilePath = "./param.nml";
    sim_bsl_vp_3d3v_cart_dd_slim_init(simPtrPtr_, paramFilePath);
    MPI_Barrier(lcomm);
    sim_bsl_vp_3d3v_cart_dd_slim_get_local_size(simPtrPtr_, localSize_.data());
    // probably?
    // std::reverse(localSize_.begin(), localSize_.end());
    sim_bsl_vp_3d3v_cart_dd_slim_get_distribution(simPtrPtr_, localDistribution_);
    assert(localDistribution_ != nullptr);
    changeDir(lcomm, true);

    // if coefficient is infty, we are dealing with a single full grid, no need to copy to dfg
    bool haveResolution = coeff_ == std::numeric_limits<combigrid::real>::max();
    if (!haveResolution) {
      dfg_ =
          new DistributedFullGrid<CombiDataType>(getDim(), getLevelVector(), lcomm, getBoundary(),
                                                 localDistribution_, p_, false, decomposition);
    }
    initialized_ = true;
    changeDir(lcomm);
    sim_bsl_vp_3d3v_cart_dd_slim_write_diagnostics_init(simPtrPtr_);
    diagnosticsInitialized_ = true;
    changeDir(lcomm, true);
  }

  /**
   * This method returns the decomposition of the grid of the specified species
   */
  std::vector<IndexVector> getDecomposition(int species) { return dfg_->getDecomposition(); }

  /**
   * This method is used to decide if the execution of the task should fail
   */
  void decideToKill() { assert(false && "not yet implemented"); }

  /**
   * Returns the path of the task
   */
  inline const std::string& getPath() const { return path_; }

  /**
   * This method writes the Selalib grid to the local checkpoint
   */
  void writeLocalDistribution(double* data) { assert(false && "not yet implemented"); }

  /*
   * Gather GENE checkpoint distributed over process group on process
   * with localRootID and convert to FullGrid fg for speciefied species. The actual full grid
   * will only be created on the process with localRootID.
   */
  void getFullGrid(FullGrid<CombiDataType>& fg, RankType lroot, CommunicatorType lcomm,
                   int species = 0) {
    assert(fg.getLevels() == dfg_->getLevels());

    dfg_->gatherFullGrid(fg, lroot);
  }

  /**
   * Returns the distributed full grid of the specified specie
   */
  DistributedFullGrid<CombiDataType>& getDistributedFullGrid(int specie) { return *dfg_; }

  /**
   * @return double* the pointer to the distribution (in Fortran allocated memory)
   */
  inline double* getLocalDistribution() { return localDistribution_; }

  /**
   * sets the dfg content to zero
   */
  void setZero() { dfg_->setZero(); }

  /**
   * Return boolean to indicate whether SelalibTask is initialized.
   */
  inline bool isInitialized() { return initialized_; }

  /**
   * Returns the current time in the simulation.
   */
  real getCurrentTime() const override { return currentNumTimeStepsRun_ * dt_; }

  void setCurrentNumTimeStepsRun(IndexType currentNumTimeStepsRun) {
    currentNumTimeStepsRun_ = currentNumTimeStepsRun;
  }

  real getCurrentNumTimeStepsRun() { return currentNumTimeStepsRun_; }

  // do task-specific postprocessing (by default: nothing)
  void doDiagnostics() override {
    assert(initialized_);
    changeDir(dfg_->getCommunicator(), false);
    // set selalib distribution from dfg
    assert(diagnosticsInitialized_);
    int32_t* iPtr = &currentNumTimeStepsRun_;
    sim_bsl_vp_3d3v_cart_dd_slim_write_diagnostics(simPtrPtr_, iPtr);
    changeDir(dfg_->getCommunicator(), true);
  }

  void receiveDiagnostics() override {}

 private:
  friend class boost::serialization::access;

  // following variables are set in manager and thus need to be included in
  // serialize function
  std::string path_;  // directory in which task should be executed
  std::vector<int> p_;

  // following variables are only accessed in worker and do not need to be
  // serialized
  std::array<int32_t, 6> localSize_;
  double* localDistribution_;
  DistributedFullGrid<CombiDataType>* dfg_;
  void* selalibSimPointer_;
  void** simPtrPtr_ = &selalibSimPointer_;

  bool initialized_;  // indicates if SelalibTask is initialized
  bool diagnosticsInitialized_;

  /*
   * simulation time specific parameters
   */
  int32_t currentNumTimeStepsRun_;  // current number of time steps already run in the simulation
  real dt_;

  /**number of time-steps in between two combinations
   */
  size_t nsteps_;

  // std::chrono::high_resolution_clock::time_point  startTimeIteration_;

  // serialize
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& boost::serialization::base_object<Task>(*this);
    ar& path_;
    ar& p_;
    ar& diagnosticsInitialized_;
    ar& currentNumTimeStepsRun_;
    ar& dt_;
    ar& nsteps_;
    // ar& localDistribution_;
  }
};

inline std::ostream& operator<<(std::ostream& os, const SelalibTask& t) {
  os << "SelalibTask:\n"
     << "\t LevelVector = " << t.getLevelVector() << "\n"
     << "\t Path = " << t.getPath();

  return os;
}

} /* namespace combigrid */

#endif /* SELALIBTASK_HPP_ */
