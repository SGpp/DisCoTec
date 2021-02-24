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

#include "sgpp/distributedcombigrid/fault_tolerance/FTUtils.hpp"
#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGrid.hpp"
#include "sgpp/distributedcombigrid/task/Task.hpp"
#include "sgpp/distributedcombigrid/utils/IndexVector.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"

extern "C" {
//  void __sll_m_collective_MOD_sll_s_boot_collective(int32_t *mpi_mode);
void sll_s_allocate_collective();
void sll_s_set_communicator_collective(MPI_Fint* mpi_comm);
void sll_s_halt_collective();

  void sim_bsl_vp_3d3v_cart_dd_slim_movingB_init(void** sim, const char *filename);
  void sim_bsl_vp_3d3v_cart_dd_slim_movingB_run(void** sim);
  void sim_bsl_vp_3d3v_cart_dd_slim_movingB_delete(void** sim);
  void sim_bsl_vp_3d3v_cart_dd_slim_movingB_get_distribution(void** sim, double *& cPtr);
  void sim_bsl_vp_3d3v_cart_dd_slim_movingB_set_distribution(void** sim, double *& cPtr);
  void sim_bsl_vp_3d3v_cart_dd_slim_movingB_get_local_size(void** sim, int32_t *cPtr);
  void sim_bsl_vp_3d3v_cart_dd_slim_movingB_advect_v(void** sim, double *delta_t);
  void sim_bsl_vp_3d3v_cart_dd_slim_movingB_advect_x(void** sim, double *delta_t);
}

namespace combigrid {

class SelalibTask;
inline std::ostream& operator<<(std::ostream& os, const SelalibTask& t);

class SelalibTask : public combigrid::Task {
 public:
  SelalibTask(DimType dim, LevelVector& l, std::vector<bool>& boundary, real coeff,
              LoadModel* loadModel, std::string& path, real dt, real combitime, size_t nsteps,
              IndexVector p = IndexVector(0),
              FaultCriterion* faultCrit = (new StaticFaults({0, IndexVector(0), IndexVector(0)})))
      : Task(dim, l, boundary, coeff, loadModel, faultCrit),
        path_(path),
        p_(p),
        localSize_(),
        localDistribution_(nullptr),
        dfg_(nullptr),
        selalibSimPointer_(nullptr),
        initialized_(false),
        currentTime_(0.0),
        currentTimestep_(0),
        dt_(dt),
        combitime_(combitime),
        nsteps_(nsteps) {}

  SelalibTask()
      : path_(),
        p_(),
        localSize_(),
        localDistribution_(nullptr),
        dfg_(nullptr),
        initialized_(false) {}

  virtual ~SelalibTask() {
    if (dfg_ != nullptr) {
      delete dfg_;
    }
    sim_bsl_vp_3d3v_cart_dd_slim_movingB_delete(simPtrPtr_);
  }

  /**
   * lcomm is the local communicator of the process group.
   */
  void run(CommunicatorType lcomm) {
    // TODO set distribution from dfg
    changeDir(lcomm);
    MASTER_EXCLUSIVE_SECTION{
      std::cout << "run " << *this << std::endl;
    }
    setLocalDistributionFromDFG();
    sim_bsl_vp_3d3v_cart_dd_slim_movingB_run(simPtrPtr_);
    setDFGfromLocalDistribution();
    changeDir(lcomm, true);
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
    assert(lcomm != MPI_COMM_NULL);
    dfg_ = new DistributedFullGrid<CombiDataType>(getDim(), getLevelVector(), lcomm, getBoundary(),
                                                  p_);
    auto f_lComm = MPI_Comm_c2f(lcomm);
    sll_s_set_communicator_collective(&f_lComm);
    changeDir(lcomm);
    auto paramFilePath = "./param.nml";
    sim_bsl_vp_3d3v_cart_dd_slim_movingB_init(simPtrPtr_, paramFilePath);
    MPI_Barrier(lcomm);
    sim_bsl_vp_3d3v_cart_dd_slim_movingB_get_local_size(simPtrPtr_, localSize_.data());
    // probably?
    // std::reverse(localSize_.begin(), localSize_.end());
    sim_bsl_vp_3d3v_cart_dd_slim_movingB_get_distribution(simPtrPtr_, localDistribution_);
    assert(localDistribution_ != nullptr);
    changeDir(lcomm, true);
    setDFGfromLocalDistribution();
    initialized_ = true;
    MASTER_EXCLUSIVE_SECTION{
      std::cout << "initialized " << *this << std::endl;
    }
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
   * Returns the time that is simulated between combinations.
   * This is only used in case we do not want to use a fixed number of timesteps
   * but a fixed period of time between combinations for each component grids.
   */
  real getCombiTime() { return combitime_; }

  /**
   * Returns the current time in the simulation. This is uses to update the time in BSL after
   * restart.
   */
  real getCurrentTime() const override { return currentTime_; }

  /**
   * Sets the current time in the simulation. This is uses to update the time in BSL after restart.
   */
  void setCurrentTime(real currentTime) { currentTime_ = currentTime; }
  /**
   * Returns the current timestep in the simulation. This is uses to update the timestep in BSL
   * after restart.
   */
  real getCurrentTimestep() { return currentTimestep_; }
  /**
   * Sets the current timestep in the simulation. This is uses to update the timestep in BSL after
   * restart.
   */
  void setCurrentTimestep(real currentTimestep) { currentTimestep_ = currentTimestep; }

 private:
  friend class boost::serialization::access;

  // following variables are set in manager and thus need to be included in
  // serialize function
  std::string path_;  // directory in which task should be executed
  IndexVector p_;

  // following variables are only accessed in worker and do not need to be
  // serialized
  std::array<int32_t, 6> localSize_;
  double* localDistribution_;
  DistributedFullGrid<CombiDataType>* dfg_;
  void* selalibSimPointer_;
  void** simPtrPtr_ = &selalibSimPointer_;

  bool initialized_;  // indicates if SelalibTask is initialized

  /*
   * simulation time specific parameters
   */
  real currentTime_;      // current time in the simulation
  real currentTimestep_;  // curent time step length in the simulation
  real dt_;
  real combitime_;  // simulation time interval between combinations

  /**number of time-steps in between two combinations (is set very large in case combitime should be
   * used); this requires equal time-steps for every component grid
   */
  size_t nsteps_;

  // std::chrono::high_resolution_clock::time_point  startTimeIteration_;

  void setDFGfromLocalDistribution() {
    auto& localDFGSize = dfg_->getLocalSizes();
    assert(dim_ == 6);
    // std::cout << localDFGSize << " " << std::vector<int>(localSize_.begin(), localSize_.end()) <<
    // std::endl;
    auto localSizeCopy = localSize_;
    sim_bsl_vp_3d3v_cart_dd_slim_movingB_get_local_size(simPtrPtr_, localSizeCopy.data());
    for (DimType d = 0; d < dim_; ++d) {
      assert(localDFGSize[d] == (localSize_[d] + 1) || localDFGSize[d] == localSize_[d]);
      assert(localSizeCopy[d] == localSize_[d]);
    }

    auto& offsets = dfg_->getLocalOffsets();
    assert(offsets[0] == 1);

    auto localDistributionIterator = localDistribution_;
    int32_t bufferSize =
        std::accumulate(localSize_.begin(), localSize_.end(), 1, std::multiplies<int32_t>());
    // assignment, leaving out the uppermost layer in each dimension, if it is not part of
    // Selalib's local distribution
    // contiguous access by Fortran ordering
    for (int i = 0; i < localSize_[5]; ++i) {
      auto offset_i = i * offsets[5];
      for (int j = 0; j < localSize_[4]; ++j) {
        auto offset_j = offset_i + j * offsets[4];
        for (int k = 0; k < localSize_[3]; ++k) {
          auto offset_k = offset_j + k * offsets[3];
          for (int l = 0; l < localSize_[2]; ++l) {
            auto offset_l = offset_k + l * offsets[2];
            for (int m = 0; m < localSize_[1]; ++m) {
              auto offset_m = offset_l + m * offsets[1];
              for (int n = 0; n < localSize_[0]; ++n) {
                auto fgIndex = offset_m + n;
                dfg_->getElementVector()[fgIndex] = *(localDistributionIterator);
                ++localDistributionIterator;
                assert((localDistributionIterator - localDistribution_) <= bufferSize);
              }
            }
          }
        }
      }
    }
    for (DimType d = 0; d < dim_; ++d) {
      dfg_->writeLowerBoundaryToUpperBoundary(d);
    }
  }

  void setLocalDistributionFromDFG() {
    // cf setDFGfromLocalDistribution
    auto& offsets = dfg_->getLocalOffsets();
    auto localDistributionIterator = localDistribution_;
    int32_t bufferSize =
        std::accumulate(localSize_.begin(), localSize_.end(), 1, std::multiplies<int32_t>());


    for (int i = 0; i < localSize_[5]; ++i) {
      auto offset_i = i * offsets[5];
      for (int j = 0; j < localSize_[4]; ++j) {
        auto offset_j = offset_i + j * offsets[4];
        for (int k = 0; k < localSize_[3]; ++k) {
          auto offset_k = offset_j + k * offsets[3];
          for (int l = 0; l < localSize_[2]; ++l) {
            auto offset_l = offset_k + l * offsets[2];
            for (int m = 0; m < localSize_[1]; ++m) {
              auto offset_m = offset_l + m * offsets[1];
              for (int n = 0; n < localSize_[0]; ++n) {
                auto fgIndex = offset_m + n;
                *(localDistributionIterator) = dfg_->getElementVector()[fgIndex];
                ++localDistributionIterator;
                assert((localDistributionIterator - localDistribution_) <= bufferSize);
              }
            }
          }
        }
      }
    }
  }

  // serialize
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& boost::serialization::base_object<Task>(*this);
    ar& path_;
    ar& dt_;
    ar& combitime_;
    ar& nsteps_;
    // ar& localDistribution_;
    ar& p_;
    ar& currentTime_;
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
