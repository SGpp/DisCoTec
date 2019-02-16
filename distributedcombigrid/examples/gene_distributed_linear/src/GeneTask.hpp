/*
 * GeneTask.hpp
 *
 *  Created on: Jul 10, 2014
 *      Author: heenemo
 */

#ifndef GENETASK_HPP_
#define GENETASK_HPP_

#include <stddef.h>
#include <cassert>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGrid.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LoadModel.hpp"
#include "sgpp/distributedcombigrid/task/Task.hpp"
#include "sgpp/distributedcombigrid/utils/IndexVector.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "GeneLocalCheckpoint.hpp"
#include "sgpp/distributedcombigrid/fault_tolerance/FTUtils.hpp"
#include "sgpp/distributedcombigrid/fullgrid/MultiArray.hpp"

namespace combigrid {

class GeneTask: public combigrid::Task {
public:
  GeneTask( DimType dim, LevelVector& l, std::vector<bool>& boundary, real coeff,
            LoadModel* loadModel, std::string& path, real dt, real combitime, size_t nsteps,
            real shat, real lx, int ky0_ind,
            IndexVector p = IndexVector(0), FaultCriterion *faultCrit = (new StaticFaults({0,IndexVector(0),IndexVector(0)})),
            IndexType numSpecies = 1, bool GENE_Global = false, bool GENE_Linear = true);

  GeneTask();

  virtual ~GeneTask();
  /**
   * This method does not acutally run GENE but just moves to folder of the GENE task.
   * lcomm is the local communicator of the process group.
   */
  void run( CommunicatorType lcomm );

  /**
   * This method changes the folder to the folder of the task
   * lcomm is the local communicator of the process group.
   */
  void changeDir(CommunicatorType lcomm);

  //void init(CommunicatorType lcomm);

  /**
   * This method initializes the task
   * lcomm is the local communicator of the process group.
   * decomposition is the spatial decomposition of the component grid
   */
  void init(CommunicatorType lcomm, std::vector<IndexVector> decomposition = std::vector<IndexVector>());

  /**
   * This method returns the decomposition of the grid of the specified species
   */
  std::vector<IndexVector> getDecomposition(int species){
      return dfgVector_[species]->getDecomposition();
  }

  /**
   * This method is used to decide if the execution of the task should fail
   */
  void decideToKill();

  /**
   * Returns the path of the task
   */
  inline const std::string& getPath() const;

  inline GeneLocalCheckpoint& getLocalCheckpoint();
  /**
   * This method writes the GENE grid to the local checkpoint
   */
  void writeLocalCheckpoint(  GeneComplex* data,
                              size_t size,
                              std::vector<size_t>& sizes,
                              std::vector<size_t>& bounds );
  /**
   * Initializes Local Checkpoint; This is only necessarry in case of a fault occured and
   * we need to reinitialize the task.
   */
  void InitLocalCheckpoint(size_t size,
      std::vector<size_t>& sizes,
      std::vector<size_t>& bounds );
  /*
   * Gather GENE checkpoint distributed over process group on process
   * with localRootID and convert to FullGrid fg for speciefied species. The actual full grid
   * will only be created on the process with localRootID.
   */
  void getFullGrid( FullGrid<CombiDataType>& fg, RankType lroot,
                    CommunicatorType lcomm, int species);
  /**
   * Returns the distributed full grid of the specified specie
   */
  DistributedFullGrid<CombiDataType>& getDistributedFullGrid(int specie);

  /*
   * Convert fg to GeneGrid and scatter over processes of pgroup. The fullgrid
   * must be available on localRootID.
   */
  void setLocalCheckpoint( FullGrid<complex>& fg,
			   CommunicatorType lcomm, RankType localRootID );

  /*
   * save a fullgrid in GENE's checkpoint format
   */
  static void saveCheckpoint( FullGrid<complex>& fg,
			      const char* filename  );

  /**
   * sets the dfg content to zero
   */
  void setZero();

  /**
   * normal initializiation of DFG at the beginning or in case the processors did not change
   */
  void initDFG( CommunicatorType comm, std::vector<IndexVector>& decomposition );

  /**
   * initializes DFG to a new version and destroys old DFG
   * necessary if set of communicators changes
   */
  void initDFG2( CommunicatorType comm, std::vector<IndexVector>& decomposition );

  /** This method copies and converts data from local checkpoint to dfg
   * This is used for importing data from GENE to our framework
   * */
  void setDFG();
  /**
   * sets boundary parameters for global simulations
   * is directly set in GENE execution before writing checkpoint
   */
  void setBoundaryParameters(double *C_y, int *size_Cy, double *q_prof, int *size_q, double *kymin ){
    C_y_ = C_y;
    size_Cy_ = *size_Cy;
    q_prof_= q_prof;
    size_q_ = *size_q;
    kymin_ = *kymin;
  }
  /* This method copies dfg data to local checkpoint
   * This is used to get the data from our framework to GENE
   * */
  void getDFG();

  /**
   * This method applies some normalization routines to the DFG
   */
  void normalizeDFG(int species);

  /**
   * This method can be used to set the nrg of the Gene simulation.
   */
  inline void setNrg(real nrg);
  /**
   * Sets the total number of timesteps computed so far. Used in case of restart of component grids
   * during fault recovery. Only valid if combitime is not used
   */
  inline void setStepsTotal( size_t stepsTotal );
  /**
   * Sets the current combination step
   */
  inline void setCombiStep(int ncombi);
  /**
   * Return boolean to indicate whether GeneTask is initialized.
   */
  inline bool isInitialized(){
      return initialized_;
  }
  /**
   * Return boolean to indicate whether checkpoint is initialized.
   */
  inline bool checkIsInitialized(){
      return checkpointInitialized_;
  }
  /**
   * This method copies the gyromatrix from GENE to our framework
   * for future reuse
   */
  void write_gyromatrix(GeneComplex* sparse_gyromatrix_buffer,
    int size);
  /**
   * This method can be used to load the gyromatrix from the memory
   * in order to initialize it in GENE
   */
  void load_gyromatrix(GeneComplex* sparse_gyromatrix_buffer,
      int size);
  /**
   * This function returns a boolean valua that indicates whether the gyromatrix
   * is already stored in memory.
   */
  bool is_gyromatrix_buffered(){
    return gyromatrix_buffered_;
  }
  /**
   * This can be used to delete the gyromatrix. This is necessary if size changes.
   * If Petsc is used size of gyromatrix changes during initialization.
   */
  void delete_gyromatrix(){
    gyromatrix_buffered_ = false;
    delete[] gyromatrix_buffer_;
  }
  /**
   * Returns the time that is simulated between combinations.
   * This is only used in case we do not want to use a fixed number of timesteps
   * but a fixed period of time between combinations for each component grids.
   */
  real getCombiTime(){
    return combitime_;
  }
  /**
   * Returns the current time in the simulation. This is uses to update the time in GENE after restart.
   */
  real getCurrentTime(){
    return currentTime_;
  }
  /**
   * Sets the current time in the simulation. This is uses to update the time in GENE after restart.
   */
  void setCurrentTime(real currentTime){
    currentTime_ = currentTime;
  }
  /**
   * Returns the current timestep in the simulation. This is uses to update the timestep in GENE after restart.
   */
  real getCurrentTimestep(){
    return currentTimestep_;
  }
  /**
   * Sets the current timestep in the simulation. This is uses to update the timestep in GENE after restart.
   */
  void setCurrentTimestep(real currentTimestep){
    currentTimestep_ = currentTimestep;
  }
private:
  friend class boost::serialization::access;

  /**
   * This method is used to apply the boundary condition to the imported GENE grid.
   * This is necessary as GENE only simulates the 2^n points exluding the right boundary of the domain
   * This "2^n+1"-st point is added by applying the boundary condition of GENE
   */
  void adaptBoundaryZ(int species);
  /**
   * This method is used if there is no parallelization in z direction
   */
  void adaptBoundaryZlocal(int species);
  /**
   * This method is used if there is parallelization in z direction
   */
  void adaptBoundaryZglobal(int species);
  /**
   * This is the z boundary kernel used in global and local methods.
   * Source specifies the data that is used to fill the target area
   * of the grid according to the boundary condition.
   * The boundary condition specific offset and factor is calculated
   * here.
   */
  void adaptBoundaryZKernel(MultiArrayRef6& sourceData, MultiArrayRef6& targetData, int species);
  /**
   * Calculates offset and multiplication factor from source to target in adapt boundary kernel
   */
  void getOffsetAndFactor( IndexType& xoffset, CombiDataType& factor, IndexType l = 1, IndexType x = 0 );

  // following variables are set in manager and thus need to be included in
  // serialize function
  std::string path_;    // directory in which task should be executed
  real dt_;
  real combitime_; //simulation time interval between combinations
  /**number of time-steps in between two combinations (is set very large in case combitime should be used);
   * this requires equal time-steps for every component grid
   */
  size_t nsteps_;
  size_t stepsTotal_; //number of time-steps done so far (there might be multiple timesteps in between two combinations)
  size_t combiStep_; //number of combinations done so far
  IndexVector p_;
  /**
   * some Gene specific physical parameters
   */
  real shat_;
  real kymin_;
  real lx_;
  real x0_;
  int ky0_ind_;
  double *C_y_;
  int size_Cy_ ;
  double *q_prof_;
  int size_q_;

  GeneComplex* gyromatrix_buffer_; //buffer for gyromatrix
  size_t gyromatrix_buffer_size_; //buffersize of gyromatrix

  bool gyromatrix_buffered_ = false; //indicates if gyromatrix is buffered
  // following variables are only accessed in worker and do not need to be
  // serialized
  GeneLocalCheckpoint checkpoint_;
  std::vector<DistributedFullGrid<CombiDataType> *> dfgVector_;
  real nrg_;

  bool initialized_; //indicates if GeneTask is initialized
  bool checkpointInitialized_; //indicates if checkpoint is initialized

  int nspecies_; //number of species
  MPI_Request * requestArray_;
  std::vector<CombiDataType *> receiveBufferArray_;
  bool _GENE_Global; //indicates whether we are doing global GENE simulations
  bool _GENE_Linear; //indicates whether we are doing linear GENE simulations
  /*
   * simulation time specific parameters
   */
  real currentTime_; //current time in the simulation
  real currentTimestep_; //curent time step length in the simulation

 // std::chrono::high_resolution_clock::time_point  startTimeIteration_;

  // serialize
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & boost::serialization::base_object<Task>(*this);
    ar & path_;
    ar & dt_;
    ar & combitime_;
    ar & nsteps_;
    ar & stepsTotal_;
    ar & combiStep_;
    ar & p_;
    ar & shat_;
    ar & kymin_;
    ar & lx_;
    ar & x0_;
    ar & ky0_ind_;
    ar & nspecies_;
    ar & _GENE_Global;
    ar & _GENE_Linear;
    ar & currentTime_;
  }
};


inline const std::string& GeneTask::getPath() const{
  return path_;
}


inline std::ostream& operator<<( std::ostream& os, const GeneTask &t ){
  os  << "GeneTask:\n"
      << "\t LevelVector = " << t.getLevelVector() << "\n"
      << "\t Path = " << t.getPath();

  return os;
}



inline void GeneTask::setStepsTotal( size_t stepsTotal ) {
    stepsTotal_ = stepsTotal;
}

inline GeneLocalCheckpoint& GeneTask::getLocalCheckpoint(){
  return checkpoint_;
}


inline void GeneTask::setNrg(real nrg){
  nrg_ = nrg;

  MASTER_EXCLUSIVE_SECTION{
    std::cout << "task " << this->getID() << " nrg = " << nrg_ << std::endl;
  }
}
inline void GeneTask::setCombiStep(int ncombi){
  combiStep_ = ncombi;
}



} /* namespace combigrid */




#endif /* GENETASK_HPP_ */
