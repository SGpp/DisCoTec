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
            real shat, real kymin, real lx, int ky0_ind,
            IndexVector p = IndexVector(0), FaultCriterion *faultCrit = (new StaticFaults({0,IndexVector(0),IndexVector(0)})),
            IndexType numSpecies = 1, bool GENE_Global = false, bool GENE_Linear = true);

  GeneTask();

  virtual ~GeneTask();

  void run( CommunicatorType lcomm );

  void changeDir(CommunicatorType lcomm);

  //void init(CommunicatorType lcomm);

  void init(CommunicatorType lcomm, std::vector<IndexVector> decomposition = std::vector<IndexVector>());

  std::vector<IndexVector> getDecomposition(int species){
      return dfgVector_[species]->getDecomposition();
  }

  void decideToKill();

  inline const std::string& getPath() const;

  inline GeneLocalCheckpoint& getLocalCheckpoint();

  void writeLocalCheckpoint(  GeneComplex* data,
                              size_t size,
                              std::vector<size_t>& sizes,
                              std::vector<size_t>& bounds );

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


  void setDFG();
  /**
   * sets boundary parameters for global simulations
   * is directly set in GENE execution before writing checkpoint
   */
  void setBoundaryParameters(double *C_y, int *size_Cy, double *q_prof, int *size_q ){
    C_y_ = C_y;
    size_Cy_ = *size_Cy;
    q_prof_= q_prof;
    size_q_ = *size_q;
  }

  void getDFG();

  void normalizeDFG(int species);

  inline void setNrg(real nrg);

  inline void setStepsTotal( size_t stepsTotal );

  inline void setCombiStep(int ncombi);
  inline bool isInitialized(){
      return initialized_;
  }

  inline bool checkIsInitialized(){
      return checkpointInitialized_;
  }
  void write_gyromatrix(GeneComplex* sparse_gyromatrix_buffer,
    int size);
  void load_gyromatrix(GeneComplex* sparse_gyromatrix_buffer,
      int size);
  bool is_gyromatrix_buffered(){
    return gyromatrix_buffered_;
  }

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
  real getCurrentTime(){
    return currentTime_;
  }
  void setCurrentTime(real currentTime){
    currentTime_ = currentTime;
  }
  real getCurrentTimestep(){
    return currentTimestep_;
  }
  void setCurrentTimestep(real currentTimestep){
    currentTimestep_ = currentTimestep;
  }
private:
  friend class boost::serialization::access;

  void adaptBoundaryZ(int species);

  void adaptBoundaryZlocal(int species);

  void adaptBoundaryZglobal(int species);

  void adaptBoundaryZKernel(MultiArrayRef6& sourceData, MultiArrayRef6& targetData, int species);

  void getOffsetAndFactor( IndexType& xoffset, CombiDataType& factor, IndexType l = 1, IndexType x = 0 );

  inline bool failNow( const int& globalRank );


  // following variables are set in manager and thus need to be included in
  // serialize function
  std::string path_;    // directory in which task should be executed
  real dt_;
  real combitime_; //simulation time interval between combinations
  size_t nsteps_;
  size_t stepsTotal_;
  size_t combiStep_;
  IndexVector p_;

  real shat_;
  real kymin_;
  real lx_;
  real x0_;
  int ky0_ind_;
  GeneComplex* gyromatrix_buffer_; //buffer for gyromatrix
  size_t gyromatrix_buffer_size_; //buffersize of gyromatrix

  bool gyromatrix_buffered_ = false; //indicates if gyromatrix is buffered
  // following variables are only accessed in worker and do not need to be
  // serialized
  GeneLocalCheckpoint checkpoint_;
  std::vector<DistributedFullGrid<CombiDataType> *> dfgVector_;
  real nrg_;

  bool initialized_;
  bool checkpointInitialized_;
  //number of species
  int nspecies_;
  MPI_Request * requestArray_;
  std::vector<CombiDataType *> receiveBufferArray_;
  bool _GENE_Global;
  bool _GENE_Linear;
  double *C_y_;
  int size_Cy_ ;
  double *q_prof_;
  int size_q_;
  real currentTime_;
  real currentTimestep_;

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
    ar & currentTimestep_;

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
