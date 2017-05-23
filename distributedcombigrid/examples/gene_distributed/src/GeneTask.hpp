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
#include "sgpp/distributedcombigrid/fault_tolerance/FaultCriterion.hpp"
#include "sgpp/distributedcombigrid/fault_tolerance/StaticFaults.hpp"

namespace combigrid {

class GeneTask: public combigrid::Task {
public:
  GeneTask( DimType dim, LevelVector& l, std::vector<bool>& boundary, real coeff,
            LoadModel* loadModel, std::string& path, real dt, size_t nsteps,
            real shat, real kymin, real lx, int ky0_ind,
            IndexVector p = IndexVector(0), FaultCriterion *faultCrit = (new StaticFaults({0,IndexVector(0),IndexVector(0)})));

  GeneTask();

  virtual ~GeneTask();

  void run( CommunicatorType lcomm );

  void changeDir();

  //void init(CommunicatorType lcomm);

  void init(CommunicatorType lcomm, std::vector<IndexVector> decomposition = std::vector<IndexVector>());

  std::vector<IndexVector> getDecomposition(){
      return dfg_->getDecomposition();
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
   * with localRootID and convert to FullGrid fg. The actual full grid
   * will only be created on the process with localRootID.
   */
  void getFullGrid( FullGrid<CombiDataType>& fg, RankType lroot,
                    CommunicatorType lcomm);

  DistributedFullGrid<CombiDataType>& getDistributedFullGrid();

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

  void getDFG();

  void normalizeDFG();

  inline void setNrg(real nrg);

  inline void setStepsTotal( size_t stepsTotal );

  inline void setCombiStep(int ncombi);
  inline bool isInitialized(){
      return initialized_;
  }

  inline bool checkIsInitialized(){
      return checkpointInitialized_;
  }

private:
  friend class boost::serialization::access;

  void adaptBoundaryZ();

  void adaptBoundaryZlocal();

  void adaptBoundaryZglobal();

  void getOffsetAndFactor( IndexType& xoffset, CombiDataType& factor );

  inline bool failNow( const int& globalRank );





  // following variables are set in manager and thus need to be included in
  // serialize function
  std::string path_;    // directory in which task should be executed
  real dt_;
  size_t nsteps_;
  size_t stepsTotal_;
  size_t combiStep_;
  IndexVector p_;

  real shat_;
  real kymin_;
  real lx_;
  int ky0_ind_;
  //fault tolerance info
  FaultCriterion *faultCriterion_;

  // following variables are only accessed in worker and do not need to be
  // serialized
  GeneLocalCheckpoint checkpoint_;
  DistributedFullGrid<CombiDataType>* dfg_;
  real nrg_;

  bool initialized_;
  bool checkpointInitialized_;

  std::chrono::high_resolution_clock::time_point  startTimeIteration_;

  // serialize
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & boost::serialization::base_object<Task>(*this);
    ar & path_;
    ar & dt_;
    ar & nsteps_;
    ar & stepsTotal_;
    ar & combiStep_;
    ar & p_;
    ar & shat_;
    ar & kymin_;
    ar & lx_;
    ar & ky0_ind_;
    ar & faultCriterion_;
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
