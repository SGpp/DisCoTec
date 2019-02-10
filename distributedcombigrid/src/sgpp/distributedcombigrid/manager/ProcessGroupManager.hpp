/*
 * ProcessGroupManager.hpp
 *
 *  Created on: Jul 17, 2014
 *      Author: heenemo
 */

#ifndef PROCESSGROUPMANAGER_HPP_
#define PROCESSGROUPMANAGER_HPP_

#include <vector>

#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/combicom/CombiCom.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupSignals.hpp"
#include "sgpp/distributedcombigrid/task/Task.hpp"
#include "sgpp/distributedcombigrid/third_level/ThirdLevelUtils.hpp"

namespace combigrid {

class ProcessGroupManager {

 public:
  ProcessGroupManager() = delete;

  explicit ProcessGroupManager( RankType pgroupRootID );

  ProcessGroupManager( ProcessGroupManager const & ) = delete;

  ProcessGroupManager& operator=( ProcessGroupManager const & ) = delete;

  bool
  runfirst(Task* t);

  bool
  runnext();

  bool
  exit();

  /* non-blocking call to retrieve status of process group */
  inline StatusType
  getStatus();

  /* blocks until process group finished computation */
  inline StatusType
  waitStatus();

  inline complex
  eval(const std::vector<real>& coords);

  inline const TaskContainer&
  getTaskContainer() const;

  bool
  combine();

  // receive shared ss data from remote, combine and send the updated back
  template<typename FG_ELEMENT>
  bool
  exchangeCommonSubspacesThirdLevel(const ThirdLevelUtils& thirdLevel, CombiParameters& params);

  // exchange shared ss data of workers with remote
  template<typename FG_ELEMENT>
  bool
  combineUniformThirdLevel(const ThirdLevelUtils& thirdLevel, CombiParameters& params);

  // distribute common ss  after thirdLevelReduce and integrate them into dfg
  bool
  integrateCommonSS();

  // retrieve shared ss data from workers and forward to remote thirdLevel manager
  // TODO: future usecase like collect on remote
  template<typename FG_ELEMENT>
  bool
  sendallCommonSS(const ThirdLevelUtils& thirdLevel_);

  // retrieve shared ss data from thirdlevel manager and overwrite those at workers.
  // TODO: future usecase like controlling
  template<typename FG_ELEMENT>
  bool
  recvallCommonSS();

  bool
  combineLocalAndGlobal();

  template<typename FG_ELEMENT>
  bool
  combineFG(FullGrid<FG_ELEMENT>& fg);

  template<typename FG_ELEMENT>
  bool
  gridEval(FullGrid<FG_ELEMENT>& fg);

  bool
  gridGather(LevelVector& leval);

  bool
  updateCombiParameters(CombiParameters& params);

  /* Check if group fault occured at this combination step using the fault simulator */
  bool
  isGroupFault();

  bool addTask( Task* );

  bool parallelEval( const LevelVector& leval, std::string& filename );

 private:
  RankType pgroupRootID_; // rank in GlobalComm of the master process of this group

  TaskContainer tasks_;

  StatusType status_;

  MPI_Request statusRequest_;

  std::vector<CombiDataType> allBetas_;

  void recvStatus();

  /* sets the rank of the process group's master in global comm. should only
   * be called by ProcessManager.
   */
  friend class ProcessManager;
  inline void setMasterRank( int pGroupRootID );
};

typedef std::shared_ptr< ProcessGroupManager > ProcessGroupManagerID;
typedef std::vector< ProcessGroupManagerID >
          ProcessGroupManagerContainer;

inline StatusType ProcessGroupManager::getStatus() {
  if( status_ == PROCESS_GROUP_WAIT )
    return PROCESS_GROUP_WAIT;

  if( status_ == PROCESS_GROUP_FAIL )
    return PROCESS_GROUP_FAIL;

  // if the process group is busy we need
  if( status_ == PROCESS_GROUP_BUSY){
    /* todo: actually MPI_TEST is not really necessary here. However, i think it
     * might be a good idea to have this here. MPI_Test might run a system
     * call which enables the OS to switch to the MPI system.
     */
    int flag;
    MPI_Test(&statusRequest_, &flag, MPI_STATUS_IGNORE);
  }

  return status_;
}


inline StatusType ProcessGroupManager::waitStatus() {
  if( status_ == PROCESS_GROUP_WAIT )
    return PROCESS_GROUP_WAIT;

  if( status_ == PROCESS_GROUP_FAIL )
    return PROCESS_GROUP_FAIL;

  if( status_ == PROCESS_GROUP_BUSY){
    /* todo: actually MPI_TEST is not really necessary here. However, i think it
     * might be a good idea to have this here. MPI_Test might run a system
     * call which enables the OS to switch to the MPI system.
     */
    MPI_Wait( &statusRequest_, MPI_STATUS_IGNORE);
  }

  return status_;
}


inline complex ProcessGroupManager::eval(const std::vector<real>& x) {
  //todo: implement
  return complex(0.0, 0.0);
}

inline const TaskContainer&
ProcessGroupManager::getTaskContainer() const {
  return tasks_;
}

template<typename FG_ELEMENT>
bool ProcessGroupManager::gridEval(FullGrid<FG_ELEMENT>& fg) {
  // can only send sync signal when in wait state, so check first
  assert(status_ == PROCESS_GROUP_WAIT);

  // send signal
  SignalType signal = GRID_EVAL;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());

  // send levelvector
  std::vector<int> tmp(fg.getLevels().begin(), fg.getLevels().end());
  MPI_Send(&tmp[0], static_cast<int>(tmp.size()), MPI_INT, pgroupRootID_, 0,
           theMPISystem()->getGlobalComm());

  return true;
}

template<typename FG_ELEMENT>
bool ProcessGroupManager::combineFG(FullGrid<FG_ELEMENT>& fg) {
  // can only send sync signal when in wait state
  assert(status_ == PROCESS_GROUP_WAIT);

  SignalType signal = COMBINE_FG;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());

  // send levelvector
  std::vector<int>& tmp = fg.getLevels();
  MPI_Send(&tmp[0], static_cast<int>(tmp.size()), MPI_INT, pgroupRootID_, 0,
           theMPISystem()->getGlobalComm());

  return true;
}

inline
bool ProcessGroupManager::gridGather(LevelVector& leval) {
  // can only send sync signal when in wait state, so check first
  assert(status_ == PROCESS_GROUP_WAIT);

  // send signal
  SignalType signal = GRID_GATHER;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());

  // send levelvector
  std::vector<int> tmp(leval.begin(), leval.end());
  MPI_Send(&tmp[0], static_cast<int>(tmp.size()), MPI_INT, pgroupRootID_, 0,
           theMPISystem()->getGlobalComm());

  return true;
}

template <typename FG_ELEMENT>
bool ProcessGroupManager::exchangeCommonSubspacesThirdLevel(const ThirdLevelUtils& thirdLevel,
                                                                  CombiParameters& params) {
  // can only send sync signal when in wait state
  assert(status_ == PROCESS_GROUP_WAIT);

  SignalType signal = EXCHANGE_COMMON_SS;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());

  // set status
  status_ = PROCESS_GROUP_BUSY;

  const std::vector<LevelVector> commonSS = params.getCommonSubspaces();
  const size_t numCommonSS = commonSS.size();
  const std::vector<CommunicatorType>& thirdLevelComms = theMPISystem()->getThirdLevelComms();

  MPI_Datatype dtype = abstraction::getMPIDatatype(
                         abstraction::getabstractionDataType<FG_ELEMENT>());

  // receive sizes
  std::vector<std::vector<int>> commonSSPartSizes;
  for (size_t p = 0; p < theMPISystem()->getNumProcs(); p++) {
    std::vector<int> sizes(numCommonSS);
    MPI_Recv(sizes.data(), static_cast<int>(numCommonSS), MPI_INT, 0, 0, thirdLevelComms[p],
        MPI_STATUS_IGNORE);
    commonSSPartSizes.push_back(std::move(sizes));
  }

  // calculate sizes of each workers part
  std::vector<size_t> partSizes(theMPISystem()->getNumProcs());
  size_t maxPartSize = 0;
  for (size_t p = 0; p < theMPISystem()->getNumProcs(); p++) {
    for (size_t ss = 0; ss < numCommonSS; ss++)
      partSizes[p] += static_cast<size_t>(commonSSPartSizes[p][ss]);
    if (maxPartSize < partSizes[p])
      maxPartSize = partSizes[p];
  }

  // reserve buffer with size of greatest part
  std::vector<FG_ELEMENT> commonSSPart;
  commonSSPart.reserve(maxPartSize);

  for (size_t p = 0; p < theMPISystem()->getNumProcs(); p++) {
    const CommunicatorType& comm = thirdLevelComms[p];

    // receive subspace data from worker
    commonSSPart.resize(partSizes[p]);
    size_t stride = 0;
    for (size_t ss = 0; ss < numCommonSS; ss++) {
      MPI_Recv(commonSSPart.data()+stride, commonSSPartSizes[p][ss], dtype, 0, 0,
          comm, MPI_STATUS_IGNORE);
      stride += static_cast<size_t>(commonSSPartSizes[p][ss]);
    }

    // send subspace data to remote
    thirdLevel.sendCommonSSPart(commonSSPart);
    // receive combined data from remote
    thirdLevel.receiveCommonSSPart(commonSSPart);

    // send combined subspace data to worker
    stride = 0;
    for (size_t ss = 0; ss < numCommonSS; ss++) {
      MPI_Send(commonSSPart.data()+stride, commonSSPartSizes[p][ss], dtype, 0, 0, comm);
      stride += static_cast<size_t>(commonSSPartSizes[p][ss]);
    }
  }

  // start non-blocking MPI_IRecv to receive status
  recvStatus();

  return true;
}

// TODO: Check Sizes
template<typename FG_ELEMENT>
bool ProcessGroupManager::combineUniformThirdLevel(const ThirdLevelUtils& thirdLevel,
                                                          CombiParameters& params) {
  // can only send sync signal when in wait state
  assert(status_ == PROCESS_GROUP_WAIT);

  SignalType signal = COMBINE_UNIFORM_THIRD_LEVEL;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());

  // set status
  status_ = PROCESS_GROUP_BUSY;

  const std::vector<LevelVector> commonSS = params.getCommonSubspaces();
  const size_t numCommonSS = commonSS.size();
  const std::vector<CommunicatorType>& thirdLevelComms = theMPISystem()->getThirdLevelComms();

  MPI_Datatype dtype = abstraction::getMPIDatatype(
                         abstraction::getabstractionDataType<FG_ELEMENT>());

  // receive sizes
  std::vector<std::vector<int>> commonSSPartSizes;
  for (size_t p = 0; p < theMPISystem()->getNumProcs(); p++) {
    std::vector<int> sizes(numCommonSS);
    MPI_Recv(sizes.data(), static_cast<int>(numCommonSS), MPI_INT, 0, 0, thirdLevelComms[p],
        MPI_STATUS_IGNORE);
    commonSSPartSizes.push_back(std::move(sizes));
  }

  // calculate sizes of each workers part
  std::vector<size_t> partSizes(theMPISystem()->getNumProcs());
  size_t maxPartSize = 0;
  for (size_t p = 0; p < theMPISystem()->getNumProcs(); p++) {
    for (size_t ss = 0; ss < numCommonSS; ss++)
      partSizes[p] += static_cast<size_t>(commonSSPartSizes[p][ss]);
    if (maxPartSize < partSizes[p])
      maxPartSize = partSizes[p];
  }
  // reserve buffer with size of greatest part
  std::vector<FG_ELEMENT> commonSSPart;
  commonSSPart.reserve(maxPartSize);

  for (size_t p = 0; p < theMPISystem()->getNumProcs(); p++) {
    const CommunicatorType& comm = thirdLevelComms[p];
    assert(comm != MPI_COMM_NULL);

    commonSSPart.resize(partSizes[p]);

    // receive data from remote
    thirdLevel.receiveCommonSSPart(commonSSPart);

    // combine subspaces sequentially
    size_t stride = 0;
    for (size_t ss = 0; ss < numCommonSS; ss++) {
      MPI_Allreduce(MPI_IN_PLACE, commonSSPart.data()+stride, commonSSPartSizes[p][ss], dtype, MPI_SUM, comm);
      stride += static_cast<size_t>(commonSSPartSizes[p][ss]);
    }

    // send data back to remote
    thirdLevel.sendCommonSSPart(commonSSPart);
  }

  // start non-blocking MPI_IRecv to receive status
  recvStatus();

  return true;
}

inline void ProcessGroupManager::setMasterRank( int pGroupRootID ){
  pgroupRootID_ = pGroupRootID;
}

} /* namespace combigrid */

#endif /* PROCESSGROUPMANAGER_HPP_ */
