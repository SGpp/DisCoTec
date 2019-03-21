/*
 * ProcessGroupManager.hpp
 *
 *  Created on: Jul 17, 2014
 *      Author: heenemo
 */

#ifndef PROCESSGROUPMANAGER_HPP_
#define PROCESSGROUPMANAGER_HPP_

#include <stdexcept>
#include <string>
#include <vector>

#include "sgpp/distributedcombigrid/combicom/CombiCom.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupSignals.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/mpi_fault_simulator/MPI-FT.h"
#include "sgpp/distributedcombigrid/task/Task.hpp"
#include "sgpp/distributedcombigrid/third_level/ThirdLevelUtils.hpp"

namespace combigrid {
class ProcessGroupManager {
 public:
  ProcessGroupManager() = delete;

  explicit ProcessGroupManager(RankType pgroupRootID);

  ProcessGroupManager(ProcessGroupManager const&) = delete;

  ProcessGroupManager& operator=(ProcessGroupManager const&) = delete;

  bool runfirst(Task* t);

  bool runnext();

  bool exit();

  /* non-blocking call to retrieve status of process group */
  inline StatusType getStatus();

  inline void setStatus(StatusType status);

  /* blocks until process group finished computation */
  inline StatusType waitStatus();

  inline complex eval(const std::vector<real>& coords);

  inline const TaskContainer& getTaskContainer() const;

  inline void removeTask(Task* t);

  /*
  * Third Level reduce: receives the remote data first, combines with local and
  * sends it back to remote.
  */
  template<typename FG_ELEMENT>
  bool
  reduceUniformThirdLevelRecvFirst(const ThirdLevelUtils& thirdLevel,
                                          CombiParameters& params);

  /*
   * Third Level reduce: sends the global reduced data to remote first and
   * recieves and integrates the remotely combined data.
   */
  template<typename FG_ELEMENT>
  bool
  reduceUniformThirdLevelSendFirst(const ThirdLevelUtils& thirdLevel,
                                         CombiParameters& params);

  // distribute common ss  after thirdLevelReduce and integrate them into dfg
  bool
  integrateCommonSS();


  bool
  combineLocalAndGlobal();

  bool combine();

  bool broadcastUpdatedDSG();

  template <typename FG_ELEMENT>
  bool combineFG(FullGrid<FG_ELEMENT>& fg);

  template <typename FG_ELEMENT>
  bool gridEval(FullGrid<FG_ELEMENT>& fg);

  bool gridGather(LevelVector& leval);

  bool updateCombiParameters(CombiParameters& params);

  bool getDSGFromProcessGroup(std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>& outbound);

  bool initiateGetDSGFromProcessGroup();

  bool getAndSetDSGInProcessGroup(std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>& inoutbound);

  bool initiateGetAndSetDSGInProcessGroup();

  size_t getDSGFromNextProcess(std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>& outbound);

  bool addDSGToProcessGroup(std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>& inbound);

  size_t addDSGToNextProcess(std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>& inbound);

  /* Check if group fault occured at this combination step using the fault simulator */
  bool isGroupFault();

  bool addTask(Task*);
  bool refreshTask(Task*);

  // resets tasks only on workers not in group manager
  bool resetTasksWorker();

  bool recompute(Task*);

  bool recoverCommunicators();

  bool parallelEval(const LevelVector& leval, std::string& filename);

 private:
  RankType pgroupRootID_;  // rank in GlobalComm of the master process of this group

  TaskContainer tasks_;

  StatusType status_;

  MPI_Request statusRequest_;

  simft::Sim_FT_MPI_Request statusRequestFT_;

  void recvStatus();

  // Helper functions for Communication with ProcessGroups
  bool storeTaskReferenceAndSendTaskToProcessGroup(Task* t, SignalType signal);

  void storeTaskReference(Task* t);

  bool sendTaskToProcessGroup(Task* t, SignalType signal);

  void sendSignalAndReceive(SignalType signal);

  void sendSignalToProcessGroup(SignalType signal);

  void sendSignalToProcess(SignalType signal, RankType rank);

  inline void setProcessGroupBusyAndReceive();

  /* sets the rank of the process group's master in global comm. should only
   * be called by ProcessManager.
   */
  friend class ProcessManager;
  inline void setMasterRank(int pGroupRootID);
};

typedef std::shared_ptr<ProcessGroupManager> ProcessGroupManagerID;
typedef std::vector<ProcessGroupManagerID> ProcessGroupManagerContainer;

inline StatusType ProcessGroupManager::getStatus() {
  if (status_ == PROCESS_GROUP_WAIT) return PROCESS_GROUP_WAIT;

  if (status_ == PROCESS_GROUP_FAIL) return PROCESS_GROUP_FAIL;
  //  std::cout << "status before test is " << status_ << " \n";

  // if the process group is busy we need
  if (status_ == PROCESS_GROUP_BUSY) {
    if (ENABLE_FT) {
      int flag = -1;
      simft::Sim_FT_MPI_Status stat;
      int err = simft::Sim_FT_MPI_Test(&statusRequestFT_, &flag, &stat);

      if (err == MPI_ERR_PROC_FAILED) status_ = PROCESS_GROUP_FAIL;
    } else {
      /* todo: actually MPI_TEST is not really necessary here. However, i think it
       * might be a good idea to have this here. MPI_Test might run a system
       * call which enables the OS to switch to the MPI system.
       */
      int flag;
      MPI_Test(&statusRequest_, &flag, MPI_STATUS_IGNORE);
    }
  }
  //  std::cout << "status is " << status_ << " \n";
  assert(status_ >= 0);  // check for invalid values
  assert(status_ <= 2);
  return status_;
}
inline void ProcessGroupManager::setStatus(StatusType status) { status_ = status; }

inline StatusType ProcessGroupManager::waitStatus() {
  if (status_ == PROCESS_GROUP_WAIT) return PROCESS_GROUP_WAIT;

  if (status_ == PROCESS_GROUP_FAIL) return PROCESS_GROUP_FAIL;

  if (status_ == PROCESS_GROUP_BUSY) {
    if (ENABLE_FT) {
      simft::Sim_FT_MPI_Status stat;
      int err = simft::Sim_FT_MPI_Wait(&statusRequestFT_, &stat);
      if (err == MPI_ERR_PROC_FAILED) status_ = PROCESS_GROUP_FAIL;
    } else {
      /* todo: actually MPI_TEST is not really necessary here. However, i think it
       * might be a good idea to have this here. MPI_Test might run a system
       * call which enables the OS to switch to the MPI system.
       */
      MPI_Wait(&statusRequest_, MPI_STATUS_IGNORE);
    }
  }
  assert(status_ >= 0);  // check for invalid values
  assert(status_ <= 2);
  return status_;
}

inline complex ProcessGroupManager::eval(const std::vector<real>& x) {
  // todo: implement
  throw std::runtime_error("Not yet implemented: " + std::string(x.begin(), x.end()) + __FILE__ +
                           " " + std::to_string(__LINE__));
  return complex(0.0, 0.0);
}

inline const TaskContainer& ProcessGroupManager::getTaskContainer() const { return tasks_; }

inline void ProcessGroupManager::removeTask(Task* t) {
  std::vector<Task*>::iterator position = std::find(tasks_.begin(), tasks_.end(), t);
  if (position != tasks_.end()) {  // == task_.end() means the element was not found
    tasks_.erase(position);
    std::cout << "Removing task" << t->getID() << " " << theMPISystem()->getWorldRank() << " !\n";

  } else {
    std::cout << "Error could not remove task" << t->getID() << " "
              << theMPISystem()->getWorldRank() << " !\n";
  }
}

template <typename FG_ELEMENT>
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

template <typename FG_ELEMENT>
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

inline bool ProcessGroupManager::gridGather(LevelVector& leval) {
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
bool ProcessGroupManager::reduceUniformThirdLevelSendFirst(const ThirdLevelUtils& thirdLevel,
                                                                  CombiParameters& params) {
  // can only send sync signal when in wait state
  assert(status_ == PROCESS_GROUP_WAIT);

  SignalType signal = COMBINE_UNIFORM_THIRD_LEVEL_SEND_FIRST;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());

  // set status
  status_ = PROCESS_GROUP_BUSY;

  const std::vector<CommunicatorType>& thirdLevelComms = theMPISystem()->getThirdLevelComms();
  assert(theMPISystem()->getNumProcs() == thirdLevelComms.size() &&
      "initialisation of third level communicator failed");

  // receive serialized grid from worker in third level comm (rank 0) and send it to remote
  IndexType numGrids = params.getNumGrids();
  for (size_t p = 0; p < theMPISystem()->getNumProcs(); p++) {
    const CommunicatorType& comm = thirdLevelComms[p];
    std::cout << "forwarding grids from worker with rank: " << p << std::endl;
    for (IndexType g = 0; g < numGrids; g++)
      thirdLevel.sendDSGUniformSerialized(recvDSGUniformSerialized(0, comm));
  }
  thirdLevel.signalReady();

  // receive combined data from remote and send it to worker
  for (size_t p = 0; p < theMPISystem()->getNumProcs(); p++) {
    const CommunicatorType& comm = thirdLevelComms[p];
    for (IndexType g = 0; g < numGrids; g++) {
      std::string serializedDSGU(thirdLevel.recvDSGUniformSerialized());
      MPI_Send(serializedDSGU.data(), static_cast<int>(serializedDSGU.size()),
          MPI_CHAR, 0, SEND_DSG_TO_MANAGER, comm);
    }
  }

  // start non-blocking MPI_IRecv to receive status
  recvStatus();

  return true;
}

/*
 * TODO Till now, after sending reduced data, this proc must wait until remote
 * process sends grid to worker, as we send dsgus in rounds.
 * This could be circumvented by sending and receiving the whole serialized grids as one part.
 * The only benefit would be that one system can terminate the simulation earlier.
 * During combinations the systems have to wait for each other anyways.
 */
template<typename FG_ELEMENT>
bool ProcessGroupManager::reduceUniformThirdLevelRecvFirst(const ThirdLevelUtils& thirdLevel,
                                                                 CombiParameters& params) {
  assert(status_ == PROCESS_GROUP_WAIT);

  SignalType signal = COMBINE_UNIFORM_THIRD_LEVEL_RECV_FIRST;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());

  // set status
  status_ = PROCESS_GROUP_BUSY;

  const std::vector<CommunicatorType>& thirdLevelComms = theMPISystem()->getThirdLevelComms();
  assert(theMPISystem()->getNumProcs() == thirdLevelComms.size() &&
      "initialisation of third level communicator failed");
  // can only send sync signal when in wait state

  // receive dsgus from remote and send them to workers
  size_t numGrids = static_cast<size_t>(params.getNumGrids());
  for (size_t p = 0; p < theMPISystem()->getNumProcs(); p++) {
    const CommunicatorType& comm = thirdLevelComms[p];
    assert(comm != MPI_COMM_NULL);

    for (size_t g = 0; g < numGrids; g++) {
      std::string serializedDSGU = thirdLevel.recvDSGUniformSerialized();
      MPI_Send(serializedDSGU.data(), static_cast<size_t>(serializedDSGU.size()),
          MPI_CHAR, 0, SEND_DSG_TO_MANAGER, comm);
    }
  }

  // receive reduced dsgus from workers and send them back to remote system
  for (size_t p = 0; p < theMPISystem()->getNumProcs(); p++) {
    const CommunicatorType& comm = theMPISystem()->getThirdLevelComms()[p];
    assert(comm != MPI_COMM_NULL);
    for (size_t g = 0; g < numGrids; g++)
      thirdLevel.sendDSGUniformSerialized(recvDSGUniformSerialized(0, comm));
  }

  thirdLevel.signalReady();

  // start non-blocking MPI_IRecv to receive status
  recvStatus();

  return true;
}

inline void ProcessGroupManager::setMasterRank(int pGroupRootID) { pgroupRootID_ = pGroupRootID; }

} /* namespace combigrid */

#endif /* PROCESSGROUPMANAGER_HPP_ */
