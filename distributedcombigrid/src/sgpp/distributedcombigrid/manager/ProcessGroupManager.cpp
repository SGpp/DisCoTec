/*
 * ProcessGroupManager.cpp
 *
 *  Created on: Jul 17, 2014
 *      Author: heenemo
*/

#include "sgpp/distributedcombigrid/manager/ProcessGroupManager.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/mpi/MPIUtils.hpp"
#include "sgpp/distributedcombigrid/mpi_fault_simulator/MPI-FT.h"

namespace combigrid {
ProcessGroupManager::ProcessGroupManager(RankType pgroupRootID)
    : pgroupRootID_(pgroupRootID),
      status_(PROCESS_GROUP_WAIT),
      statusRequest_(MPI_Request()),
      statusRequestFT_(nullptr) {}

bool ProcessGroupManager::runfirst(Task* t) {
  return storeTaskReferenceAndSendTaskToProcessGroup(t, RUN_FIRST);
}

bool ProcessGroupManager::storeTaskReferenceAndSendTaskToProcessGroup(Task* t, SignalType signal) {
  // first check status
  // tying to add a task to a busy group is an invalid operation
  // and should be avoided
  if (status_ != PROCESS_GROUP_WAIT) return false;

  storeTaskReference(t);
  return sendTaskToProcessGroup(t, signal);
}

void ProcessGroupManager::storeTaskReference(Task* t) {
  // add task to list of tasks managed by this pgroup
  tasks_.push_back(t);
}

bool ProcessGroupManager::sendTaskToProcessGroup(Task* t, SignalType signal) {
  // send signal to pgroup
  sendSignalToProcessGroup(signal);

  // send task
  Task::send(&t, pgroupRootID_, theMPISystem()->getGlobalComm());

  setProcessGroupBusyAndReceive();

  // only return true if task successfully sent to pgroup
  return true;
}

void ProcessGroupManager::sendSignalAndReceive(SignalType signal) {
  sendSignalToProcessGroup(signal);
  setProcessGroupBusyAndReceive();
}

void ProcessGroupManager::sendSignalToProcessGroup(SignalType signal) {
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());
}

void ProcessGroupManager::sendSignalToProcess(SignalType signal, RankType rank) { //TODO send only to process in this pgroup
  MPI_Send(&signal, 1, MPI_INT, rank, signalTag, theMPISystem()->getGlobalComm());
}

inline void ProcessGroupManager::setProcessGroupBusyAndReceive() {
  // set status
  status_ = PROCESS_GROUP_BUSY;

  // start non-blocking MPI_IRecv to receive status
  recvStatus();
}

bool ProcessGroupManager::runnext() {
  // first check status
  // trying to send a command to a busy group is an invalid operation
  // and should be avoided
  assert(status_ == PROCESS_GROUP_WAIT);

  if (tasks_.size() == 0) return false;

  sendSignalAndReceive(RUN_NEXT);

  return true;
}

bool ProcessGroupManager::exit() {
  // can only send exit signal when in wait state
  if (status_ != PROCESS_GROUP_WAIT) return false;

  sendSignalAndReceive(EXIT);
  return true;
}

bool ProcessGroupManager::combine() {
  // can only send sync signal when in wait state
  assert(status_ == PROCESS_GROUP_WAIT);

  sendSignalAndReceive(COMBINE);

  return true;
}

bool ProcessGroupManager::sendDSGUniformToRemote(const ThirdLevelUtils& thirdLevel,
                                                       CombiParameters& params) {
  sendSignalSendDSGU();
  forwardDSGUFromPGToRemote(thirdLevel, params);

  // start non-blocking MPI_IRecv to receive status
  recvStatus();
  return true;
}

bool ProcessGroupManager::recvDSGUniformFromRemote(const ThirdLevelUtils& thirdLevel,
                                                         CombiParameters& params) {
  sendSignalRecvDSGU();
  forwardDSGUFromRemoteToPG(thirdLevel, params);

  // start non-blocking MPI_IRecv to receive status
  recvStatus();
  return true;
}

bool ProcessGroupManager::recvAndAddDSGUniformFromRemote(const ThirdLevelUtils& thirdLevel, 
                                                               CombiParameters& params) {
  sendSignalRecvAndAdd();
  forwardDSGUFromRemoteToPG(thirdLevel, params);

  // start non-blocking MPI_IRecv to receive status
  recvStatus();
  return true;
}


bool ProcessGroupManager::sendSignalSendDSGU() {
  assert(status_ == PROCESS_GROUP_WAIT);

  SignalType signal = SEND_DSGU_TO_MANAGER;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());

  // set status
  status_ = PROCESS_GROUP_BUSY;
  return true;
}

bool ProcessGroupManager::sendSignalRecvDSGU() {
  assert(status_ == PROCESS_GROUP_WAIT);

  SignalType signal = RECV_DSGU_FROM_MANAGER;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());

  // set status
  status_ = PROCESS_GROUP_BUSY;
  return true;
}

bool ProcessGroupManager::sendSignalRecvAndAdd() {
  assert(status_ == PROCESS_GROUP_WAIT);

  SignalType signal = RECV_AND_ADD_DSGU_FROM_MANAGER;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());

  // set status
  status_ = PROCESS_GROUP_BUSY;
  return true;
}

bool ProcessGroupManager::forwardDSGUFromRemoteToPG(const ThirdLevelUtils& thirdLevel,
                                                          CombiParameters& params) {
  const std::vector<CommunicatorType>& thirdLevelComms = theMPISystem()->getThirdLevelComms();
  assert(theMPISystem()->getNumGroups() == thirdLevelComms.size() &&
      "initialisation of third level communicator failed");

  // receive dsgus from remote and forward them to workers
  size_t numGrids = static_cast<size_t>(params.getNumGrids());
  for (size_t p = 0; p < theMPISystem()->getNumProcs(); p++) {
    const CommunicatorType& comm = thirdLevelComms[p];
    assert(comm != MPI_COMM_NULL);

    for (size_t g = 0; g < numGrids; g++) {
      std::string serializedDSGU = thirdLevel.recvDSGUniformSerialized();
      MPI_Send(serializedDSGU.data(), static_cast<int>(serializedDSGU.size()),
          MPI_CHAR, 0, sendDSGTag, comm);
    }
  }
  thirdLevel.signalReady();
  return true;
}

bool ProcessGroupManager::forwardDSGUFromPGToRemote(const ThirdLevelUtils& thirdLevel,
                                                          CombiParameters& params) {
  const std::vector<CommunicatorType>& thirdLevelComms = theMPISystem()->getThirdLevelComms();
  assert(theMPISystem()->getNumGroups() == thirdLevelComms.size() &&
      "initialisation of third level communicator failed");

  // receive grid from worker and send it to remote
  IndexType numGrids = params.getNumGrids();
  for (size_t p = 0; p < theMPISystem()->getNumProcs(); p++) {
    const CommunicatorType& comm = thirdLevelComms[p];
    std::cout << "forwarding grids from worker with rank: " << p << std::endl;
    for (IndexType g = 0; g < numGrids; g++)
      thirdLevel.sendDSGUniformSerialized(recvDSGUniformSerialized(0, comm));
  }
  thirdLevel.signalReady();
  return true;
}


bool ProcessGroupManager::combineLocalAndGlobal() {
  // can only send sync signal when in wait state
  assert(status_ == PROCESS_GROUP_WAIT);

  SignalType signal = COMBINE_LOCAL_AND_GLOBAL;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());

  // set status
  status_ = PROCESS_GROUP_BUSY;

  // start non-blocking MPI_IRecv to receive status
  recvStatus();

  return true;
}

bool ProcessGroupManager::integrateCombinedDSGUniform() {
  // can only send sync signal when in wait state
  assert(status_ == PROCESS_GROUP_WAIT);

  SignalType signal = INTEGRATE_COMBINED_DSGU;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());

  // set status
  status_ = PROCESS_GROUP_BUSY;

  // start non-blocking MPI_IRecv to receive status
  recvStatus();

  return true;
}


bool ProcessGroupManager::updateCombiParameters(CombiParameters& params) {
  // can only send sync signal when in wait state
  assert(status_ == PROCESS_GROUP_WAIT);

  sendSignalToProcessGroup(UPDATE_COMBI_PARAMETERS);

  // send combiparameters
  // std::cout << "sending class \n";
  MPIUtils::sendClass(&params, pgroupRootID_, theMPISystem()->getGlobalComm());

  setProcessGroupBusyAndReceive();
  // std::cout << "manager received status \n";
  return true;
}

void getDSGFrom(std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>& outbound, size_t from) {
  size_t numGrids = outbound.size();
  for (size_t g = 0; g < numGrids; ++g) {
    outbound[g].reset(recvDSGUniform<CombiDataType>(from, theMPISystem()->getWorldComm()));
  }
}

bool ProcessGroupManager::addTask(Task* t) {
  return storeTaskReferenceAndSendTaskToProcessGroup(t, ADD_TASK);
}

bool ProcessGroupManager::refreshTask(Task* t) {
  // first check status
  // tying to add a task to a busy group is an invalid operation
  // and should be avoided
  if (status_ != PROCESS_GROUP_WAIT) {
    std::cout << "refreshing failed! \n";
    return false;
  }

  sendSignalToProcessGroup(ADD_TASK);

  // send task
  Task::send(&t, pgroupRootID_, theMPISystem()->getGlobalComm());

  setProcessGroupBusyAndReceive();

  // only return true if task successfully send to pgroup
  return true;
}

bool ProcessGroupManager::resetTasksWorker() {
  // first check status
  // tying to reset tasks of a busy group is an invalid operation
  // and should be avoided
  if (status_ != PROCESS_GROUP_WAIT) {
    assert(false);
    // return false;
  }

  // add task to list of tasks managed by this pgroup
  // tasks_.clear(); we do not clear group manager tasks

  sendSignalAndReceive(RESET_TASKS);

  return true;
}

bool ProcessGroupManager::recompute(Task* t) {
  storeTaskReferenceAndSendTaskToProcessGroup(t, RECOMPUTE);

  // only return true if task successfully send to pgroup
  return true;
}

bool ProcessGroupManager::parallelEval(const LevelVector& leval, std::string& filename) {
  // can only send sync signal when in wait state, so check first
  assert(status_ == PROCESS_GROUP_WAIT);

  sendSignalToProcessGroup(PARALLEL_EVAL);

  // send levelvector
  std::vector<int> tmp(leval.begin(), leval.end());
  MPI_Send(&tmp[0], static_cast<int>(tmp.size()), MPI_INT, pgroupRootID_, 0,
           theMPISystem()->getGlobalComm());

  // send filename
  MPIUtils::sendClass(&filename, pgroupRootID_, theMPISystem()->getGlobalComm());

  setProcessGroupBusyAndReceive();

  return true;
}

void ProcessGroupManager::recvStatus() {
  // start non-blocking call to receive status
  if (ENABLE_FT) {
    simft::Sim_FT_MPI_Irecv(&status_, 1, MPI_INT, pgroupRootID_, statusTag,
                            theMPISystem()->getGlobalCommFT(), &statusRequestFT_);
  } else {
    MPI_Irecv(&status_, 1, MPI_INT, pgroupRootID_, statusTag, theMPISystem()->getGlobalComm(),
              &statusRequest_);
  }
}

bool ProcessGroupManager::recoverCommunicators() {
  assert(status_ == PROCESS_GROUP_WAIT);

  sendSignalToProcessGroup(RECOVER_COMM);

  return true;
}

} /* namespace combigrid */
