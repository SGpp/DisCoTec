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

bool ProcessGroupManager::getDSGFromProcessGroup(std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>& outbound) {
  // this is not really useful, as it immediately overwrites the same outbound grid, so just here for illustration on the idea
  initiateGetDSGFromProcessGroup();
  for (size_t i = 0; i < theMPISystem()->getNumProcs(); ++i) {
    getDSGFromNextProcess(outbound);
  }
  return true;
}

bool ProcessGroupManager::initiateGetDSGFromProcessGroup() {
  // can only send sync signal when in wait state
  assert(status_ == PROCESS_GROUP_WAIT);
  status_ = PROCESS_GROUP_BUSY;

  sendSignalToProcessGroup(SEND_DSG_AND_CONTINUE);

  return true;
}

bool ProcessGroupManager::getAndSetDSGInProcessGroup(std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>& outbound) {
  // this is not really useful, as it immediately overwrites the same outbound grid, so just here for illustration on the idea
  initiateGetAndSetDSGInProcessGroup();
  for (size_t i = 0; i < theMPISystem()->getNumProcs(); ++i) {
    getDSGFromNextProcess(outbound);
  }
  addDSGToProcessGroup(outbound);
  return true;
}

bool ProcessGroupManager::initiateGetAndSetDSGInProcessGroup() {
  // can only send sync signal when in wait state
  assert(status_ == PROCESS_GROUP_WAIT);
  status_ = PROCESS_GROUP_BUSY;

  sendSignalToProcessGroup(SEND_DSG_AND_WAIT_FOR_ADD);

  return true;
}


void getDSGFrom(std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>& outbound, size_t from) {
  size_t numGrids = outbound.size();
  for (size_t g = 0; g < numGrids; ++g) {
    outbound[g].reset(recvDSGUniform<CombiDataType>(from, theMPISystem()->getWorldComm()));
  }
}

size_t ProcessGroupManager::getDSGFromNextProcess(std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>& outbound) {
  // static variable to iterate the processes in the first group
  static size_t from = 0;
  getDSGFrom(outbound, from);
  size_t processesLeft = (theMPISystem()->getNumProcs() - 1) - from;

  from = (from+1) % theMPISystem()->getNumProcs(); //TODO get actual indices in case of not first process group
  if(processesLeft == 0){
    status_ = PROCESS_GROUP_WAIT; //TODO
  }
  
  return processesLeft;
}

bool ProcessGroupManager::addDSGToProcessGroup(std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>& inbound){
  // this is not really useful, as it overwrites with the same inbound grid, so just here for illustration on the idea
  for (size_t i = 0; i < theMPISystem()->getNumProcs(); ++i) {
    addDSGToNextProcess(inbound);
  }
  return true;
}

void addDSGTo(std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>& inbound, size_t to) {
  size_t numGrids = inbound.size();
  for (size_t g = 0; g < numGrids; ++g) {
    sendDSGUniform<CombiDataType>(inbound[g].get(), to, theMPISystem()->getWorldComm());
  }
}

size_t ProcessGroupManager::addDSGToNextProcess(std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>& inbound) {
  // static variable to iterate the processes in the first group
  static size_t to = 0;
  // std::cerr << "adding inbound to " << to << std::endl;

  addDSGTo(inbound, to);

  size_t processesLeft = (theMPISystem()->getNumProcs() - 1) - to;
  if(processesLeft == 0){
    status_ = PROCESS_GROUP_WAIT; //TODO
  }

  to = (to+1) % theMPISystem()->getNumProcs(); //TODO get actual indices in case of not first process group
  return processesLeft;
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
