/*
 * ProcessGroupManager.cpp
 *
 *  Created on: Jul 17, 2014
 *      Author: heenemo
 */

#include "sgpp/distributedcombigrid/manager/ProcessGroupManager.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/mpi/MPIUtils.hpp"

namespace combigrid {
ProcessGroupManager::ProcessGroupManager( RankType pgroupRootID ) :
        pgroupRootID_(pgroupRootID),
        status_(PROCESS_GROUP_WAIT),
        statusRequest_(MPI_Request())
{
}

bool ProcessGroupManager::runfirst(Task* t) {
  // first check status
  // tying to add a task to a busy group is an invalid operation
  // and should be avoided
  if (status_ != PROCESS_GROUP_WAIT)
    return false;

  // add task to list of tasks managed by this pgroup
  tasks_.push_back(t);

  // send runfirst_signal to pgroup
  SignalType signal = RUN_FIRST;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());

  // send task
  Task::send( &t, pgroupRootID_, theMPISystem()->getGlobalComm() );

  // set status
  status_ = PROCESS_GROUP_BUSY;

  // start non-blocking MPI_IRecv to receive status
  recvStatus();

  // only return true if task successfully send to pgroup
  return true;
}

bool ProcessGroupManager::runnext() {
  // first check status
  // trying to send a command to a busy group is an invalid operation
  // and should be avoided
  assert(status_ == PROCESS_GROUP_WAIT);

  if (tasks_.size() == 0)
    return false;

  // send runnext signal
  SignalType signal = RUN_NEXT;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());

  status_ = PROCESS_GROUP_BUSY;

  // start non-blocking MPI_IRecv to receive status
  recvStatus();

  return true;
}

bool ProcessGroupManager::exit() {
  // can only send exit signal when in wait state
  if (status_ != PROCESS_GROUP_WAIT)
    return false;

  // send exit signal
  SignalType signal = EXIT;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());

  status_ = PROCESS_GROUP_BUSY;

  // start non-blocking MPI_IRecv to receive status
  recvStatus();

  return true;
}

bool ProcessGroupManager::combine() {
  // can only send sync signal when in wait state
  assert(status_ == PROCESS_GROUP_WAIT);

  SignalType signal = COMBINE;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());

  // set status
  status_ = PROCESS_GROUP_BUSY;

  // start non-blocking MPI_IRecv to receive status
  recvStatus();

  return true;
}

bool ProcessGroupManager::updateCombiParameters(CombiParameters& params) {
  // can only send sync signal when not in busy state
  assert(status_ != PROCESS_GROUP_BUSY);

  SignalType signal = UPDATE_COMBI_PARAMETERS;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());

  // send combiparameters
  MPIUtils::sendClass(&params, pgroupRootID_, theMPISystem()->getGlobalComm());
  // set status
  status_ = PROCESS_GROUP_BUSY;

  // start non-blocking MPI_IRecv to receive status
  recvStatus();
  return true;
}


bool ProcessGroupManager::addTask( Task* t ) {
  // first check status
  // tying to add a task to a busy group is an invalid operation
  // and should be avoided
  if (status_ != PROCESS_GROUP_WAIT)
    return false;

  // add task to list of tasks managed by this pgroup
  tasks_.push_back(t);

  // send add task signal to pgroup
  SignalType signal = ADD_TASK;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());

  // send task
  Task::send(&t, pgroupRootID_, theMPISystem()->getGlobalComm());

  // set status
  status_ = PROCESS_GROUP_BUSY;

  // start non-blocking MPI_IRecv to receive status
  recvStatus();

  // only return true if task successfully send to pgroup
  return true;
}


bool ProcessGroupManager::recompute( Task* t ) {
  // first check status
  // tying to add a task to a busy group is an invalid operation
  // and should be avoided
  if (status_ != PROCESS_GROUP_WAIT)
    return false;

  // add task to list of tasks managed by this pgroup
  tasks_.push_back(t);

  // send add task signal to pgroup
  SignalType signal = RECOMPUTE;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());

  // send task
  Task::send( &t, pgroupRootID_, theMPISystem()->getGlobalComm() );

  // set status
  status_ = PROCESS_GROUP_BUSY;

  // start non-blocking MPI_IRecv to receive status
  recvStatus();

  // only return true if task successfully send to pgroup
  return true;
}


void ProcessGroupManager::recvStatus(){
  // start non-blocking call to receive status
  MPI_Irecv(&status_, 1, MPI_INT, pgroupRootID_, statusTag, theMPISystem()->getGlobalComm(),
              &statusRequest_);
}


bool ProcessGroupManager::recoverCommunicators(){
  assert( status_ == PROCESS_GROUP_WAIT );

  // send signal to pgroup
  SignalType signal = RECOVER_COMM;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());

  return true;
}

bool ProcessGroupManager::searchSDC( SDCMethodType method ){

  if (status_ != PROCESS_GROUP_WAIT)
    return false;

  // send signal to pgroup
  SignalType signal = SEARCH_SDC;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());

  // send method to pgroup
  MPI_Send(&method, 1, MPI_INT, pgroupRootID_, infoTag, theMPISystem()->getGlobalComm());

  status_ = PROCESS_GROUP_BUSY;

  // start non-blocking MPI_IRecv to receive status
  recvStatus();

  return true;
}
} /* namespace combigrid */
