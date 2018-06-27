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
ProcessGroupManager::ProcessGroupManager( RankType pgroupRootID ) :
        pgroupRootID_(pgroupRootID),
        status_(PROCESS_GROUP_WAIT),
        statusRequest_(MPI_Request()),
        statusRequestFT_(nullptr)
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
  // can only send sync signal when in wait state
  assert(status_ == PROCESS_GROUP_WAIT);

  SignalType signal = UPDATE_COMBI_PARAMETERS;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());
  std::cout << "sending class \n";
  // send combiparameters
  MPIUtils::sendClass(&params, pgroupRootID_, theMPISystem()->getGlobalComm());

  // set status
  status_ = PROCESS_GROUP_BUSY;

  // start non-blocking MPI_IRecv to receive status
  recvStatus();
  std::cout << "manager received status \n";
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

bool ProcessGroupManager::refreshTask( Task* t ) {
  // first check status
  // tying to add a task to a busy group is an invalid operation
  // and should be avoided
  if (status_ != PROCESS_GROUP_WAIT){
    std::cout << "refreshing failed! \n";
    return false;

  }

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

bool ProcessGroupManager::resetTasksWorker() {
  // first check status
  // tying to reset tasks of a busy group is an invalid operation
  // and should be avoided
  if (status_ != PROCESS_GROUP_WAIT){
    assert(false);
    //return false;
  }

  // add task to list of tasks managed by this pgroup
  //tasks_.clear(); we do not clear group manager tasks

  // send add task signal to pgroup
  SignalType signal = RESET_TASKS;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());


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

bool ProcessGroupManager::parallelEval( const LevelVector& leval,
                                        std::string& filename ) {
  // can only send sync signal when in wait state, so check first
  assert(status_ == PROCESS_GROUP_WAIT);

  // send signal
  SignalType signal = PARALLEL_EVAL;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());

  // send levelvector
  std::vector<int> tmp( leval.begin(), leval.end() );
  MPI_Send(&tmp[0], static_cast<int>(tmp.size()), MPI_INT, pgroupRootID_, 0,
           theMPISystem()->getGlobalComm());

  // send filename
  MPIUtils::sendClass( &filename, pgroupRootID_,
                       theMPISystem()->getGlobalComm() );

  // set status
  status_ = PROCESS_GROUP_BUSY;

  // start non-blocking MPI_IRecv to receive status
  recvStatus();

  return true;
}

void ProcessGroupManager::recvStatus(){
  // start non-blocking call to receive status
  if( ENABLE_FT){
    simft::Sim_FT_MPI_Irecv( &status_, 1, MPI_INT, pgroupRootID_, statusTag,
                             theMPISystem()->getGlobalCommFT(), &statusRequestFT_ );
  } else{
    MPI_Irecv(&status_, 1, MPI_INT, pgroupRootID_, statusTag, theMPISystem()->getGlobalComm(),
              &statusRequest_);
  }
}


bool ProcessGroupManager::recoverCommunicators(){
  assert( status_ == PROCESS_GROUP_WAIT );

  // send signal to pgroup
  SignalType signal = RECOVER_COMM;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());

  return true;
}

void ProcessGroupManager::startBestExpansion(){
	sendSignal(BEST_EXPANSION);
}

std::pair<double, LevelVector> ProcessGroupManager::getBestExpansion(DimType dim){

	double error = -1;
	LevelVector expansion (dim);
	MPI_Recv(&error, 1, MPI_DOUBLE, pgroupRootID_, 1234, theMPISystem()->getGlobalComm(), MPI_STATUS_IGNORE);
	MPI_Recv(expansion.data(), expansion.size(), MPI_LONG, pgroupRootID_, 1235, theMPISystem()->getGlobalComm(), MPI_STATUS_IGNORE);

	return std::make_pair(error, expansion);
}

void ProcessGroupManager::sendTaskToProc(const std::map<int, int>& taskToProc){
	sendSignal(TASK_TO_PROC);
	MPIUtils::sendClass(&taskToProc, pgroupRootID_, theMPISystem()->getGlobalComm());
}

void ProcessGroupManager::addExpansion(const LevelVector& expansion){
	constexpr int addExpansionTag = 1235;
	sendSignal(ADD_EXPANSION);
	MPI_Send(expansion.data(), expansion.size(), MPI_INT, pgroupRootID_, addExpansionTag, theMPISystem()->getGlobalComm());
}


} /* namespace combigrid */
