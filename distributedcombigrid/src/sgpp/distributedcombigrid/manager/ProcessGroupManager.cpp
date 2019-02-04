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


template <typename FG_ELEMENT>
bool ProcessGroupManager::exchangeCommonSubspacesThirdLevel(const ThirdLevelUtils& thirdLevel_,
                                                                  CombiParameters& params) {
  // can only send sync signal when in wait state
  assert(status_ == PROCESS_GROUP_WAIT);

  SignalType signal = EXCHANGE_COMMON_SUBSPACES_THIRD_LEVEL;
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, signalTag, theMPISystem()->getGlobalComm());

  // set status
  status_ = PROCESS_GROUP_BUSY;

  const std::vector<LevelVector> commonSS = params.getCommonSubspaces();
  const size_t numCommonSS = commonSS.size();
  const std::vector<CommunicatorType>& thirdLevelComms = theMPISystem()->getThirdLevelComms();

  MPI_Datatype dtype = abstraction::getMPIDatatype(
                         abstraction::getabstractionDataType<FG_ELEMENT>());
  for (size_t p = 0; p < theMPISystem()->getNumProcs(); p++) {
    const CommunicatorType& comm = thirdLevelComms[p];

    // receive sizes
    std::vector<int> commonSSPartSizes(numCommonSS);
    MPI_Recv(commonSSPartSizes.data(), static_cast<int>(numCommonSS), MPI_INT, 
        0, 0, comm, MPI_STATUS_IGNORE);

    // create buffer
    size_t buffsize = 0;
    for (size_t ss = 0; ss < numCommonSS; ss++)
      buffsize += static_cast<size_t>(commonSSPartSizes[ss]);
    std::vector<FG_ELEMENT> commonSSPart(buffsize, FG_ELEMENT(0));

    // receive subspace data from worker
    size_t stride = 0;
    for (size_t ss = 0; ss < numCommonSS; ss++) {
      MPI_Recv(commonSSPart.data()+stride, commonSSPartSizes[ss], dtype, 0, 0,
          comm, MPI_STATUS_IGNORE);
      stride += static_cast<size_t>(commonSSPartSizes[ss]);
    }

    // send subspace data to remote
    thirdLevel_.sendCommonSubspaces(commonSSPart);
    // receive combined data from remote
    thirdLevel_.receiveCommonSubspaces(commonSSPart);

    // send combined subspace data to worker
    stride = 0;
    for (size_t ss = 0; ss < numCommonSS; ss++) {
      MPI_Send(commonSSPart.data()+stride, dtype, 0, 0, comm);
      stride += static_cast<size_t>(commonSSPartSizes[ss]);
    }
  }

  // start non-blocking MPI_IRecv to receive status
  recvStatus();

  return true;
}

template<typename FG_ELEMENT>
bool ProcessGroupManager::combineUniformThirdLevel(const ThirdLevelUtils& thirdLevel_,
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
  for (size_t p = 0; p < theMPISystem()->getNumProcs(); p++) {
    const CommunicatorType& comm = thirdLevelComms[p];
    assert(comm != MPI_COMM_NULL);

    // receive sizes
    std::vector<int> commonSSPartSizes(numCommonSS);
    MPI_Recv(commonSSPartSizes.data(), static_cast<int>(numCommonSS), MPI_INT, 
        0, 0, comm, MPI_STATUS_IGNORE);

    // create buffer
    size_t buffsize = 0;
    for (size_t ss = 0; ss < numCommonSS; ss++)
      buffsize += static_cast<size_t>(commonSSPartSizes[ss]);
    std::vector<FG_ELEMENT> commonSSPart(buffsize, FG_ELEMENT(0));

    // receive data from remote
    thirdLevel_.receiveCommonSubspaces(commonSSPart, commonSSPartSizes[p]);

    // combine subspaces sequentially
    size_t stride = 0;
    for (size_t ss = 0; ss < numCommonSS; ss++) {
      MPI_Allreduce(MPI_IN_PLACE, commonSSPart.data()+stride, commonSSPartSizes[ss], dtype, MPI_SUM, comm);
      stride += static_cast<size_t>(commonSSPartSizes[ss]);
    }

    // send data back to remote
    thirdLevel_.sendCommonSubspaces(commonSSPart, commonSSPartSizes[p]);
  }

  // start non-blocking MPI_IRecv to receive status
  recvStatus();

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

bool ProcessGroupManager::updateCombiParameters(CombiParameters& params) {
  // can only send sync signal when in wait state
  assert(status_ == PROCESS_GROUP_WAIT);

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


void ProcessGroupManager::recvStatus(){
  // start non-blocking call to receive status
  MPI_Irecv(&status_, 1, MPI_INT, pgroupRootID_, statusTag, theMPISystem()->getGlobalComm(),
              &statusRequest_);
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


} /* namespace combigrid */
