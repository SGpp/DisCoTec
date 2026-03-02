#include "manager/ProcessGroupManager.hpp"

#include <complex>

#include "manager/CombiParameters.hpp"
#include "mpi/MPIUtils.hpp"
#include "mpi_fault_simulator/MPI-FT.h"

namespace combigrid {
template <typename CombiDataType>
ProcessGroupManager<CombiDataType>::ProcessGroupManager(RankType pgroupRootID)
    : pgroupRootID_(pgroupRootID),
      status_(PROCESS_GROUP_WAIT),
      statusRequest_(MPI_Request()),
      statusRequestFT_(nullptr) {}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::runfirst(Task<CombiDataType>* t) {
  return storeTaskReferenceAndSendTaskToProcessGroup(t, RUN_FIRST);
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::storeTaskReferenceAndSendTaskToProcessGroup(
    Task<CombiDataType>* t, SignalType signal) {
  // first check status
  // tying to add a task to a busy group is an invalid operation
  // and should be avoided
  if (status_ != PROCESS_GROUP_WAIT) return false;

  storeTaskReference(t);
  return sendTaskToProcessGroup(t, signal);
}

template <typename CombiDataType>
void ProcessGroupManager<CombiDataType>::storeTaskReference(Task<CombiDataType>* t) {
  // add task to list of tasks managed by this pgroup
  tasks_.push_back(t);
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::sendTaskToProcessGroup(Task<CombiDataType>* t,
                                                                SignalType signal) {
  // send signal to pgroup
  sendSignalToProcessGroup(signal);

  // send task
  Task<CombiDataType>::send(t, pgroupRootID_, theMPISystem()->getGlobalComm());

  setProcessGroupBusyAndReceive();

  // only return true if task successfully sent to pgroup
  return true;
}

template <typename CombiDataType>
void ProcessGroupManager<CombiDataType>::sendSignalAndReceive(SignalType signal) {
  sendSignalToProcessGroup(signal);
  setProcessGroupBusyAndReceive();
}

template <typename CombiDataType>
void ProcessGroupManager<CombiDataType>::sendSignalToProcessGroup(SignalType signal) {
  MPI_Send(&signal, 1, MPI_INT, pgroupRootID_, TRANSFER_SIGNAL_TAG,
           theMPISystem()->getGlobalComm());
}

template <typename CombiDataType>
inline void ProcessGroupManager<CombiDataType>::setProcessGroupBusyAndReceive() {
  // set status
  status_ = PROCESS_GROUP_BUSY;

  // start non-blocking MPI_IRecv to receive status
  recvStatus();
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::runnext() {
  // first check status
  // trying to send a command to a busy group is an invalid operation
  // and should be avoided
  assert(status_ == PROCESS_GROUP_WAIT);

  if (tasks_.size() == 0) return false;

  sendSignalAndReceive(RUN_NEXT);

  return true;
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::exit() {
  // can only send exit signal when in wait state
  if (status_ != PROCESS_GROUP_WAIT) return false;

  sendSignalAndReceive(EXIT);
  return true;
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::combine() {
  // can only send sync signal when in wait state
  assert(status_ == PROCESS_GROUP_WAIT);

  sendSignalAndReceive(COMBINE);

  return true;
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::initDsgus() {
  // can only send sync signal when in wait state
  assert(status_ == PROCESS_GROUP_WAIT);

  sendSignalAndReceive(INIT_DSGUS);

  return true;
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::combineThirdLevel(const ThirdLevelUtils& thirdLevel,
                                                           CombiParameters& params,
                                                           bool isSendingFirst) {
  // can only send sync signal when in wait state
  assert(status_ == PROCESS_GROUP_WAIT);

  sendSignalAndReceive(COMBINE_THIRD_LEVEL);

  exchangeDsgus(thirdLevel, params, isSendingFirst);

  return true;
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::combineThirdLevelFileBased(
    const std::string& filenamePrefixToWrite, const std::string& writeCompleteTokenFileName,
    const std::string& filenamePrefixToRead, const std::string& startReadingTokenFileName) {
  assert(waitStatus() == PROCESS_GROUP_WAIT);
  sendSignalAndReceive(COMBINE_THIRD_LEVEL_FILE);

  // send filenames
  MPIUtils::sendClass(&filenamePrefixToWrite, this->pgroupRootID_, theMPISystem()->getGlobalComm());
  MPIUtils::sendClass(&writeCompleteTokenFileName, this->pgroupRootID_,
                      theMPISystem()->getGlobalComm());
  MPIUtils::sendClass(&filenamePrefixToRead, this->pgroupRootID_, theMPISystem()->getGlobalComm());
  MPIUtils::sendClass(&startReadingTokenFileName, this->pgroupRootID_,
                      theMPISystem()->getGlobalComm());
  return this->waitStatus() == PROCESS_GROUP_WAIT;
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::combineThirdLevelFileBasedWrite(
    const std::string& filenamePrefixToWrite, const std::string& writeCompleteTokenFileName) {
  assert(waitStatus() == PROCESS_GROUP_WAIT);
  sendSignalAndReceive(COMBINE_WRITE_DSGS);

  // send filenames
  MPIUtils::sendClass(&filenamePrefixToWrite, this->pgroupRootID_, theMPISystem()->getGlobalComm());
  MPIUtils::sendClass(&writeCompleteTokenFileName, this->pgroupRootID_,
                      theMPISystem()->getGlobalComm());
  auto waiting = this->waitStatus() == PROCESS_GROUP_WAIT;
  assert(waiting);
  return waiting;
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::combineThirdLevelFileBasedReadReduce(
    const std::string& filenamePrefixToRead, const std::string& startReadingTokenFileName) {
  assert(waitStatus() == PROCESS_GROUP_WAIT);
  sendSignalAndReceive(COMBINE_READ_DSGS_AND_REDUCE);

  // send filenames
  MPIUtils::sendClass(&filenamePrefixToRead, this->pgroupRootID_, theMPISystem()->getGlobalComm());
  MPIUtils::sendClass(&startReadingTokenFileName, this->pgroupRootID_,
                      theMPISystem()->getGlobalComm());
  auto waiting = this->waitStatus() == PROCESS_GROUP_WAIT;
  assert(waiting);
  return waiting;
}

template <typename CombiDataType>
void recvDsguFromWorker(std::vector<CombiDataType>& dsguData, RankType r, CommunicatorType comm) {
  MPI_Datatype dataType = getMPIDatatype(abstraction::getabstractionDataType<CombiDataType>());
  auto dsguSize = dsguData.size();
  if (dsguSize == 0) {
    throw std::runtime_error("dsguData is empty");
  }
  size_t sentRecvd = 0;
  while ((dsguSize - sentRecvd) / INT_MAX > 0) {
    MPI_Recv(dsguData.data() + sentRecvd, (int)INT_MAX, dataType, r, TRANSFER_DSGU_DATA_TAG, comm,
             MPI_STATUS_IGNORE);
    sentRecvd += INT_MAX;
  }
  MPI_Recv(dsguData.data() + sentRecvd, (int)(dsguSize - sentRecvd), dataType, r,
           TRANSFER_DSGU_DATA_TAG, comm, MPI_STATUS_IGNORE);
}

template <typename CombiDataType>
void sendDsguToWorker(std::vector<CombiDataType>& dsguData, RankType r, CommunicatorType comm) {
  MPI_Datatype dataType = getMPIDatatype(abstraction::getabstractionDataType<CombiDataType>());
  auto dsguSize = dsguData.size();
  if (dsguSize == 0) {
    throw std::runtime_error("dsguData is empty");
  }
  size_t sentRecvd = 0;
  while ((dsguSize - sentRecvd) / INT_MAX > 0) {
    MPI_Send(dsguData.data() + sentRecvd, (int)INT_MAX, dataType, r, TRANSFER_DSGU_DATA_TAG, comm);
    sentRecvd += INT_MAX;
  }
  MPI_Send(dsguData.data() + sentRecvd, (int)(dsguSize - sentRecvd), dataType, r,
           TRANSFER_DSGU_DATA_TAG, comm);
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::pretendCombineThirdLevelForWorkers(
    CombiParameters& params) {
  // can only send sync signal when in wait state
  assert(status_ == PROCESS_GROUP_WAIT);

  sendSignalAndReceive(COMBINE_THIRD_LEVEL);

  const std::vector<CommunicatorType>& thirdLevelComms = theMPISystem()->getThirdLevelComms();
  assert(theMPISystem()->getNumGroups() == thirdLevelComms.size() &&
         "initialisation of third level communicator failed");
  const CommunicatorType& comm = thirdLevelComms[params.getThirdLevelPG()];

  // exchange dsgus
  auto numGrids = params.getNumGrids();
  std::vector<CombiDataType> dsguData;
  for (size_t g = 0; g < numGrids; g++) {
    for (RankType p = 0; p < (RankType)theMPISystem()->getNumProcs(); p++) {
      // we assume here that all dsgus have the same size otherwise size collection must change
      size_t dsguSize = (size_t)(dsguDataSizePerWorker_[(size_t)p] / numGrids);
      assert(dsguSize > 0);

      // recv dsgu from worker
      dsguData.resize(dsguSize);
      recvDsguFromWorker(dsguData, p, comm);

      // send back to worker immediately
      sendDsguToWorker(dsguData, p, comm);
    }
  }

  return true;
}

/*
 * Differs from third level reduce since we have enough memory space to collect
 * the subspace sizes from all dsgs of all procs in the third level pg in a single
 * MPI_Gather call.
 */
template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::reduceLocalAndRemoteSubspaceSizes(
    const ThirdLevelUtils& thirdLevel, bool isSendingFirst, bool thirdLevelExtraSparseGrid) {
  // tell workers to perform reduce
  if (thirdLevelExtraSparseGrid) {
    sendSignalAndReceive(REDUCE_SUBSPACE_SIZES_TL_AND_ALLOCATE_EXTRA_SG);
  } else {
    sendSignalAndReceive(REDUCE_SUBSPACE_SIZES_TL);
  }

  // prepare buffers
  std::vector<SubspaceSizeType> sendBuff;
  std::vector<SubspaceSizeType> recvBuff;
  size_t buffSize;
  std::vector<int> numSubspacesPerWorker;

  // gather subspace sizes from workers
  collectSubspaceSizes(thirdLevel, sendBuff, buffSize, numSubspacesPerWorker);
  recvBuff.resize(buffSize);

  /* TODO can be easily parallelized by removing the condition and call send and
     receive in a separate thread*/
  if (isSendingFirst) {
    // send subspace sizes to remote
    thirdLevel.sendData(sendBuff.data(), buffSize);
    // receive remote subspace sizes
    thirdLevel.recvData(recvBuff.data(), buffSize);
  } else {
    // receive remote subspace sizes
    thirdLevel.recvData(recvBuff.data(), buffSize);
    // send subspace sizes to remote
    thirdLevel.sendData(sendBuff.data(), buffSize);
  }

  // set accumulated dsgu sizes per worker
  formerDsguDataSizePerWorker_.resize(numSubspacesPerWorker.size());
  auto from = sendBuff.begin();
  auto to = sendBuff.begin();
  assert(numSubspacesPerWorker.size() > 0);
  for (size_t w = 0; w < numSubspacesPerWorker.size(); ++w) {
    std::advance(to, numSubspacesPerWorker[w]);
    formerDsguDataSizePerWorker_[w] = std::accumulate(from, to, 0);
    from = to;
  }
  assert(to == sendBuff.end());

  if (thirdLevelExtraSparseGrid) {
    // perform min reduce
    for (size_t i = 0; i < buffSize; ++i) {
      assert(recvBuff[i] == sendBuff[i] || recvBuff[i] == 0 || sendBuff[i] == 0);
      sendBuff[i] = std::min(sendBuff[i], recvBuff[i]);
    }
  } else {
    // perform max reduce
    for (size_t i = 0; i < buffSize; ++i) {
      assert(recvBuff[i] == sendBuff[i] || recvBuff[i] == 0 || sendBuff[i] == 0);
      sendBuff[i] = std::max(sendBuff[i], recvBuff[i]);
    }
  }

  // scatter data back to workers
  distributeSubspaceSizes(thirdLevel, sendBuff, buffSize, numSubspacesPerWorker);

  // set accumulated dsgu sizes per worker
  dsguDataSizePerWorker_.resize(numSubspacesPerWorker.size());
  from = sendBuff.begin();
  to = sendBuff.begin();
  assert(numSubspacesPerWorker.size() > 0);
  for (size_t w = 0; w < numSubspacesPerWorker.size(); ++w) {
    std::advance(to, numSubspacesPerWorker[w]);
    dsguDataSizePerWorker_[w] = std::accumulate(from, to, 0);
    from = to;
  }
  assert(to == sendBuff.end());
  return true;
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::pretendReduceLocalAndRemoteSubspaceSizes(
    const ThirdLevelUtils& thirdLevel) {
  sendSignalAndReceive(REDUCE_SUBSPACE_SIZES_TL);

  // prepare buffers
  std::vector<SubspaceSizeType> sendBuff;
  std::vector<SubspaceSizeType> recvBuff;
  size_t buffSize;
  std::vector<int> numSubspacesPerWorker;

  // gather subspace sizes from workers
  collectSubspaceSizes(thirdLevel, sendBuff, buffSize, numSubspacesPerWorker);
  recvBuff.resize(buffSize);

  // don't send subspace sizes to remote
  // don't receive remote subspace sizes
  // instead, just return zeros to process group
  std::memset(recvBuff.data(), 0, recvBuff.size() * sizeof(SubspaceSizeType));

  // set accumulated dsgu sizes per worker
  formerDsguDataSizePerWorker_.resize(numSubspacesPerWorker.size());
  auto from = sendBuff.begin();
  auto to = sendBuff.begin();
  assert(numSubspacesPerWorker.size() > 0);
  for (size_t w = 0; w < numSubspacesPerWorker.size(); ++w) {
    std::advance(to, numSubspacesPerWorker[w]);
    formerDsguDataSizePerWorker_[w] = std::accumulate(from, to, 0);
    from = to;
  }
  assert(to == sendBuff.end());
  for ([[maybe_unused]] const auto& dataSize : formerDsguDataSizePerWorker_) {
    assert(dataSize > 0);
  }

  // perform max reduce
  for (size_t i = 0; i < buffSize; ++i) sendBuff[i] = std::max(sendBuff[i], recvBuff[i]);

  // scatter data back to workers
  distributeSubspaceSizes(thirdLevel, sendBuff, buffSize, numSubspacesPerWorker);

  // set accumulated dsgu sizes per worker
  dsguDataSizePerWorker_.resize(numSubspacesPerWorker.size());
  from = sendBuff.begin();
  to = sendBuff.begin();
  assert(numSubspacesPerWorker.size() > 0);
  for (size_t w = 0; w < numSubspacesPerWorker.size(); ++w) {
    std::advance(to, numSubspacesPerWorker[w]);
    dsguDataSizePerWorker_[w] = std::accumulate(from, to, 0);
    from = to;
  }
  assert(to == sendBuff.end());
  assert(waitStatus() == PROCESS_GROUP_WAIT);
  return true;
}

template <typename CombiDataType>
void ProcessGroupManager<CombiDataType>::exchangeDsgus(const ThirdLevelUtils& thirdLevel,
                                                       CombiParameters& params,
                                                       bool isSendingFirst) {
  const std::vector<CommunicatorType>& thirdLevelComms = theMPISystem()->getThirdLevelComms();
  assert(theMPISystem()->getNumGroups() == thirdLevelComms.size() &&
         "initialisation of third level communicator failed");
  const CommunicatorType& comm = thirdLevelComms[params.getThirdLevelPG()];

  // exchange dsgus
  auto numGrids = params.getNumGrids();
  std::vector<CombiDataType> dsguData;
  for (size_t g = 0; g < numGrids; g++) {
    for (RankType p = 0; p < (RankType)theMPISystem()->getNumProcs(); p++) {
      // we assume here that all dsgus have the same size otherwise size collection must change
      size_t dsguSize = (size_t)(dsguDataSizePerWorker_[(size_t)p] / numGrids);

      // recv dsgu from worker
      dsguData.resize(dsguSize);
      recvDsguFromWorker(dsguData, p, comm);

      if (isSendingFirst) {
        // send dsgu to remote
        thirdLevel.sendData(dsguData.data(), dsguSize);
        // recv combined dsgu from remote
        thirdLevel.recvData(dsguData.data(), dsguSize);
      } else {
        // recv and combine dsgu from remote
        thirdLevel.recvAndAddToData(dsguData.data(), dsguSize);
        // send combined solution to remote
        thirdLevel.sendData(dsguData.data(), dsguSize);
      }
      // send to worker
      sendDsguToWorker(dsguData, p, comm);
    }
  }
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::collectSubspaceSizes(
    const ThirdLevelUtils& thirdLevel, std::vector<SubspaceSizeType>& buff, size_t& buffSize,
    std::vector<int>& numSubspacesPerWorker) {
  // prepare args of MPI_Gather
  const CommunicatorType& comm = theMPISystem()->getThirdLevelComms()[(size_t)pgroupRootID_];
  size_t nprocs = theMPISystem()->getNumProcs();
  std::vector<int> recvCounts(nprocs + 1);  // includes master
  RankType thirdLevelManagerRank = theMPISystem()->getThirdLevelManagerRank();
  int dummy = 0;

  // gather number of subspaces in all dsgus per worker for upcoming MPI_Gatherv
  // for now all workers should have the same number of subspaces
  MPI_Gather(&dummy, 1, MPI_INT, recvCounts.data(), (int)1, MPI_INT, thirdLevelManagerRank, comm);

  buffSize = std::accumulate(recvCounts.begin(), recvCounts.end(), 0U);

  std::vector<SubspaceSizeType> mdBuff(buffSize);
  assert(buffSize < INT_MAX &&
         "bufSize is larger than what we can send in a "
         "single mpi call");

  // prepare displacements for MPI_Gatherv
  std::vector<int> displacements(nprocs + 1);  // includes master
  int disp = 0;
  for (size_t i = 0; i < displacements.size(); ++i) {
    displacements[i] = disp;
    disp += recvCounts[i];
  }
  // perform gather of subspace sizes
  MPI_Datatype dataType = getMPIDatatype(abstraction::getabstractionDataType<SubspaceSizeType>());
  MPI_Gatherv(&dummy, 0, dataType, mdBuff.data(), recvCounts.data(), displacements.data(), dataType,
              thirdLevelManagerRank, comm);

  // remove master
  numSubspacesPerWorker = recvCounts;
  numSubspacesPerWorker.pop_back();

  // create machine independent buffer
  buff.resize(buffSize);
  buff.assign(mdBuff.begin(), mdBuff.end());

  return true;
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::distributeSubspaceSizes(
    const ThirdLevelUtils& thirdLevel, const std::vector<SubspaceSizeType>& buff, size_t buffSize,
    const std::vector<int>& numSubspacesPerWorker) {
  // prepare args of MPI_Scatterv
  const CommunicatorType& comm = theMPISystem()->getThirdLevelComms()[(size_t)pgroupRootID_];
  RankType thirdLevelManagerRank = theMPISystem()->getThirdLevelManagerRank();
  std::vector<int> sendCounts(numSubspacesPerWorker);
  sendCounts.push_back(0);  // append manager
  std::vector<int> displacements(sendCounts.size());
  int disp = 0;
  for (size_t i = 0; i < displacements.size(); ++i) {
    displacements[i] = disp;
    disp += sendCounts[i];
  }

  // create machine dependent buffer
  std::vector<SubspaceSizeType> mdBuff(buffSize);
  mdBuff.assign(buff.begin(), buff.end());

  // perform scatter of subspace sizes
  MPI_Datatype dataType = getMPIDatatype(abstraction::getabstractionDataType<SubspaceSizeType>());
  MPI_Scatterv(mdBuff.data(), sendCounts.data(), displacements.data(), dataType, nullptr, 0,
               dataType, thirdLevelManagerRank, comm);

  return true;
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::combineSystemWide() {
  // can only send sync signal when in wait state
  assert(status_ == PROCESS_GROUP_WAIT);

  sendSignalAndReceive(COMBINE_LOCAL_AND_GLOBAL);
  return true;
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::waitForThirdLevelCombiResult() {
  // can only send sync signal when in wait state
  assert(this->waitStatus() == PROCESS_GROUP_WAIT);

  sendSignalAndReceive(WAIT_FOR_TL_COMBI_RESULT);
  return true;
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::waitForThirdLevelSizeUpdate() {
  // can only send sync signal when in wait state
  assert(status_ == PROCESS_GROUP_WAIT);

  sendSignalAndReceive(WAIT_FOR_TL_SIZE_UPDATE);
  return true;
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::updateCombiParameters(CombiParameters& params) {
  // can only send sync signal when in wait state
  assert(status_ == PROCESS_GROUP_WAIT);

  sendSignalToProcessGroup(UPDATE_COMBI_PARAMETERS);

  // send combiparameters
  MPIUtils::sendClass(&params, pgroupRootID_, theMPISystem()->getGlobalComm());

  setProcessGroupBusyAndReceive();
  return true;
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::addTask(Task<CombiDataType>* t) {
  return storeTaskReferenceAndSendTaskToProcessGroup(t, ADD_TASK);
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::refreshTask(Task<CombiDataType>* t) {
  // first check status
  // tying to add a task to a busy group is an invalid operation
  // and should be avoided
  if (status_ != PROCESS_GROUP_WAIT) {
    std::cout << "refreshing failed! \n";
    return false;
  }

  sendSignalToProcessGroup(ADD_TASK);

  // send task
  Task<CombiDataType>::send(t, pgroupRootID_, theMPISystem()->getGlobalComm());

  setProcessGroupBusyAndReceive();

  // only return true if task successfully send to pgroup
  return true;
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::resetTasksWorker() {
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

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::recompute(Task<CombiDataType>* t) {
  storeTaskReferenceAndSendTaskToProcessGroup(t, RECOMPUTE);

  // only return true if task successfully send to pgroup
  return true;
}

void sendLevelVector(const LevelVector& leval, RankType pgroupRootID) {
  std::vector<int> tmp(leval.begin(), leval.end());
  MPI_Send(&tmp[0], static_cast<int>(tmp.size()), MPI_INT, pgroupRootID, TRANSFER_LEVAL_TAG,
           theMPISystem()->getGlobalComm());
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::parallelEval(const LevelVector& leval,
                                                      const std::string& filename) {
  // can only send sync signal when in wait state, so check first
  assert(status_ == PROCESS_GROUP_WAIT);

  sendSignalToProcessGroup(PARALLEL_EVAL);

  sendLevelVector(leval, pgroupRootID_);

  // send filename
  MPIUtils::sendClass(&filename, pgroupRootID_, theMPISystem()->getGlobalComm());

  setProcessGroupBusyAndReceive();

  return true;
}

template <typename CombiDataType>
void ProcessGroupManager<CombiDataType>::writeSparseGridMinMaxCoefficients(
    const std::string& filename) {
  this->sendSignalToProcessGroup(WRITE_DSG_MINMAX_COEFFICIENTS);

  // send filename
  MPIUtils::sendClass(&filename, this->pgroupRootID_, theMPISystem()->getGlobalComm());

  this->setProcessGroupBusyAndReceive();
}

template <typename CombiDataType>
void ProcessGroupManager<CombiDataType>::doDiagnostics(size_t taskID) {
  [[maybe_unused]] auto status = waitStatus();
  assert(status == PROCESS_GROUP_WAIT);
  for (auto task : tasks_) {
    if (task->getID() == taskID) {
      sendSignalToProcessGroup(DO_DIAGNOSTICS);
      // send task ID to do postprocessing on
      MPI_Send(&taskID, 1,
               abstraction::getMPIDatatype(abstraction::getabstractionDataType<decltype(taskID)>()),
               this->pgroupRootID_, 0, theMPISystem()->getGlobalComm());
      return;
    }
  }
  assert(false && "Task was not present in this process group");
}

std::vector<double> receiveThreeNorms(RankType pgroupRootID) {
  std::vector<double> norms;
  for (int i = 0; i < 3; ++i) {
    double recvbuf;

    MPI_Recv(&recvbuf, 1, MPI_DOUBLE, pgroupRootID, TRANSFER_NORM_TAG,
             theMPISystem()->getGlobalComm(), MPI_STATUS_IGNORE);

    norms.push_back(recvbuf);
  }
  return norms;
}

template <typename CombiDataType>
void ProcessGroupManager<CombiDataType>::getLpNorms(int p, std::map<size_t, double>& norms) {
  SignalType signal = GET_L1_NORM;
  if (p == 2) {
    signal = GET_L2_NORM;
  } else if (p == 1) {
    signal = GET_L1_NORM;
  } else if (p == 0) {
    signal = GET_MAX_NORM;
  } else {
    assert(false && "please implement signal for this norm");
  }

  this->sendSignalToProcessGroup(signal);

  std::vector<double> recvbuf;
  auto numTasks = this->getTaskContainer().size();
  recvbuf.resize(numTasks);

  MPI_Recv(recvbuf.data(), static_cast<int>(numTasks), MPI_DOUBLE, pgroupRootID_, TRANSFER_NORM_TAG,
           theMPISystem()->getGlobalComm(), MPI_STATUS_IGNORE);

  for (size_t j = 0; j < numTasks; ++j) {
    norms[this->getTaskContainer()[j]->getID()] = recvbuf[j];
  }

  this->setProcessGroupBusyAndReceive();
}

template <typename CombiDataType>
std::vector<double> ProcessGroupManager<CombiDataType>::evalAnalyticalOnDFG(
    const LevelVector& leval) {
  sendSignalToProcessGroup(EVAL_ANALYTICAL_NORM);
  sendLevelVector(leval, pgroupRootID_);

  auto norms = receiveThreeNorms(pgroupRootID_);
  setProcessGroupBusyAndReceive();
  return norms;
}

template <typename CombiDataType>
std::vector<double> ProcessGroupManager<CombiDataType>::evalErrorOnDFG(const LevelVector& leval) {
  sendSignalToProcessGroup(EVAL_ERROR_NORM);
  sendLevelVector(leval, pgroupRootID_);

  auto norms = receiveThreeNorms(pgroupRootID_);
  setProcessGroupBusyAndReceive();
  return norms;
}

template <typename CombiDataType>
void ProcessGroupManager<CombiDataType>::interpolateValues(
    const std::vector<real>& interpolationCoordsSerial, std::vector<CombiDataType>& values,
    MPI_Request* request, const std::string& filenamePrefix) {
  assert(interpolationCoordsSerial.size() < static_cast<size_t>(std::numeric_limits<int>::max()) &&
         "needs chunking!");
  for ([[maybe_unused]] const auto& coord : interpolationCoordsSerial) {
    assert(coord >= 0.0 && coord <= 1.0);
  }
  // if no request was passed, assume we are not waiting for a reply from that group
  if (request == nullptr && values.empty() && filenamePrefix == "") {
    sendSignalToProcessGroup(INTERPOLATE_VALUES);
  } else if (request != nullptr && !values.empty() && filenamePrefix == "") {
    sendSignalToProcessGroup(INTERPOLATE_VALUES_AND_SEND_BACK);
  } else if (request == nullptr && values.empty() && filenamePrefix != "") {
    sendSignalToProcessGroup(INTERPOLATE_VALUES_AND_WRITE_SINGLE_FILE);
    MPIUtils::sendClass(&filenamePrefix, pgroupRootID_, theMPISystem()->getGlobalComm());
  } else {
    throw std::runtime_error("invalid combination of arguments");
  }
  MPI_Request dummyRequest;
  // send interpolation coordinates to group
  MPI_Isend(interpolationCoordsSerial.data(), static_cast<int>(interpolationCoordsSerial.size()),
            abstraction::getMPIDatatype(abstraction::getabstractionDataType<real>()), pgroupRootID_,
            TRANSFER_INTERPOLATION_TAG, theMPISystem()->getGlobalComm(), &dummyRequest);
  if (!values.empty()) {
    MPI_Request_free(&dummyRequest);
  } else {
    MPI_Wait(&dummyRequest, MPI_STATUS_IGNORE);
  }
  if (request != nullptr) {
    MPI_Irecv(values.data(), static_cast<int>(values.size()),
              abstraction::getMPIDatatype(abstraction::getabstractionDataType<CombiDataType>()),
              pgroupRootID_, TRANSFER_INTERPOLATION_TAG, theMPISystem()->getGlobalComm(), request);
  }
  setProcessGroupBusyAndReceive();
}

template <typename CombiDataType>
void ProcessGroupManager<CombiDataType>::writeInterpolatedValuesPerGrid(
    const std::vector<real>& interpolationCoordsSerial, const std::string& filenamePrefix) {
  sendSignalToProcessGroup(WRITE_INTERPOLATED_VALUES_PER_GRID);
  // send filename prefix to group
  MPIUtils::sendClass(&filenamePrefix, pgroupRootID_, theMPISystem()->getGlobalComm());
  assert(interpolationCoordsSerial.size() < static_cast<size_t>(std::numeric_limits<int>::max()) &&
         "needs chunking!");
  MPI_Send(interpolationCoordsSerial.data(), static_cast<int>(interpolationCoordsSerial.size()),
           abstraction::getMPIDatatype(abstraction::getabstractionDataType<real>()), pgroupRootID_,
           TRANSFER_INTERPOLATION_TAG, theMPISystem()->getGlobalComm());
  setProcessGroupBusyAndReceive();
  assert(waitStatus() == PROCESS_GROUP_WAIT);
}

template <typename CombiDataType>
void ProcessGroupManager<CombiDataType>::recvStatus() {
  // start non-blocking call to receive status
  if (ENABLE_FT) {
    simft::Sim_FT_MPI_Irecv(&status_, 1, MPI_INT, pgroupRootID_, TRANSFER_STATUS_TAG,
                            theMPISystem()->getGlobalCommFT(), &statusRequestFT_);
  } else {
    MPI_Irecv(&status_, 1, MPI_INT, pgroupRootID_, TRANSFER_STATUS_TAG,
              theMPISystem()->getGlobalComm(), &statusRequest_);
  }
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::recoverCommunicators() {
  assert(status_ == PROCESS_GROUP_WAIT);

  sendSignalToProcessGroup(RECOVER_COMM);

  return true;
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::rescheduleAddTask(Task<CombiDataType>* task) {
  return storeTaskReferenceAndSendTaskToProcessGroup(task, RESCHEDULE_ADD_TASK);
}

template <typename CombiDataType>
Task<CombiDataType>* ProcessGroupManager<CombiDataType>::rescheduleRemoveTask(
    const LevelVector& lvlVec) {
  for (typename std::vector<Task<CombiDataType>*>::size_type i = 0; i < this->tasks_.size(); ++i) {
    Task<CombiDataType>* currentTask = this->tasks_[i];
    if (currentTask->getLevelVector() == lvlVec) {
      // if the task has been found send remove signal and return the task
      Task<CombiDataType>* removedTask;
      auto taskID = currentTask->getID();
      sendSignalToProcessGroup(RESCHEDULE_REMOVE_TASK);
      MPI_Send(&taskID, 1,
               abstraction::getMPIDatatype(abstraction::getabstractionDataType<decltype(taskID)>()),
               this->pgroupRootID_, 0, theMPISystem()->getGlobalComm());
      Task<CombiDataType>::receive(&removedTask, this->pgroupRootID_,
                                   theMPISystem()->getGlobalComm());
      setProcessGroupBusyAndReceive();

      tasks_.erase(tasks_.begin() + i);
      delete currentTask;

      return removedTask;
    }
  }
  return nullptr;
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::writeCombigridsToVTKPlotFile() {
  // can only send sync signal when in wait state
  assert(waitStatus() == PROCESS_GROUP_WAIT);

  sendSignalAndReceive(WRITE_DFGS_TO_VTK);
  return true;
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::writeDSGsToDisk(const std::string& filenamePrefix) {
  assert(waitStatus() == PROCESS_GROUP_WAIT);
  sendSignalAndReceive(WRITE_DSGS_TO_DISK);
  MPIUtils::sendClass(&filenamePrefix, pgroupRootID_, theMPISystem()->getGlobalComm());
  return true;
}

template <typename CombiDataType>
bool ProcessGroupManager<CombiDataType>::readDSGsFromDisk(const std::string& filenamePrefix) {
  assert(waitStatus() == PROCESS_GROUP_WAIT);
  sendSignalAndReceive(READ_DSGS_FROM_DISK);
  MPIUtils::sendClass(&filenamePrefix, pgroupRootID_, theMPISystem()->getGlobalComm());
  return true;
}

// explicit instantiations
template class ProcessGroupManager<double>;
template class ProcessGroupManager<std::complex<float>>;

template void recvDsguFromWorker<double>(std::vector<double>&, RankType, CommunicatorType);
template void recvDsguFromWorker<std::complex<float>>(std::vector<std::complex<float>>&, RankType,
                                                      CommunicatorType);

template void sendDsguToWorker<double>(std::vector<double>&, RankType, CommunicatorType);
template void sendDsguToWorker<std::complex<float>>(std::vector<std::complex<float>>&, RankType,
                                                    CommunicatorType);

} /* namespace combigrid */
