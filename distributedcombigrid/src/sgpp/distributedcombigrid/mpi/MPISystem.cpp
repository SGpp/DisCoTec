#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupManager.hpp"
#include "sgpp/distributedcombigrid/utils/Stats.hpp"

#include <iostream>

namespace {
using namespace combigrid;

std::string getMinMaxAvg(RankType r, int size, std::string timerName, bool isTimer, MPI_Comm comm) {
  double value;
  /*    if( isTimer )
        value = theStatsContainer()->getDuration( timerName );
      else
        value = theStatsContainer()->getValue( timerName );
  */
  double min, max, sum;
  MPI_Reduce(&value, &min, 1, MPI_DOUBLE, MPI_MIN, r, comm);
  MPI_Reduce(&value, &max, 1, MPI_DOUBLE, MPI_MAX, r, comm);
  MPI_Reduce(&value, &sum, 1, MPI_DOUBLE, MPI_SUM, r, comm);

  double avg = sum / static_cast<double>(size);

  std::stringstream ss;
  ss << timerName << "\t\t" << min << "\t" << max << "\t" << avg;

  return ss.str();
}
}  // anonymous namespace

namespace combigrid {

/*!\brief Constructor for the MPISystem class.
 //
 // \exception std::runtime_error MPI was not initialized
 //
 // Constructor for the MPI System. The default global communicator and local communicator is
 MPI_COMM_WORLD.
 // The total number of MPI processes and the rank of the MPI process in is determined from
 // the communicator. The default root process has rank zero.
 */
MPISystem::MPISystem()
    : initialized_(false),
      worldComm_(MPI_COMM_NULL),
      globalComm_(MPI_COMM_NULL),
      localComm_(MPI_COMM_NULL),
      globalReduceComm_(MPI_COMM_NULL),
      worldCommFT_(nullptr),
      globalCommFT_(nullptr),
      spareCommFT_(nullptr),
      localCommFT_(nullptr),
      globalReduceCommFT_(nullptr),
      worldRank_(MPI_UNDEFINED),
      globalRank_(MPI_UNDEFINED),
      localRank_(MPI_UNDEFINED),
      globalReduceRank_(MPI_UNDEFINED),
      managerRank_(MPI_UNDEFINED),
      managerRankWorld_(MPI_UNDEFINED),
      masterRank_(MPI_UNDEFINED),
      reusableRanks_(std::vector<RankType>(0)) {
  // check if MPI was initialized (e.g. by MPI_Init or similar)
  int mpiInitialized(0);
  MPI_Initialized(&mpiInitialized);
  assert(mpiInitialized && "MPI is not initialized! Call MPI_Init first.");
}

/*!\brief Destructor for the MPISystem class.
 */
MPISystem::~MPISystem() {
  // todo: the fault tolerant communicator are initialized with new -> delete
}

void MPISystem::initSystemConstants(size_t ngroup, size_t nprocs, CommunicatorType worldComm = MPI_COMM_WORLD, bool reusable = false) {
  assert(reusable || !initialized_ && "MPISystem already initialized!");
  
  ngroup_ = ngroup;
  nprocs_ = nprocs;

  /* init worldComm
   * the manager has highest rank here
   */
  worldComm_ = worldComm;  
  int worldSize = getCommSize(worldComm_);
  assert(worldSize == int(ngroup_ * nprocs_ + 1));

  worldRank_ = getCommRank(worldComm_);
  managerRankWorld_ = worldSize - 1;
  managerRankFT_ = managerRankWorld_;
  
  if (ENABLE_FT) {
    worldCommFT_ = simft::Sim_FT_MPI_COMM_WORLD;
    MPI_Comm worldCommdup;
    MPI_Comm_dup(worldComm_, &worldCommdup);
    createCommFT(&spareCommFT_, worldCommdup);
    // spareCommFT_ = simft::Sim_FT_MPI_COMM_WORLD;
  }
  // std::cout << "Global rank of root is" << worldCommFT_->Root_Rank << "\n";
  // worldCommFT_->Root_Rank = worldSize - 1;
}

void MPISystem::init(size_t ngroup, size_t nprocs) {
  initSystemConstants(ngroup, nprocs);

  /* init localComm
   * lcomm is the local communicator of its own process group for each worker process.
   * for manager, lcomm is a group which contains only manager process and can be ignored
   */
  initLocalComm();

  /* create global communicator which contains only the manager and the master
   * process of each process group
   * the master processes of the process groups are the processes which have
   * rank 0 in lcomm
   * this communicator is used for communication between master processes of the
   * process groups and the manager and the master processes to each other
   */
  initGlobalComm();

  initGlobalReduceCommm();

  initialized_ = true;
}

/*  here the local communicator has already been created by the application */
void MPISystem::init(size_t ngroup, size_t nprocs, CommunicatorType lcomm) {
  initSystemConstants(ngroup, nprocs);

  storeLocalComm(lcomm);

  /* create global communicator which contains only the manager and the master
   * process of each process group
   * the master processes of the process groups are the processes which have
   * rank 0 in lcomm
   * this communicator is used for communication between master processes of the
   * process groups and the manager and the master processes to each other
   */
  initGlobalComm();

  initGlobalReduceCommm();

  initialized_ = true;
}

/* overload for initialization with given wold communicator
 * this method can be called multiple times (needed for tests)
 */
void MPISystem::initWorldReusable(CommunicatorType wcomm, size_t ngroup, size_t nprocs) {
  
  initSystemConstants(ngroup, nprocs, wcomm, true);

  /* init localComm
   * lcomm is the local communicator of its own process group for each worker process.
   * for manager, lcomm is a group which contains only manager process and can be ignored
   */
  initLocalComm();

  /* create global communicator which contains only the manager and the master
   * process of each process group
   * the master processes of the process groups are the processes which have
   * rank 0 in lcomm
   * this communicator is used for communication between master processes of the
   * process groups and the manager and the master processes to each other
   */
  initGlobalComm();

  initGlobalReduceCommm();

  initialized_ = true;
}

void MPISystem::initLocalComm() {
  int color = worldRank_ / int(nprocs_);
  int key = worldRank_ - color * int(nprocs_);
  MPI_Comm_split(worldComm_, color, key, &localComm_);

  /* set group number in Stats. this is necessary for postprocessing */
  Stats::setAttribute("group", std::to_string(color));

  storeLocalComm();
}

void MPISystem::storeLocalComm(CommunicatorType lcomm_optional /*= MPI_COMM_NULL*/){
  /* store localComm
   * lcomm is the local communicator of its own process group for each worker process.
   * for manager, lcomm is a group which contains only manager process and can be ignored
   */
  // manager is not supposed to have a localComm
  if (worldRank_ == managerRankWorld_)
    localComm_ = MPI_COMM_NULL;
  else {
    // duplicate to localComm_ only if lcomm_optional was given, otherwise assume already initialized
    if(lcomm_optional != MPI_COMM_NULL)  MPI_Comm_dup(lcomm_optional, &localComm_);
    // todo: think through which side effects changing the master rank would have
    // in principle this does not have to be 0
    masterRank_ = 0;

    int localSize = getCommSize(localComm_);
    assert(masterRank_ < localSize);
    localRank_ = getCommRank(localComm_);
  }

  if (ENABLE_FT) {
    if (localComm_ != MPI_COMM_NULL) {
      createCommFT(&localCommFT_, localComm_);
    }
  }
}

/**
* create global communicator which contains only the manager and the master
* process of each process group
* the master processes of the process groups are the processes which have
* rank 0 in lcomm
* this communicator is used for communication between master processes of the
* process groups and the manager and the master processes to each other
*/
void MPISystem::initGlobalComm() {
  MPI_Group worldGroup;
  MPI_Comm_group(worldComm_, &worldGroup);

  std::vector<int> ranks(ngroup_ + 1);
  for (size_t i = 0; i < ngroup_; i++) {
    ranks[i] = int(i * nprocs_);
  }
  ranks.back() = managerRankWorld_;

  MPI_Group globalGroup;
  MPI_Group_incl(worldGroup, int(ranks.size()), &ranks[0], &globalGroup);

  MPI_Comm_create(worldComm_, globalGroup, &globalComm_);

  if (globalComm_ != MPI_COMM_NULL) {
    int globalSize = getCommSize( globalComm_);

    managerRank_ = globalSize - 1;
    std::cout << "new manager rank: " << managerRank_ << " \n";
    globalRank_ = getCommRank(globalComm_);
  }

  if (ENABLE_FT) {
    if (globalComm_ != MPI_COMM_NULL) {
      createCommFT(&globalCommFT_, globalComm_);
    }
  }
  /* mark master processes and manager process in Stats. this is necessary for postprocessing */
  Stats::setAttribute("group_manager", std::to_string(globalComm_ != MPI_COMM_NULL));
}

void MPISystem::initGlobalReduceCommm() {

  if (worldRank_ != managerRankWorld_) {
    int workerID = getCommRank(worldComm_);

    //      MPI_Comm globalReduceComm;
    int color = workerID % int(nprocs_);
    int key = workerID / int(nprocs_);
    MPI_Comm_split(worldComm_, color, key, &globalReduceComm_);

    if (ENABLE_FT) {
      createCommFT(&globalReduceCommFT_, globalReduceComm_);
    }
    int size = getCommSize(globalReduceComm_);
    //std::cout << "size if global reduce comm " << size << "\n";
    MPI_Barrier(globalReduceComm_);
  } else {
    MPI_Comm_split(worldComm_, MPI_UNDEFINED, -1, &globalReduceComm_);
  }
}

void MPISystem::createCommFT(simft::Sim_FT_MPI_Comm* commFT, CommunicatorType comm) {
  *commFT = new simft::Sim_FT_MPI_Comm_struct;

  (*commFT)->c_comm = comm;
  simft::Sim_FT_Initialize_new_comm(commFT, true);
}

std::vector<RankType> MPISystem::getFailedRanks(
    int numFailedProcs) {  // world comm still contains failed ranks at this point!
  std::vector<RankType> failedProcs(numFailedProcs);
  for (int i = 0; i < numFailedProcs; i++) {
    int failedRank;
    MPI_Recv(&failedRank, 1, MPI_INT, MPI_ANY_SOURCE, FT_FAILED_RANK_TAG,
             theMPISystem()->getWorldComm(), MPI_STATUS_IGNORE);
    failedProcs[i] = failedRank;
  }
  return failedProcs;
}

std::vector<RankType> MPISystem::getReusableRanks(int remainingProcs) {  // toDo matching send
  std::vector<RankType> reusableProcs(remainingProcs);
  for (int i = 0; i < remainingProcs; i++) {
    int reusableRank;
    MPI_Recv(&reusableRank, 1, MPI_INT, MPI_ANY_SOURCE, FT_REUSABLE_RANK_TAG,
             theMPISystem()->getWorldComm(), MPI_STATUS_IGNORE);
    reusableProcs[i] = reusableRank;
    std::cout << "reusable rank id " << reusableRank << "\n";
  }
  return reusableProcs;
}

void MPISystem::getReusableRanksSpare(
    std::vector<RankType>& reusableRanks) {  // update ranks after shrink
  for (size_t i = 0; i < reusableRanks.size(); i++) {
    std::cout << "getting new ranks \n";
    int reusableRank;
    MPI_Recv(&reusableRank, 1, MPI_INT, MPI_ANY_SOURCE, FT_REUSABLE_RANK_TAG,
             theMPISystem()->getSpareCommFT()->c_comm, MPI_STATUS_IGNORE);
    reusableRanks[i] = reusableRank;
    std::cout << "new rank id " << reusableRank << "\n";
  }
}

bool MPISystem::sendRankIds(std::vector<RankType>& failedRanks,
                            std::vector<RankType>& reusableRanks) {  // toDo matching send
  std::vector<MPI_Request> requests(failedRanks.size());
  for (size_t i = 0; i < failedRanks.size(); i++) {
    std::cout << "sending rank: " << failedRanks[i] << " to rank: " << reusableRanks[i] << "\n";
    MPI_Isend(&failedRanks[i], 1, MPI_INT, reusableRanks[i], FT_NEW_RANK_TAG,
              theMPISystem()->getSpareCommFT()->c_comm, &requests[i]);
  }
  bool success = true;  // so far no time out failing possible
  
  for (size_t i = 0; i < failedRanks.size(); i++) {
    MPI_Wait(&requests[i],
             MPI_STATUS_IGNORE);  // toDo implement timeout and remove time-outed ranks immedeately
  }
  auto lastIndex = failedRanks.size();  // needs to be increased in case of ranks with time-out
  bool recoveryFailed;
  if (success) {  // send recovery status
    recoveryFailed = false;
    std::vector<MPI_Request> requests2(failedRanks.size());
    for (size_t i = 0; i < failedRanks.size(); i++) {
      MPI_Isend(&recoveryFailed, 1, MPI::BOOL, reusableRanks[i], FT_RECOVERY_STATUS_TAG,
                theMPISystem()->getSpareCommFT()->c_comm, &requests2[i]);
    }
    for (size_t i = 0; i < failedRanks.size(); i++) {
      MPI_Wait(&requests2[i], MPI_STATUS_IGNORE);
    }
    reusableRanks.erase(reusableRanks.begin(), reusableRanks.begin() + lastIndex);
  }
  return recoveryFailed;  // so far failing not implemented
}

void MPISystem::sendRecoveryStatus(bool failedRecovery,
                                   std::vector<RankType>& newReusableRanks) {  // toDo matching send
  std::vector<MPI_Request> requests(newReusableRanks.size());
  for (size_t i = 0; i < newReusableRanks.size(); i++) {
    MPI_Isend(&failedRecovery, 1, MPI::BOOL, newReusableRanks[i], FT_RECOVERY_STATUS_TAG,
              theMPISystem()->getWorldComm(), &requests[i]);
  }
  for (size_t i = 0; i < newReusableRanks.size(); i++) {
    MPI_Wait(&requests[i], MPI_STATUS_IGNORE);  // toDo implement timeout
  }
}

void MPISystem::sendExcludeSignal(std::vector<RankType>& reusableRanks) {
  std::vector<MPI_Request> requests(reusableRanks.size());

  for (size_t i = 0; i < reusableRanks.size(); i++) {
    int excludeData;  // not used so far
    MPI_Isend(&excludeData, 0, MPI_INT, reusableRanks[i], FT_EXCLUDE_TAG,
              theMPISystem()->getSpareCommFT()->c_comm, &requests[i]);
  }

  for (size_t i = 0; i < reusableRanks.size(); i++) {
    MPI_Wait(&requests[i], MPI_STATUS_IGNORE);  // toDo implement timeout
  }
}

void MPISystem::sendShrinkSignal(std::vector<RankType>& reusableRanks) {
  std::vector<MPI_Request> requests(reusableRanks.size());

  for (size_t i = 0; i < reusableRanks.size(); i++) {
    int shrinkData;  // not used so far
    MPI_Isend(&shrinkData, 0, MPI_INT, reusableRanks[i], FT_SHRINK_TAG,
              theMPISystem()->getSpareCommFT()->c_comm, &requests[i]);
    std::cout << "sending shrink signal to " << reusableRanks[i] << "\n";
  }

  for (size_t i = 0; i < reusableRanks.size(); i++) {
    MPI_Wait(&requests[i], MPI_STATUS_IGNORE);  // toDo implement timeout
  }
}

void MPISystem::sendReusableSignal() {
  std::cout << "sending reusable signal \n";
  MPI_Send(&worldRank_, 1, MPI_INT, managerRankWorld_, FT_REUSABLE_RANK_TAG,
           theMPISystem()->getWorldComm());
}

void MPISystem::sendReusableSignalSpare() {
  MPI_Send(&worldRank_, 1, MPI_INT, managerRankFT_, FT_REUSABLE_RANK_TAG,
           theMPISystem()->getSpareCommFT()->c_comm);
}

void MPISystem::sendFailedSignal() {  // worldComm still includes failed ranks at this point!
  MPI_Request sendReq;
  MPI_Isend(&worldRank_, 1, MPI_INT, managerRankWorld_, FT_FAILED_RANK_TAG,
            theMPISystem()->getWorldComm(), &sendReq);
}

bool MPISystem::receiveRecoverStatus() {
  bool recoveryState;

  // receive recovery state
  MPI_Recv(&recoveryState, 1, MPI::BOOL, managerRankWorld_, FT_RECOVERY_STATUS_TAG,
           theMPISystem()->getWorldComm(), MPI_STATUS_IGNORE);
  return recoveryState;
}
/**
 * Routine that handles waiting for idle procs that can be reused.
 * Processors stay in this routine until the can be reused + they contribute in all shrink and split
 * operations.
 */
void MPISystem::waitForReuse() {
  int excludeFlag = 0;
  int newRankFlag = 0;
  int shrinkFlag = 0;
  MPI_Status shrinkStatus;
  MPI_Status excludeStatus;
  MPI_Status newRankStatus;
  while (true) {
    MPI_Iprobe(managerRankFT_, FT_NEW_RANK_TAG, theMPISystem()->getSpareCommFT()->c_comm,
               &newRankFlag, &newRankStatus);
    if (newRankFlag) {
      int rankID;
      // receive rank ID
      MPI_Recv(&rankID, 1, MPI_INT, managerRankFT_, FT_NEW_RANK_TAG,
               theMPISystem()->getSpareCommFT()->c_comm, MPI_STATUS_IGNORE);
      std::cout << "spare rank " << worldRank_ << " receives new rank: " << rankID << "\n";

      newRankFlag = 0;
      // receive recovery state to check if recovery really succeeded
      bool recoveryFailed;
      MPI_Recv(&recoveryFailed, 1, MPI::BOOL, managerRankFT_, FT_RECOVERY_STATUS_TAG,
               theMPISystem()->getSpareCommFT()->c_comm, MPI_STATUS_IGNORE);
      if (!recoveryFailed) {
        worldRank_ = rankID;
        return;
      }
      // otherwise wait for recovery again
    }
    MPI_Iprobe(managerRankFT_, FT_EXCLUDE_TAG, theMPISystem()->getSpareCommFT()->c_comm,
               &excludeFlag, &excludeStatus);
    if (excludeFlag) {
      int excludeData;
      excludeFlag = 0;
      // receive exclude data ; so far nothing there
      std::cout << "spare rank " << worldRank_ << " gets excluded \n";

      MPI_Recv(&excludeData, 0, MPI_INT, managerRankFT_, FT_EXCLUDE_TAG,
               theMPISystem()->getSpareCommFT()->c_comm, MPI_STATUS_IGNORE);
      // perform split with color undefined
      int color = MPI_UNDEFINED;
      int key = getCommSize(spareCommFT_->c_comm);
      MPI_Comm_split(spareCommFT_->c_comm, color, key, &worldComm_);
    }
    MPI_Iprobe(managerRankFT_, FT_SHRINK_TAG, theMPISystem()->getSpareCommFT()->c_comm, &shrinkFlag,
               &shrinkStatus);
    if (shrinkFlag) {
      shrinkFlag = 0;
      int shrinkData;  // not used so far
      // receive exclude data ; so far nothing there
      std::cout << "spare rank " << worldRank_ << " performs shrink \n";
      MPI_Recv(&shrinkData, 0, MPI_INT, managerRankFT_, FT_SHRINK_TAG,
               theMPISystem()->getSpareCommFT()->c_comm, MPI_STATUS_IGNORE);
      // perform shrink
      simft::Sim_FT_MPI_Comm newSpareCommFT;
      MPI_Comm_shrink(
          theMPISystem()->getSpareCommFT(),
          &newSpareCommFT);  // remove dead processors from spareComm(worldComm + reusable ranks)
      deleteCommFTAndCcomm(&spareCommFT_);
      createCommFT(&spareCommFT_, newSpareCommFT->c_comm);  // removes dead processors from
                                                            // worldComm
      // delete temporary communicator
      deleteCommFT(&newSpareCommFT);
      // adjust manger rank in spareComm as it has changed during shrink
      int ftCommSize = getCommSize(spareCommFT_->c_comm);
      managerRankFT_ = ftCommSize - 1;
      worldRank_ = getCommRank(spareCommFT_->c_comm);
      sendReusableSignalSpare();
    }
  }
}

void MPISystem::deleteCommFT(simft::Sim_FT_MPI_Comm* commFT) {  // does not delete c_comm ->
                                                                // important if c_comm is used in
                                                                // other FT comm (e.g. through
                                                                // creatCommFT)
  if (commFT != NULL && *commFT != NULL) {  // delete old communicator of exists
    if ((*commFT)->c_comm != MPI_COMM_WORLD && (*commFT)->c_comm != MPI_COMM_NULL) {
      simft::Sim_FT_MPI_Comm_free2(commFT);  // do not delete c_comm
    }
  }
}
void MPISystem::deleteCommFTAndCcomm(simft::Sim_FT_MPI_Comm* commFT, CommunicatorType* ccommCopy) {
  deleteCommFTAndCcomm(commFT);
  if (*ccommCopy != MPI_COMM_NULL) {
    *ccommCopy = MPI_COMM_NULL;
  }
}

void MPISystem::deleteCommFTAndCcomm(simft::Sim_FT_MPI_Comm* commFT) {
  if (commFT != NULL && *commFT != NULL) {  // delete old communicator of exists
    if ((*commFT)->c_comm != MPI_COMM_WORLD && (*commFT)->c_comm != MPI_COMM_NULL) {
      simft::Sim_FT_MPI_Comm_free(commFT);  // do not delete c_comm
    }
  }
}

bool MPISystem::recoverCommunicators(bool groupAlive,
                                     std::vector<std::shared_ptr<ProcessGroupManager>>
                                         failedGroups) {  // toDo fix multiple failed groups
  assert(ENABLE_FT && "this funtion is only availabe if FT enabled!");
  // std::cout << "start recovery \n";
  // revoke commmworld
  // theStatsContainer()->setTimerStart("recoverComm-revoke");
  // WORLD_MANAGER_EXCLUSIVE_SECTION{
  // MPI_Comm_revoke( theMPISystem()->getWorldCommFT() );
  //}
  // theStatsContainer()->setTimerStop("recoverComm-revoke");
  // std::cout << "revoked MPI comm \n";
  // shrink world
  // theStatsContainer()->setTimerStart("recoverComm-shrink");
  simft::Sim_FT_MPI_Comm newSpareCommFT;
  simft::Sim_FT_MPI_Comm newWorldCommFT;
  WORLD_MANAGER_EXCLUSIVE_SECTION {
    // indicate shrink to reusable ranks
    sendShrinkSignal(reusableRanks_);
  }
  // shrink of all processors including reusable ones
  MPI_Comm_shrink(
      theMPISystem()->getSpareCommFT(),
      &newSpareCommFT);  // remove dead processors from spareComm(worldComm + reusable ranks)
  deleteCommFTAndCcomm(&spareCommFT_);

  createCommFT(&spareCommFT_, newSpareCommFT->c_comm);  // removes dead processors from worldComm
  deleteCommFT(&newSpareCommFT);

  // adjust manger rank in spareComm as it has changed durin shrink
  int ftCommSize = getCommSize( spareCommFT_->c_comm);
  managerRankFT_ = ftCommSize - 1;
  std::vector<RankType> newReusableRanks;

  MPI_Comm_shrink(theMPISystem()->getWorldCommFT(),
                  &newWorldCommFT);  // remove dead processors from current worldComm

  bool failedRecovery = true;
  int sizeNew = getCommSize( newWorldCommFT->c_comm);
  // deleteing tompary world comm
  deleteCommFTAndCcomm(&newWorldCommFT);
  WORLD_MANAGER_EXCLUSIVE_SECTION {  // get failed ranks
    int sizeOld = getWorldSize();
    int sizeSpare = getCommSize( spareCommFT_->c_comm);
    std::cout << "size old = " << sizeOld << "\n";
    std::cout << "size new = " << sizeNew << "\n";
    std::cout << "size spare = " << sizeSpare << "\n";

    int numFailedRanks = sizeOld - sizeNew;
    std::vector<RankType> failedRanks =
        getFailedRanks(numFailedRanks);  // has to be solved differently with ULFM
    std::cout << "nprocs - numFailed " << failedGroups.size() * nprocs_ - numFailedRanks << "\n";
    newReusableRanks = getReusableRanks(static_cast<int>(failedGroups.size() * nprocs_ - numFailedRanks));
    getReusableRanksSpare(reusableRanks_);  // update ranks of reusable ranks
    // toDO reusableRanks might be outdated due to new failures there
    bool enoughSpareProcs = sizeSpare - sizeNew >= numFailedRanks;
    std::cout << "enoughSpareProcs: " << enoughSpareProcs << "\n";
    // std::cout << "failedRanks: " << failedRanks[0] << failedRanks.size() << " ;reusable Ranks: "
    // << newReusableRanks[0] << newReusableRanks.size() << "\n";
    bool failedSendingRankIds = false;
    if (enoughSpareProcs) {  // send failed ranks to reusable procs so they can set their worldRank
                             // accordingly
      failedSendingRankIds =
          sendRankIds(failedRanks,
                      reusableRanks_);  // check with timeout if reusable ranks are still available;
      std::cout << "finished assigning ranks \n";
      // otherwise fail recovery of broken processgroup
      // remove failed spare procs from reusableRanks_
      if (!failedSendingRankIds) {
        failedRecovery = false;
      }
    }
    if (!enoughSpareProcs or failedSendingRankIds) {
      failedRecovery = true;
    }
    std::cout << "failed recovery: " << failedRecovery << "\n";
    std::cout << "sending recovery status \n";
    sendRecoveryStatus(failedRecovery, newReusableRanks);
    // order shrink with color 0 to reusable ranks as they need to be excluded for this recovery
    // process
    std::cout << "sending exclude signal \n";
    sendExcludeSignal(reusableRanks_);
    if (failedRecovery) {  // add in case of failure newreusable ranks to vector
      getReusableRanksSpare(newReusableRanks);  // update ranks to value in spareft communicator

      reusableRanks_.insert(reusableRanks_.end(), newReusableRanks.begin(), newReusableRanks.end());
    }
  }

  int color;
  if (!groupAlive) {  // check if group was recovered and mark as reusable
    sendReusableSignal();
    bool recoveryFailed = receiveRecoverStatus();  // check if exclude signal
    std::cout << "recovery failed: " << recoveryFailed << "\n";
    deleteCommFTAndCcomm(&localCommFT_, &localComm_);

    if (recoveryFailed) {
      color = MPI_UNDEFINED;
      int key = worldRank_;
      // update worldRank_ to spare communicator
      worldRank_ = getCommRank(spareCommFT_->c_comm);
      sendReusableSignalSpare();  // send rank in ft communicator to master
      deleteCommFTAndCcomm(&worldCommFT_, &worldComm_);
      MPI_Comm_split(spareCommFT_->c_comm, color, key, &worldComm_);
      // delete communicators
      deleteCommFTAndCcomm(&globalCommFT_, &globalComm_);
      deleteCommFTAndCcomm(&globalReduceCommFT_, &globalReduceComm_);
      // sendReusableSignal(theMPISystem()->getSpareCommFT()->c_comm);
      waitForReuse();  // does only leave this routine when reused later //participates in future
                       // shrinks
      color = 1;
    } else {
      deleteCommFTAndCcomm(&worldCommFT_, &worldComm_);
    }
  }
  // theStatsContainer()->setTimerStop("recoverComm-shrink");
  // std::cout << "shrinked communicator \n";
  // split off alive procs. this will be the new WorldComm
  // processes of dead groups set color to MPI_UNDEFINED. in this case
  // MPI_Comm_split returns MPI_COMMM_NULL
  color = 1;  // all processors at this point belong to alive or recovered process groups
  int key = worldRank_;
  if (groupAlive) {
    deleteCommFTAndCcomm(&worldCommFT_, &worldComm_);
  }
  std::cout << "performing MPI split \n";

  MPI_Comm_split(spareCommFT_->c_comm, color, key, &worldComm_);
  
  int worldSize = getWorldSize();
  assert((worldSize - 1) % nprocs_ == 0);
  ngroup_ = (worldSize - 1) / nprocs_;
  MPI_Comm_rank(worldComm_, &worldRank_); 
  managerRankWorld_ = worldSize - 1;

  if (worldComm_ != MPI_COMM_NULL) {
    createCommFT(&worldCommFT_, worldComm_);
  }
  if (!groupAlive) {
    // toDo:init local comm?
    std::cout << "initialize local comm \n";
    initLocalComm();
    std::cout << "initialized local comm \n";
  } else {
    CommunicatorType tmp;
    color = MPI_UNDEFINED;
    key = -1;
    MPI_Comm_split(worldComm_, color, key, &tmp);
  }
  std::cout << "initializing global comm \n";
  deleteCommFTAndCcomm(&globalCommFT_, &globalComm_);
  initGlobalComm();
  std::cout << "initializing global reduce comm \n";
  deleteCommFTAndCcomm(&globalReduceCommFT_, &globalReduceComm_);
  initGlobalReduceCommm();

  // toDo return fixed process group IDs
  std::cout << "returning \n";
  //  MPI_Barrier( worldComm_ );

  return failedRecovery;
}

}  // namespace combigrid
