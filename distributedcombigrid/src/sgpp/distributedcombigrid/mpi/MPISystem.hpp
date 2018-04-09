#ifndef MPISYSTEM_HPP
#define MPISYSTEM_HPP

#include <mpi.h>
#include <ostream>
#include <vector>
#include <assert.h>

#include "sgpp/distributedcombigrid/mpi/MPISystemID.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/mpi_fault_simulator/MPI-FT.h"

#define MASTER_EXCLUSIVE_SECTION if( combigrid::theMPISystem()->isMaster() )

#define GLOBAL_MANAGER_EXCLUSIVE_SECTION if( combigrid::theMPISystem()->isGlobalManager() )

#define WORLD_MANAGER_EXCLUSIVE_SECTION  if( combigrid::theMPISystem()->isWorldManager() )

namespace combigrid {

/** MPI Communication System
 *
 * This class encapsulates the access to the different communicators being used.
 *
 * Each process belongs to different MPI communicators which can overlap. If a
 * process does not belong to a specific communicator the corresponding get function
 * returns MPI_COMM_NULL.
 *
 * CommWorld: contains all processes. usually equal to MPI_COMM_WORLD. if fault
 * tolerance is enabled the processes of groups detected as failed will be removed
 *
 * GlobalComm: contains the manager and the master process of each process group
 *
 * SpareCommFT: contains inactive but reusable ranks. Only used for FT simulations.
 *
 * LocalComm: contains the processes of the group a process is assigned to. this
 * is MPI_COMM_NULL for the manager process. this is the communicator the
 * application MUST use for the computations. the application MUST NEVER use
 * MPI_COMM_WORLD internally.
 *
 * GlobalReduceComm: contains the processes for the global reduction. in the
 * case of equally sized process groups (at the moment this is the only option)
 * this communicator contains all processes which have the same rank in
 * LocalComm. it is MPI_COMM_NULL on the manager.
 *
 * getXXXCommFT returns the fault tolerant equivalent of Communicator XXX
 *
 * getXXXRank return the rank of the process in the Communicator XXX. If the
 * process is not part of Communicator XXX it returns MPI_UNDEFINED.
 *
 * * getManagerRankWorld returns the rank of the manager process in the
 * WorldComm.
 *
 * getManagerRank returns the rank of the manager process in the Communicator
 * GlobalComm. if the process is not part of GlobalComm it returns MPI_UNDEFINED
 *
 * getMasterRank returns the rank of the master process in the Communicator
 * LocalComm. if the process is not part of LocalComm it returns MPI_UNDEFINED
 *
 * isManager returns true if the process is the manager process
 *
 * isLocalMaster returns true if the calling process is the master process of
 * its process group
 */
class ProcessGroupManager;
class MPISystem {
 public:
  ~MPISystem();

  MPISystem( MPISystem const & ) = delete;

  MPISystem& operator=( MPISystem const & ) = delete;

  /**initializes MPI system including local, global and allreduce communicator
   * for specified number of groups and number of processors per group
   */
  void init( size_t ngroups, size_t nprocs );

  /**
   *initializes MPI system including global and allreduce communicator
   *local communicator is given here and not set by the init procedure
   */
  void init( size_t ngroups, size_t nprocs, CommunicatorType lcomm );

  /**
   * returns the world communicator which contains all ranks (excluding spare ranks)
   */
  inline const CommunicatorType& getWorldComm() const;

  /**
   * returns the world communicator which contains all manager and master ranks
   */
  inline const CommunicatorType& getGlobalComm() const;

  /**
   * returns the local communicator which contains all ranks within the process group of caller
   */
  inline const CommunicatorType& getLocalComm() const;

  /**
   * returns the global reduce communicator which contains all ranks with wich the rank needs to
   * communicate in global allreduce step
   * All of these ranks are responsible for the same area in the domain
   */
  inline const CommunicatorType& getGlobalReduceComm() const;

  /**
   * returns the fault tolerant version of the world comm (excluding spare ranks)
   */
  inline simft::Sim_FT_MPI_Comm getWorldCommFT();

  /**
   * returns the communicator also containing spare processors (or ranks)
   */
  inline simft::Sim_FT_MPI_Comm getSpareCommFT();

  /**
   * returns the fault tolerant version of the global comm
   */
  inline simft::Sim_FT_MPI_Comm getGlobalCommFT();

  /**
   * returns the fault tolerant version of the local comm
   */
  inline simft::Sim_FT_MPI_Comm getLocalCommFT();

  /**
   * returns the fault tolerant version of the allreduce comm
   */
  inline simft::Sim_FT_MPI_Comm getGlobalReduceCommFT();

  /**
   * returns MPI rank number in world comm
   */
  inline const RankType& getWorldRank() const;

  /**
   * returns MPI rank number in global comm
   */
  inline const RankType& getGlobalRank() const;

  /**
   * returns MPI rank number in local comm
   */
  inline const RankType& getLocalRank() const;

  /**
   * returns MPI rank number of manager in world comm
   */
  inline const RankType& getManagerRankWorld() const;

  /**
   * returns MPI rank number of manager in global comm
   */
  inline const RankType& getManagerRank() const;

  /**
   * returns MPI rank number of master in local comm
   */
  inline const RankType& getMasterRank() const;

  /**
   * returns boolean that indicates if caller is manager in world comm
   */
  inline bool isWorldManager() const;

  /**
   * returns boolean that indicates if caller is manager in global comm
   */
  inline bool isGlobalManager() const;

  /**
   * returns boolean that indicates if caller is master in local comm
   */
  inline bool isMaster() const;

  /**
   * returns the number of process groups
   */
  inline size_t getNumGroups() const;

  /**
   * returns the number of processors per process group
   */
  inline size_t getNumProcs() const;

  /**
   * returns boolean that indicates if MPISystem is initialized
   */
  inline bool isInitialized() const;

  /**
   * This routine starts the recovery procesdure.
   * groupAlive indicates if the process group of the calling rank is alive
   * failedGroups is a vector of the failed process groups
   */
  bool recoverCommunicators( bool groupAlive, std::vector< std::shared_ptr< ProcessGroupManager >> failedGroups = std::vector< std::shared_ptr< ProcessGroupManager >>(0) );

  /**
   * This routine frees the specified fault tolerant MPI communicator
   * The corresponding non-fault tolerant communicator associated with the FT-communicator
   * is not destroyed!
   */
  void deleteCommFT(simft::Sim_FT_MPI_Comm *comm);

  /**
   * This routine frees the specified fault tolerant MPI communicator
   * The corresponding non-fault tolerant communicator associated with the FT-communicator
   * is also destroyed!
   */
  void deleteCommFTAndCcomm(simft::Sim_FT_MPI_Comm *commFT, CommunicatorType *ccommCopy);

  /**
   * This routine frees the specified fault tolerant MPI communicator
   * The corresponding non-fault tolerant communicator associated with the FT-communicator
   * is also destroyed!
   */
  void deleteCommFTAndCcomm(simft::Sim_FT_MPI_Comm *comm);

  /**
   * sends a message to the manager that this rank has failed -> used for FT simulator
   */
  void sendFailedSignal();

 private:
  explicit MPISystem();

  friend MPISystemID theMPISystem();
  /**
   * checks if initialized
   */
  inline void checkPreconditions() const;
  /**
   * checks if initialized + if FT enabled
   */
  inline void checkPreconditionsFT() const;

  /* create the global communicators for the global reduce.
   * all processes which have local rank i in their process group will be grouped in a
   * distinct communicator.
   * this will only work if all pgroups have same size and dsg uses the same assignment
   * of subspaces to processes in each pgroup
   */
  void initGlobalReduceCommm();

  /**
   * creates a FT communicator associated with comm
   */
  void createCommFT( simft::Sim_FT_MPI_Comm* commFT, CommunicatorType comm );

  /**
   * initializes local comm + FT version if FT_ENABLED
   */
  void initLocalComm();

  /**
   * initializes global comm + FT version if FT_ENABLED
   */
  void initGlobalComm();

  /**
   * Send signal to manager that this rank is reusable -> becomes spare process
   * The message is send in world comm
   */
  void sendReusableSignal();

  /**
   * Send signal to manager that this rank is reusable -> becomes spare processor
   * The message is send in spare comm
   */
  void sendReusableSignalSpare();

  /**
   * a process enters this routine if it becomes a spare process and waits for reuse
   * all shrinking operations of spare communicator are performed here too
   */
  void waitForReuse();

  /**
   * receive the status if failed process group could be recoverd
   * (only relevant for procs in a failed group that have not failed themselves)
   */
  bool receiveRecoverStatus();

  /**
   * start shrinking procedure at all ranks waiting for reuse
   */
  void sendShrinkSignal(std::vector<RankType>& reusableRanks);

  /**
   * When splitting the world comm after fault all spare processors need to be excluded from world comm.
   * This routine sends a signal to all waiting spare processes so that they start the MPI split and exclude themselves.
   */
  void sendExcludeSignal(std::vector<RankType>& reusableRanks);

  /**
   * sends to all processes that can be reused if the group recovery failed
   */
  void sendRecoveryStatus(bool failedRecovery, std::vector<RankType>& newReusableRanks );

  /**
   * In case of successful recovery send the vector of failed ranks to the reusable ranks.
   * The reusable ranks will copy the corresponding failed ID and will then replace this process.
   */
  bool sendRankIds(std::vector<RankType>& failedRanks, std::vector<RankType>& reusableRanks );

  /**
   * returns the vector of the rank numbers of the newly added spare procs due to process fault
   * which could not be recovered
   */
  std::vector<RankType> getReusableRanks( int remainingProcs );

  /**
   * returns the vector of the rank numbers in the spare communicator
   * of the newly added spare procs due to process fault which could not be recovered
   */
  void getReusableRanksSpare(std::vector<RankType>& reusableRanks);

  /**
   * returns the vector of the failed ranks that caused the current recovery procedure
   */
  std::vector<RankType> getFailedRanks( int numFailedProcs );


  bool initialized_; //is MPISystem initialized?

  CommunicatorType worldComm_; //contains all processes that are active

  //contains alive procs from dead process groups and manager
  simft::Sim_FT_MPI_Comm spareCommFT_;

  CommunicatorType globalComm_; //contains the manager and master processes

  CommunicatorType localComm_; //contains all processes in process group

  /**
   * contains all processes that share same domain in other process groups
   * -> only communicate with these ranks during allreduce
   */
  CommunicatorType globalReduceComm_;

  simft::Sim_FT_MPI_Comm worldCommFT_; //FT version of world comm

  simft::Sim_FT_MPI_Comm globalCommFT_; //FT version of global comm

  simft::Sim_FT_MPI_Comm localCommFT_; //FT version of local comm

  simft::Sim_FT_MPI_Comm globalReduceCommFT_; //FT version of global reduce comm

  RankType worldRank_; //rank number in world comm

  RankType globalRank_; //rank number in global comm

  RankType localRank_; //rank number in local comm

  RankType globalReduceRank_; //rank number in global reduce comm

  RankType managerRank_; //rank number of manager in global comm

  RankType managerRankWorld_; //rank number of manager in world comm

  RankType managerRankFT_; //rank number of manager in spare comm

  RankType masterRank_; //rank number of master in local comm

  size_t ngroup_; //number of process groups

  size_t nprocs_; //number of processes per process group

  //ranks that er still functional but not assigned to any process group
  std::vector<RankType> reusableRanks_;
};

/*!\name MPI communication system setup functions */
//@{
inline MPISystemID theMPISystem();
//@}

/** Returns a handle to the MPI communication system.
 *
 * This function returns a handle to the MPI communication system. This handle
 * can be used to configure the communication system or to acquire
 * the current settings. The function expects that MPI has already been properly
 * initialized (e.g. via the MPI_Init() or any similar function). In case MPI
 * was not initialized, a \a std::runtime_error exception is thrown.
 *
 * \return Handle to the MPI communication system.
 * \exception std::runtime_error MPI system is not initialized.
 */
inline MPISystemID theMPISystem() {
  static MPISystemID system(new MPISystem());
  return system;
}


inline void MPISystem::checkPreconditions() const{
  assert( initialized_ && "MPI System not initialized!");
}


inline void MPISystem::checkPreconditionsFT() const{
  checkPreconditions();
  assert( ENABLE_FT && "Fault Tolerance not enabled!" );
}


inline const CommunicatorType& MPISystem::getWorldComm() const {
  checkPreconditions();

  return worldComm_;
}


inline const CommunicatorType& MPISystem::getGlobalComm() const{
  checkPreconditions();

  return globalComm_;
}


inline const CommunicatorType& MPISystem::getLocalComm() const{
  checkPreconditions();

  return localComm_;
}


inline const CommunicatorType& MPISystem::getGlobalReduceComm() const{
  checkPreconditions();

  return globalReduceComm_;
}

inline simft::Sim_FT_MPI_Comm MPISystem::getWorldCommFT(){
  checkPreconditionsFT();

  return worldCommFT_;
}

inline simft::Sim_FT_MPI_Comm MPISystem::getSpareCommFT(){
  checkPreconditionsFT();

  return spareCommFT_;
}

inline simft::Sim_FT_MPI_Comm MPISystem::getGlobalCommFT(){
  checkPreconditionsFT();

  return globalCommFT_;
}


inline simft::Sim_FT_MPI_Comm MPISystem::getLocalCommFT(){
  checkPreconditionsFT();

  return localCommFT_;
}


inline simft::Sim_FT_MPI_Comm MPISystem::getGlobalReduceCommFT(){
  checkPreconditionsFT();

  return globalReduceCommFT_;
}


inline const RankType& MPISystem::getWorldRank() const{
  checkPreconditions();

  return worldRank_;
}


inline const RankType& MPISystem::getGlobalRank() const{
  checkPreconditions();

  return globalRank_;
}


inline const RankType& MPISystem::getLocalRank() const{
  checkPreconditions();

  return localRank_;
}


inline const RankType& MPISystem::getManagerRankWorld() const{
  checkPreconditions();

  return managerRankWorld_;
}


inline const RankType& MPISystem::getManagerRank() const{
  checkPreconditions();

  return managerRank_;
}


inline const RankType& MPISystem::getMasterRank() const{
  checkPreconditions();

  return masterRank_;
}


inline bool MPISystem::isWorldManager() const{
  return ( worldRank_ == managerRankWorld_ );
}


inline bool MPISystem::isGlobalManager() const{
  return ( globalRank_ == managerRank_ );
}


inline bool MPISystem::isMaster() const{
  return ( localRank_ == masterRank_ );
}


inline size_t MPISystem::getNumGroups() const{
  checkPreconditions();

  return ngroup_;
}

inline size_t MPISystem::getNumProcs() const{
  checkPreconditions();

  return nprocs_;
}


inline bool MPISystem::isInitialized() const{
  return initialized_;
}

/*
// operators
std::ostream& operator<<(std::ostream& os, const MPISystem& ms);
std::ostream& operator<<(std::ostream& os, const MPISystemID& ms);
std::ostream& operator<<(std::ostream& os, const ConstMPISystemID& ms);
*/

} //namespace combigrid




#endif // MPISYSTEM_HPP
