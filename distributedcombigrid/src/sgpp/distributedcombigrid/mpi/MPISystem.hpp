#ifndef MPISYSTEM_HPP
#define MPISYSTEM_HPP

#include <mpi.h>
#include <ostream>
#include <vector>
#include <type_traits>
#include <assert.h>

#include "sgpp/distributedcombigrid/mpi/MPISystemID.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/mpi_fault_simulator/MPI-FT.h"

#define TEAM_LEADER_EXCLUSIVE_SECTION if( combigrid::theMPISystem()->isTeamLeader() )

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

  void init( size_t ngroups, size_t nprocs );

  void init( size_t ngroups, std::vector<size_t> nprocsByGroup );

  void init( size_t ngroups, size_t nprocs, CommunicatorType lcomm );

  inline const CommunicatorType& getWorldComm() const;

  inline const CommunicatorType& getGlobalComm() const;

  inline const CommunicatorType& getLocalComm() const;

  inline const CommunicatorType& getTeamComm() const;

  inline const CommunicatorType& getGlobalReduceComm() const;

  inline simft::Sim_FT_MPI_Comm getWorldCommFT();

  inline simft::Sim_FT_MPI_Comm getSpareCommFT();

  inline simft::Sim_FT_MPI_Comm getGlobalCommFT();

  inline simft::Sim_FT_MPI_Comm getLocalCommFT();

  inline simft::Sim_FT_MPI_Comm getTeamCommFT();

  inline simft::Sim_FT_MPI_Comm getGlobalReduceCommFT();

  inline void deleteReduceCommsFTAndComm();

  inline const RankType& getWorldRank() const;

  inline const RankType& getGlobalRank() const;

  inline const RankType& getLocalRank() const;

  inline const RankType& getManagerRankWorld() const;

  inline const RankType& getManagerRank() const;

  inline const RankType& getMasterRank() const;

  inline bool isWorldManager() const;

  inline bool isGlobalManager() const;

  inline bool isMaster() const;

  inline bool isTeamLeader() const;

  inline const RankType& getTeamColor() const;

  inline size_t getNumGroups() const;

  inline size_t getNumProcs( GroupType group ) const;

  inline RankType getGroupBaseWorldRank( GroupType group ) const;

  inline bool isInitialized() const;

  bool recoverCommunicators( bool groupAlive, std::vector< std::shared_ptr< ProcessGroupManager >> failedGroups = std::vector< std::shared_ptr< ProcessGroupManager >>(0) );
  void deleteCommFT(simft::Sim_FT_MPI_Comm *comm);

  void deleteCommFTAndCcomm(simft::Sim_FT_MPI_Comm *comm, CommunicatorType *ccomm);

  void sendFailedSignal();

 private:
  enum class InitializationStage {
    // Nothing is initialized yet
    PRE_INIT,
    // The world communicator is initialized
    WORLD_INIT,
    // The local communicator is initialized
    LOCAL_INIT,
    // The team communicator is initilalized
    TEAM_INIT,
    // The global communicator is initialized
    GLOBAL_INIT,
    // The reduce communicator is initialized
    ALL_INIT
  };

  explicit MPISystem();

  friend MPISystemID theMPISystem();

  inline void checkPreconditions(InitializationStage requiredStage = InitializationStage::ALL_INIT) const;

  inline void checkPreconditionsFT(InitializationStage requiredStage = InitializationStage::ALL_INIT) const;

  void initWorldComm( std::vector<size_t> nprocsByGroup );

  /* create the global communicators for the global reduce.
   * all processes which have local rank i in their process group will be grouped in a
   * distinct communicator.
   * this will only work if all pgroups have same size and dsg uses the same assignment
   * of subspaces to processes in each pgroup
   */
  void initGlobalReduceCommm();

  void createCommFT( simft::Sim_FT_MPI_Comm* commFT, CommunicatorType comm );

  void initLocalComm();

  void initLocalComm(CommunicatorType lcomm);

  void initGlobalComm();

  void initTeamComm();

  void sendReusableSignal();

  void sendReusableSignalSpare();

  void waitForReuse();

  bool receiveRecoverStatus();

  void sendShrinkSignal(std::vector<RankType>& reusableRanks);

  void sendExcludeSignal(std::vector<RankType>& reusableRanks);

  void sendRecoveryStatus(bool failedRecovery, std::vector<RankType>& newReusableRanks );

  bool sendRankIds(std::vector<RankType>& failedRanks, std::vector<RankType>& reusableRanks );

  std::vector<RankType> getReusableRanks( int remainingProcs );

  void getReusableRanksSpare(std::vector<RankType>& reusableRanks);

  std::vector<RankType> getFailedRanks( int numFailedProcs );



  InitializationStage initialized_;

  CommunicatorType worldComm_;

  CommunicatorType globalComm_;

  CommunicatorType localComm_;

  CommunicatorType teamComm_;

  CommunicatorType globalReduceComm_;

  simft::Sim_FT_MPI_Comm worldCommFT_;

  simft::Sim_FT_MPI_Comm globalCommFT_;

  //contains alive procs from dead process groups and manager
  simft::Sim_FT_MPI_Comm spareCommFT_;

  simft::Sim_FT_MPI_Comm localCommFT_;

  simft::Sim_FT_MPI_Comm teamCommFT_;

  simft::Sim_FT_MPI_Comm globalReduceCommFT_;

  RankType worldRank_;

  RankType globalRank_;

  RankType localRank_;

  RankType teamRank_;

  RankType teamColor_;

  RankType globalReduceRank_;

  RankType managerRank_;

  RankType managerRankWorld_;

  RankType managerRankFT_;

  RankType masterRank_;

  RankType teamLeaderRank_;

    // Each group inhabits the ranks range
    // [groupBaseRank[g], groupBaseRank[g] + nprocsByGroup_[g]) in worldComm
  std::vector<size_t> nprocsByGroup_;
  std::vector<RankType> groupBaseWorldRank_;
  GroupType group_;

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


inline void MPISystem::checkPreconditions(InitializationStage requiredStage) const {
  using OrdinalType = typename std::underlying_type<InitializationStage>::type;
  auto requiredStageNr = static_cast<OrdinalType>(requiredStage);
  auto initializedStageNr = static_cast<OrdinalType>(initialized_);
  assert( requiredStageNr <= initializedStageNr && "MPI System not initialized!");
}


inline void MPISystem::checkPreconditionsFT(InitializationStage requiredStage) const{
  checkPreconditions( requiredStage );
  assert( ENABLE_FT && "Fault Tolerance not enabled!" );
}


inline const CommunicatorType& MPISystem::getWorldComm() const {
  checkPreconditions( InitializationStage::WORLD_INIT );

  return worldComm_;
}


inline const CommunicatorType& MPISystem::getGlobalComm() const{
  checkPreconditions( InitializationStage::GLOBAL_INIT );

  return globalComm_;
}


inline const CommunicatorType& MPISystem::getLocalComm() const{
  checkPreconditions( InitializationStage::LOCAL_INIT );

  return localComm_;
}


inline const CommunicatorType& MPISystem::getTeamComm() const{
  checkPreconditions( InitializationStage::TEAM_INIT );

  return teamComm_;
}


inline const CommunicatorType& MPISystem::getGlobalReduceComm() const{
  checkPreconditions( InitializationStage::ALL_INIT );

  return globalReduceComm_;
}

inline simft::Sim_FT_MPI_Comm MPISystem::getWorldCommFT(){
  checkPreconditionsFT( InitializationStage::WORLD_INIT );

  return worldCommFT_;
}

inline simft::Sim_FT_MPI_Comm MPISystem::getSpareCommFT(){
  checkPreconditionsFT( InitializationStage::WORLD_INIT );

  return spareCommFT_;
}

inline simft::Sim_FT_MPI_Comm MPISystem::getGlobalCommFT(){
  checkPreconditionsFT( InitializationStage::GLOBAL_INIT );

  return globalCommFT_;
}


inline simft::Sim_FT_MPI_Comm MPISystem::getLocalCommFT(){
  checkPreconditionsFT( InitializationStage::LOCAL_INIT );

  return localCommFT_;
}


inline simft::Sim_FT_MPI_Comm MPISystem::getTeamCommFT(){
  checkPreconditionsFT( InitializationStage::TEAM_INIT );

  return teamCommFT_;
}


inline simft::Sim_FT_MPI_Comm MPISystem::getGlobalReduceCommFT(){
  checkPreconditionsFT( InitializationStage::ALL_INIT );

  return globalReduceCommFT_;
}


inline const RankType& MPISystem::getWorldRank() const{
  checkPreconditions( InitializationStage::WORLD_INIT );

  return worldRank_;
}


inline const RankType& MPISystem::getGlobalRank() const{
  checkPreconditions( InitializationStage::GLOBAL_INIT );

  return globalRank_;
}


inline const RankType& MPISystem::getLocalRank() const{
  checkPreconditions( InitializationStage::LOCAL_INIT );

  return localRank_;
}


inline const RankType& MPISystem::getManagerRankWorld() const{
  checkPreconditions( InitializationStage::WORLD_INIT );

  return managerRankWorld_;
}


inline const RankType& MPISystem::getManagerRank() const{
  checkPreconditions( InitializationStage::GLOBAL_INIT );

  return managerRank_;
}


inline const RankType& MPISystem::getMasterRank() const{
  checkPreconditions( InitializationStage::LOCAL_INIT );

  return masterRank_;
}


inline bool MPISystem::isWorldManager() const{
  return ( worldRank_ == managerRankWorld_ );
}


inline bool MPISystem::isGlobalManager() const{
  return ( globalRank_ == managerRank_ && managerRank_ != MPI_UNDEFINED );
}


inline bool MPISystem::isMaster() const{
  return ( localRank_ == masterRank_ && masterRank_ != MPI_UNDEFINED );
}


inline bool MPISystem::isTeamLeader() const{
  return ( teamRank_ == teamLeaderRank_ && teamLeaderRank_ != MPI_UNDEFINED );
}


inline const RankType& MPISystem::getTeamColor() const{
  return teamColor_;
}


inline size_t MPISystem::getNumGroups() const{
  checkPreconditions( InitializationStage::WORLD_INIT );

  return nprocsByGroup_.size();
}


inline size_t MPISystem::getNumProcs( GroupType group ) const{
  checkPreconditions( InitializationStage::WORLD_INIT );

  return nprocsByGroup_[group];
}


inline RankType MPISystem::getGroupBaseWorldRank( GroupType group ) const {
  checkPreconditions( InitializationStage::WORLD_INIT );

  return groupBaseWorldRank_[group];
}


inline bool MPISystem::isInitialized() const{
  return initialized_ == InitializationStage::ALL_INIT;
}

inline void MPISystem::deleteReduceCommsFTAndComm() {
  deleteCommFTAndCcomm(&globalReduceCommFT_, &globalReduceComm_);
}

/*
// operators
std::ostream& operator<<(std::ostream& os, const MPISystem& ms);
std::ostream& operator<<(std::ostream& os, const MPISystemID& ms);
std::ostream& operator<<(std::ostream& os, const ConstMPISystemID& ms);
*/

} //namespace combigrid




#endif // MPISYSTEM_HPP
