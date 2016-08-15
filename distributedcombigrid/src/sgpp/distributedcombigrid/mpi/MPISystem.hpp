#ifndef MPISYSTEM_HPP
#define MPISYSTEM_HPP

#include <mpi.h>
#include <ostream>
#include <vector>
#include <assert.h>

#include "sgpp/distributedcombigrid/mpi/MPISystemID.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"

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
class MPISystem {
 public:
  ~MPISystem();

  MPISystem( MPISystem const & ) = delete;

  MPISystem& operator=( MPISystem const & ) = delete;

  void init( size_t ngroups, size_t nprocs );

  inline const CommunicatorType& getWorldComm() const;

  inline const CommunicatorType& getGlobalComm() const;

  inline const CommunicatorType& getLocalComm() const;

  inline const CommunicatorType& getGlobalReduceComm() const;

  inline const RankType& getWorldRank() const;

  inline const RankType& getGlobalRank() const;

  inline const RankType& getLocalRank() const;

  inline const RankType& getManagerRankWorld() const;

  inline const RankType& getManagerRank() const;

  inline const RankType& getMasterRank() const;

  inline bool isWorldManager() const;

  inline bool isGlobalManager() const;

  inline bool isMaster() const;

  inline size_t getNumGroups() const;

  inline size_t getNumProcs() const;

 private:
  explicit MPISystem();

  friend MPISystemID theMPISystem();

  inline void checkPreconditions() const;

  inline void checkPreconditionsFT() const;

  /* create the global communicators for the global reduce.
   * all processes which have local rank i in their process group will be grouped in a
   * distinct communicator.
   * this will only work if all pgroups have same size and dsg uses the same assignment
   * of subspaces to processes in each pgroup
   */
  void initGlobalReduceCommm();

  void initLocalComm();

  void initGlobalComm();

  bool initialized_;

  CommunicatorType worldComm_;

  CommunicatorType globalComm_;

  CommunicatorType localComm_;

  CommunicatorType globalReduceComm_;

  RankType worldRank_;

  RankType globalRank_;

  RankType localRank_;

  RankType globalReduceRank_;

  RankType managerRank_;

  RankType managerRankWorld_;

  RankType masterRank_;

  size_t ngroup_;

  size_t nprocs_;
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


/*
// operators
std::ostream& operator<<(std::ostream& os, const MPISystem& ms);
std::ostream& operator<<(std::ostream& os, const MPISystemID& ms);
std::ostream& operator<<(std::ostream& os, const ConstMPISystemID& ms);
*/

} //namespace combigrid




#endif // MPISYSTEM_HPP
