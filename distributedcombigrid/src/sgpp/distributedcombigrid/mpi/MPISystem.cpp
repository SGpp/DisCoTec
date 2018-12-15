/*
 * MPISetup.cpp
 *
 *  Created on: Jan 23, 2013
 *      Author: mh
 *
 *  Partially copied from the pe Physics Engine class MPISystem
 */

#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/utils/Stats.hpp"
#include <iostream>

namespace combigrid {

/*!\brief Constructor for the MPISystem class.
 //
 // \exception std::runtime_error MPI was not initialized
 //
 // Constructor for the MPI System. The default global communicator and local communicator is MPI_COMM_WORLD.
 // The total number of MPI processes and the rank of the MPI process in is determined from
 // the communicator. The default root process has rank zero.
 */
MPISystem::MPISystem() :
    initialized_(false),
    worldComm_(MPI_COMM_NULL),
    globalComm_(MPI_COMM_NULL),
    localComm_(MPI_COMM_NULL),
    globalReduceComm_(MPI_COMM_NULL),
    worldRank_(MPI_UNDEFINED),
    globalRank_(MPI_UNDEFINED),
    localRank_(MPI_UNDEFINED),
    globalReduceRank_(MPI_UNDEFINED),
    managerRank_(MPI_UNDEFINED),
    managerRankWorld_(MPI_UNDEFINED),
    masterRank_(MPI_UNDEFINED)
{
  // check if MPI was initialized (e.g. by MPI_Init or similar)
  int mpiInitialized(0);
  MPI_Initialized(&mpiInitialized);
  assert( mpiInitialized && "MPI is not initialized! Call MPI_Init first." );
}

/*!\brief Destructor for the MPISystem class.
 */
MPISystem::~MPISystem() {
  // todo: the fault tolerant communicator are initialized with new -> delete
}


void MPISystem::init( size_t ngroup, size_t nprocs ){
  assert( !initialized_ && "MPISystem already initialized!" );

  ngroup_ = ngroup;
  nprocs_ = nprocs;

  worldComm_ = MPI_COMM_WORLD;

  /* init worldComm
   * the manager has highest rank here
   */
  int worldSize;
  MPI_Comm_size( worldComm_, &worldSize );
  assert( worldSize == int(ngroup_ * nprocs_ + 1) );

  MPI_Comm_rank( worldComm_, &worldRank_ );
  managerRankWorld_ = worldSize - 1;

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
void MPISystem::init( size_t ngroup, size_t nprocs, CommunicatorType lcomm ){
  assert( !initialized_ && "MPISystem already initialized!" );

  ngroup_ = ngroup;
  nprocs_ = nprocs;

  worldComm_ = MPI_COMM_WORLD;

  /* init worldComm
   * the manager has highest rank here
   */
  int worldSize;
  MPI_Comm_size( worldComm_, &worldSize );
  assert( worldSize == int(ngroup_ * nprocs_ + 1) );

  MPI_Comm_rank( worldComm_, &worldRank_ );
  managerRankWorld_ = worldSize - 1;

  /* init localComm
   * lcomm is the local communicator of its own process group for each worker process.
   * for manager, lcomm is a group which contains only manager process and can be ignored
   */
  // manager is not supposed to have a localComm
  if( worldRank_ == managerRankWorld_ )
    localComm_ = MPI_COMM_NULL;
  else{
    localComm_ = lcomm;

    // todo: think through which side effects changing the master rank would have
    // in principle this does not have to be 0
    const int masterRank = 0;

    int localSize;
    MPI_Comm_size( localComm_, &localSize );
    assert( masterRank < localSize );

    masterRank_ = masterRank;

    MPI_Comm_rank( localComm_, &localRank_ );
  }

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
void MPISystem::init( CommunicatorType wcomm, size_t ngroup, size_t nprocs ){
  ngroup_ = ngroup;
  nprocs_ = nprocs;

  worldComm_ = wcomm;

  /* init worldComm
   * the manager has highest rank here
   */
  int worldSize;
  MPI_Comm_size( worldComm_, &worldSize );
  assert( worldSize == int(ngroup_ * nprocs_ + 1) );

  MPI_Comm_rank( worldComm_, &worldRank_ );
  managerRankWorld_ = worldSize - 1;

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

void MPISystem::initLocalComm(){
  int color = worldRank_ / int(nprocs_);
  int key = worldRank_ - color * int(nprocs_);
  MPI_Comm_split( worldComm_, color, key, &localComm_ );

  /* set group number in Stats. this is necessary for postprocessing */
  Stats::setAttribute("group", std::to_string(color));

  // manager is not supposed to have a localComm
  if( worldRank_ == managerRankWorld_ )
    localComm_ = MPI_COMM_NULL;
  else{
    // todo: think through which side effects changing the master rank would have
    // in principle this does not have to be 0
    const int masterRank = 0;

    int localSize;
    MPI_Comm_size( localComm_, &localSize );
    assert( masterRank < localSize );

    masterRank_ = masterRank;

    MPI_Comm_rank( localComm_, &localRank_ );
  }
}


void MPISystem::initGlobalComm(){
  MPI_Group worldGroup;
  MPI_Comm_group( worldComm_, &worldGroup);

  std::vector<int> ranks( ngroup_ + 1 );
  for (size_t i = 0; i < ngroup_; i++) {
    ranks[i] = int( i * nprocs_ );
  }
  ranks.back() = managerRankWorld_;

  MPI_Group globalGroup;
  MPI_Group_incl( worldGroup, int( ranks.size() ), &ranks[0], &globalGroup );

  MPI_Comm_create( worldComm_, globalGroup, &globalComm_ );

  if( globalComm_ != MPI_COMM_NULL ) {
    int globalSize;
    MPI_Comm_size( globalComm_, &globalSize );

    managerRank_ = globalSize - 1;

    MPI_Comm_rank( globalComm_, &globalRank_ );
  }

  /* mark master processes and manager process in Stats. this is necessary for postprocessing */
  Stats::setAttribute("group_manager", std::to_string(globalComm_ != MPI_COMM_NULL));
}

void MPISystem::initGlobalReduceCommm() {
  // create communicator which only contains workers
  MPI_Comm workerComm;
  {
    int color = ( worldRank_ != managerRankWorld_ ) ? 1 : 0;
    int key = (worldRank_ != managerRankWorld_ ) ? worldRank_ : 0;
    MPI_Comm_split( worldComm_, color, key, &workerComm);
  }

  if( worldRank_ != managerRankWorld_ ) {
    int workerID;
    MPI_Comm_rank(workerComm, &workerID);

    MPI_Comm globalReduceComm;
    int color = workerID % int(nprocs_);
    int key = workerID / int(nprocs_);
    MPI_Comm_split(workerComm, color, key, &globalReduceComm);

    globalReduceComm_ = globalReduceComm;

    MPI_Comm_rank( globalReduceComm_, &globalReduceRank_ );
  }
}

} // namespace combigrid


