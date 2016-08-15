/*
 * MPISetup.cpp
 *
 *  Created on: Jan 23, 2013
 *      Author: mh
 *
 *  Partially copied from the pe Physics Engine class MPISystem
 */

#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/utils/StatsContainer.hpp"
#include <iostream>

namespace{
  using namespace combigrid;

  std::string getMinMaxAvg( RankType r, int size, std::string timerName, bool isTimer,
                            MPI_Comm comm ){
    double value;
    if( isTimer )
      value = theStatsContainer()->getDuration( timerName );
    else
      value = theStatsContainer()->getValue( timerName );

    double min, max, sum;
    MPI_Reduce( &value, &min, 1, MPI_DOUBLE, MPI_MIN, r, comm );
    MPI_Reduce( &value, &max, 1, MPI_DOUBLE, MPI_MAX, r, comm );
    MPI_Reduce( &value, &sum, 1, MPI_DOUBLE, MPI_SUM, r, comm );

    double avg = sum / static_cast<double>( size );

    std::stringstream ss;
    ss << timerName << "\t\t" << min << "\t" << max << "\t" << avg;

    return ss.str();
  }
} // anonymous namespace


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
  theStatsContainer()->setTimerStart("init-local");
  initLocalComm();
  theStatsContainer()->setTimerStop("init-local");

  /* create global communicator which contains only the manager and the master
   * process of each process group
   * the master processes of the process groups are the processes which have
   * rank 0 in lcomm
   * this communicator is used for communication between master processes of the
   * process groups and the manager and the master processes to each other
   */
  theStatsContainer()->setTimerStart("init-global");
  initGlobalComm();
  theStatsContainer()->setTimerStop("init-global");

  theStatsContainer()->setTimerStart("init-global-reduce");
  initGlobalReduceCommm();
  theStatsContainer()->setTimerStop("init-global-reduce");

  initialized_ = true;
}


void MPISystem::initLocalComm(){
  theStatsContainer()->setTimerStart("init-local-split");
  int color = worldRank_ / int(nprocs_);
  int key = worldRank_ - color * int(nprocs_);
  MPI_Comm_split( worldComm_, color, key, &localComm_ );
  theStatsContainer()->setTimerStop("init-local-split");

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
  }
}


void MPISystem::recoverCommunicators( bool groupAlive ){
  // revoke commmworld
  theStatsContainer()->setTimerStart("recoverComm-revoke");
  //WORLD_MANAGER_EXCLUSIVE_SECTION{
    MPI_Comm_revoke( theMPISystem()->getWorldComm() );
  //}
  theStatsContainer()->setTimerStop("recoverComm-revoke");

  // shrink world
  theStatsContainer()->setTimerStart("recoverComm-shrink");
  MPI_Comm newCommWorld;
  MPI_Comm_shrink( theMPISystem()->getWorldComm(), &newCommWorld );
  theStatsContainer()->setTimerStop("recoverComm-shrink");

  // split off alive procs. this will be the new WorldComm
  // processes of dead groups set color to MPI_UNDEFINED. in this case
  // MPI_Comm_split returns MPI_COMMM_NULL
  int color = (groupAlive) ? 1 : MPI_UNDEFINED;
  int key = worldRank_;
  MPI_Comm_split( newCommWorld, color, key, &worldComm_ );

  // early exit for procs in process group that has been split off
  if( worldComm_ == MPI_COMM_NULL)
    return;

  // todo: remove
  // output new commWorld
  {
    int newRank, newSize;
    MPI_Comm_rank( worldComm_, &newRank );
    MPI_Comm_size( worldComm_, &newSize );

    if( newRank == 0 )
      std::cout << "new WorldComm:" << std::endl;
    for( auto r=0; r < newSize; ++r ){
        if( newRank == r ){
          std::cout << "rank " << theMPISystem()->getWorldRank()
                    << " new rank " << newRank
                    << " new size " << newSize
                    << std::endl;
        }
        MPI_Barrier( worldComm_ );
    }
  }

  int worldSize;
  MPI_Comm_size( worldComm_, &worldSize );
  assert( (worldSize - 1) % nprocs_ == 0 );
  ngroup_ = (worldSize - 1) / nprocs_;

  MPI_Comm_rank( worldComm_, &worldRank_ );
  managerRankWorld_ = worldSize - 1;

  //initLocalComm();

  initGlobalComm();

  initGlobalReduceCommm();

  /* print stats */
   int ngroup( theMPISystem()->getNumGroups() );
   int nprocs( theMPISystem()->getNumProcs() );
   std::string t_revoke = getMinMaxAvg( theMPISystem()->getManagerRankWorld(),
                                        ngroup*nprocs+1,
                                      "recoverComm-revoke", true,
                                      theMPISystem()->getWorldComm() );
   std::string t_shrink = getMinMaxAvg( theMPISystem()->getManagerRankWorld(),
                                        ngroup*nprocs+1,
                                      "recoverComm-shrink", true,
                                      theMPISystem()->getWorldComm() );
   WORLD_MANAGER_EXCLUSIVE_SECTION{ std::cout << t_revoke << std::endl; }
   WORLD_MANAGER_EXCLUSIVE_SECTION{ std::cout << t_shrink << std::endl; }
}

} // namespace combigrid



