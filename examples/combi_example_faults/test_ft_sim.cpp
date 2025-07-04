
#include <cassert>
#include <iostream>
#include <unistd.h>
#include "../../include/discotec/mpi/MPISystem.hpp"

using namespace combigrid;


void recoverComm(){
  std::cout << "rank " << theMPISystem()->getWorldRank() << " "
            << "recover communicator " << std::endl;

  // revoke commmworld
  MPI_Comm_revoke( theMPISystem()->getWorldCommFT() );

  // shrink world
  simft::Sim_FT_MPI_Comm newCommWorldFT;
  MPI_Comm_shrink( theMPISystem()->getWorldCommFT(), &newCommWorldFT );

  MPI_Comm newCommWorld = newCommWorldFT->c_comm;

  int newRank, newSize;
  MPI_Comm_rank( newCommWorld, &newRank );
  MPI_Comm_size( newCommWorld, &newSize );

  for( auto r=0; r < newSize; ++r ){
      if( newRank == r ){

        std::cout << "rank " << theMPISystem()->getWorldRank()
                  << "new rank " << newRank
                  << "new size " << newSize
                  << std::endl;
      }
      MPI_Barrier( newCommWorld );
  }
}

//start with e.g. mpirun.mpich -n 5 ./test_ft_sim 2
int main(int argc, char** argv) {
  assert( argc == 2 );

  //TODO this thread https://stackoverflow.com/questions/2642996/why-does-mpi-init-accept-pointers-to-argc-and-argv
  // indicates that this might not make sense before MPI_Init
  int rank_fail = atoi(argv[1]);

  //MPI_Init(&argc, &argv);
  simft::Sim_FT_MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  assert( size == 5 );

  theMPISystem()->init( 2, 2 );

  // worldComm
  // r 0,1 -> group 1
  // r 2,3 -> group 2
  // r 4 -> manager

  // globalComm
  // r 0 -> master 1
  // r 1 -> master 2
  // r 2 -> manager

  // manager code
  if( theMPISystem()->getWorldRank() == theMPISystem()->getManagerRankWorld() ){
    // wait for status from master 1
    int signal0;
    simft::Sim_FT_MPI_Request req0;
    simft::Sim_FT_MPI_Irecv( &signal0, 1, MPI_INT, 0, 0,
                             theMPISystem()->getGlobalCommFT(), &req0 );

    // wait for status from master 2
    int signal1;
    simft::Sim_FT_MPI_Request req1;
    simft::Sim_FT_MPI_Irecv( &signal1, 1, MPI_INT, 1, 0,
                             theMPISystem()->getGlobalCommFT(), &req1 );

    simft::Sim_FT_MPI_Status stat0;
    int err0 = simft::Sim_FT_MPI_Wait( &req0, &stat0 );
    if( err0 != 0 )
      signal0 = err0;

    simft::Sim_FT_MPI_Status stat1;
    int err1 = simft::Sim_FT_MPI_Wait( &req1, &stat1 );
    if( err1 != 0 )
      signal1 = err1;

    std::cout << "group 0 signal = " << signal0 << std::endl;
    std::cout << "group 1 signal = " << signal1 << std::endl;

    // if error send signal to recover communicators
    if( signal0 == 0 ){
      int command = 100;
      simft::Sim_FT_MPI_Send( &command, 1, MPI_INT, 0
                                    , 0, theMPISystem()->getGlobalCommFT() );
    }

    if( signal1 == 0 ){
      int command = 100;
      simft::Sim_FT_MPI_Send( &command, 1, MPI_INT, 1
                                    , 0, theMPISystem()->getGlobalCommFT() );
    }

    recoverComm();
  }

  // worker code
  else{
    if( theMPISystem()->getWorldRank() == rank_fail )
      simft::Sim_FT_kill_me();

    int err = simft::Sim_FT_MPI_Barrier( theMPISystem()->getLocalCommFT() );

    // send status signal
    MASTER_EXCLUSIVE_SECTION{
      simft::Sim_FT_MPI_Send( &err, 1, MPI_INT, theMPISystem()->getManagerRank()
                              , 0, theMPISystem()->getGlobalCommFT() );
    }

    // when local error detected sleep until recovery, nothing better to
    // do anyway
    if( err != 0){
      recoverComm();

    } else {
      // recv command to recover communicator
      int command;
      MASTER_EXCLUSIVE_SECTION{
        simft::Sim_FT_MPI_Status stat;
        simft::Sim_FT_MPI_Recv( &command, 1, MPI_INT, theMPISystem()->getManagerRank()
                                      , 0, theMPISystem()->getGlobalCommFT(), &stat );
      }

      simft::Sim_FT_MPI_Bcast(
          &command, 1, MPI_INT, theMPISystem()->getMasterRank(),
          theMPISystem()->getLocalCommFT() );

      recoverComm();
    }
  }

  simft::Sim_FT_MPI_Finalize();
}
