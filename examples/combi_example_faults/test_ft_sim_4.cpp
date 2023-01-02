#include "mpi.h"
#include <cassert>
#include <iostream>
#include <unistd.h>
#include "mpi_fault_simulator/MPI-FT.h"

int main(int argc, char** argv) {
  assert( argc == 2 );
  int rank_fail = atoi(argv[1]);

  //MPI_Init(&argc, &argv);
  simft::Sim_FT_MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  assert( size > 1 );

  simft::Sim_FT_MPI_Comm cworld = simft::Sim_FT_MPI_COMM_WORLD;

  // process 0 waits for message from each other process
  if( rank == 0 ){
    std::vector< simft::Sim_FT_MPI_Request > requests( size );
    std::vector< int > signals( size );
    for( auto r = 1; r < size; ++r ){
      simft::Sim_FT_MPI_Irecv( &signals[r], 1, MPI_INT, r, 0, cworld,
                               &requests[r] );
    }

    sleep(1);

    std::vector< int > flags( size );
    std::vector< int > errors( size );
    std::vector< simft::Sim_FT_MPI_Status > status( size );
    for( auto r = 1; r < size; ++r ){
      errors[r] = simft::Sim_FT_MPI_Test( &requests[r], &flags[r],
                                              &status[r] );
    }

    for( auto r = 1; r < size; ++r ){
      std::cout << "r = " << r << " "
                << "signal = " << signals[r] << " "
                << "flag = " << flags[r] << " "
                << "error = " << errors[r] << " "
                << std::endl;
    }

    // check for errors and revoke communicator if necessary
    int err = 0;
    for( auto r = 1; r < size; ++r ){
      if( errors[r] != 0 ){
        err = 1;
        MPI_Comm_revoke( cworld );
        break;
      }
    }

    // notify alive process about erroneous comm. this will trigger shrinking
    // the communicator
    for( auto r = 1; r < size; ++r ){
      if( errors[r] == 0 )
        MPI_Send( &err, 1, MPI_INT, r, 0, cworld->c_comm );
    }

    if( err == 1 ){
      simft::Sim_FT_MPI_Comm newcomm;
      MPI_Comm_shrink( cworld, &newcomm );
      cworld = newcomm;
    }

    MPI_Barrier( cworld->c_comm );
  }

  if( rank > 0 ){
    int signal = 100;

    if( rank == rank_fail )
      simft::Sim_FT_kill_me();

    simft::Sim_FT_MPI_Send( &signal, 1, MPI_INT, 0, 0, cworld );

    // recv error signal
    int err = -1;
    MPI_Recv( &err, 1, MPI_INT, 0, 0, cworld->c_comm, MPI_STATUS_IGNORE );

    if( err == 1 ){
      simft::Sim_FT_MPI_Comm newcomm;
      MPI_Comm_shrink( cworld, &newcomm );
      cworld = newcomm;
    }

    MPI_Barrier( cworld->c_comm );
  }

  int nrank, nsize;
  MPI_Comm_rank( cworld->c_comm, &nrank );
  MPI_Comm_size( cworld->c_comm, &nsize );

  for( auto r=0; r < nsize; ++r ){
    if( nrank == r ){
      std::cout << "old rank = " << rank << " "
                << "new rank = " << nrank << " "
                << "old size = " << size << " "
                << "new size = " << nsize
                << std::endl;
    }

    MPI_Barrier( cworld->c_comm );
  }

  simft::Sim_FT_MPI_Finalize();
}
