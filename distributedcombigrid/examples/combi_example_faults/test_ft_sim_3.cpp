#include "mpi.h"
#include <cassert>
#include <iostream>
#include <unistd.h>
#include "sgpp/distributedcombigrid/mpi_fault_simulator/MPI-FT.h"

int main(int argc, char** argv) {
  //MPI_Init(&argc, &argv);
  simft::Sim_FT_MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  assert( size == 4 );

  // with simft::Sim_FT_MPI_COMM_WORLD everything is fine
  {
    if( rank == 0 ){
      int signal = -1;
      simft::Sim_FT_MPI_Request req;

      simft::Sim_FT_MPI_Irecv( &signal, 1, MPI_INT, 1, 0,
                               simft::Sim_FT_MPI_COMM_WORLD, &req );

      sleep(1);

      int flag = -1;
      simft::Sim_FT_MPI_Status stat;
      int err = -1;
      while( flag != 1 ) err = simft::Sim_FT_MPI_Test( &req, &flag, &stat );

      std::cout << "rank 0 recv: "
                << "signal = " << signal << "\t"
                << "flag = " << flag << "\t"
                << "err = " << err << std::endl;
    }

    if( rank == 1 ){
      int signal = 100;

      //simft::Sim_FT_kill_me();

      simft::Sim_FT_MPI_Send( &signal, 1, MPI_INT, 0, 0,
                              simft::Sim_FT_MPI_COMM_WORLD );
    }
  }


  //now lets try a custom communicator
  {
    // split comm world into two group
    int nprocs = size / 2;
    int color = rank / nprocs;
    int key = rank - color * nprocs;
    MPI_Comm lcomm;
    int lrank;
    MPI_Comm_split(MPI_COMM_WORLD, color, key, &lcomm);
    MPI_Comm_rank( lcomm, &lrank );

    // create wrapper for lcomm
    simft::Sim_FT_MPI_Comm lcomm_ft = new simft::Sim_FT_MPI_Comm_struct;
    lcomm_ft->c_comm = lcomm;
    Sim_FT_Initialize_new_comm( &lcomm_ft, true );

    for( auto r=0; r<size; ++r ){
      if( rank == r ){
        std::cout << "rank = " << rank << "\t"
                  << "lrank = " << lrank << "\t"
                  << std::endl;
      }
      MPI_Barrier( MPI_COMM_WORLD );
    }

    if( lrank == 0 ){
      int signal = -1;
      simft::Sim_FT_MPI_Request req;

      simft::Sim_FT_MPI_Irecv( &signal, 1, MPI_INT, 1, 0,
                               lcomm_ft, &req );

      sleep(1);

      int flag = -1;
      simft::Sim_FT_MPI_Status stat;
      int err = -1;
      while( flag != 1 ) err = simft::Sim_FT_MPI_Test( &req, &flag, &stat );

      std::cout << "lrank 0 recv: "
                << "signal = " << signal << "\t"
                << "flag = " << flag << "\t"
                << "err = " << err << std::endl;
    }

    if( lrank == 1 ){
      int signal = 100;

      //simft::Sim_FT_kill_me();

      simft::Sim_FT_MPI_Send( &signal, 1, MPI_INT, 0, 0, lcomm_ft );
    }

    MPI_Barrier( MPI_COMM_WORLD );

    // the other way around, lrank 0 sends to lrank 1
    if( lrank == 1 ){
      int signal = -1;
      simft::Sim_FT_MPI_Request req;

      simft::Sim_FT_MPI_Irecv( &signal, 1, MPI_INT, 0, 0,
                               lcomm_ft, &req );

      sleep(1);

      int flag = -1;
      simft::Sim_FT_MPI_Status stat;
      int err = -1;
      while( flag != 1 ) err = simft::Sim_FT_MPI_Test( &req, &flag, &stat );

      std::cout << "lrank 1 recv: "
                << "signal = " << signal << "\t"
                << "flag = " << flag << "\t"
                << "err = " << err << std::endl;
    }

    if( lrank == 0 ){
      int signal = 100;

      simft::Sim_FT_kill_me();

      simft::Sim_FT_MPI_Send( &signal, 1, MPI_INT, 1, 0, lcomm_ft );
    }
  }




  simft::Sim_FT_MPI_Finalize();
}
