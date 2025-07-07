#include "mpi.h"
#include <cassert>
#include <iostream>
#include <unistd.h>
#include "discotec/MPI-FT.h"

void createCommunicators( size_t ngroup, size_t nprocs, int globalID, int globalSize,
    int& managerIDgcomm, int& gglobalID, int& lglobalID, MPI_Comm& gcomm, MPI_Comm& lcomm){
  /* determine global globalID of each process
   * the manager process always has the highest globalID
   * all other processes are worker processes */


  /* create a local communicator for each process group
   * lcomm is the local communicator of its own process group for each worker process
   * for manager, lcomm is a group which contains only manager process and can be ignored
   */
  int color = globalID / nprocs;
  int key = globalID - color * nprocs;
  MPI_Comm_split(MPI_COMM_WORLD, color, key, &lcomm);
  MPI_Comm_rank(lcomm, &lglobalID);
  const int managerIDworld = globalSize - 1;

  /* create global communicator which contains only the manager and the master
   * process of each process group
   * the master processes of the process groups are the processes which have
   * globalID 0 in lcomm
   * this communicator is used for communication between master processes of the
   * process groups and the manager and the master processes to each other
   */
  MPI_Group worldGroup;
  MPI_Comm_group(MPI_COMM_WORLD, &worldGroup);

  std::vector<int> globalIDs(ngroup + 1);
  for (size_t i = 0; i < ngroup; i++) {
    globalIDs[i] = i * nprocs;
  }
  globalIDs.back() = managerIDworld;

  MPI_Group rootGroup;
  MPI_Group_incl(worldGroup, (int) globalIDs.size(), &globalIDs[0], &rootGroup);

  MPI_Comm_create(MPI_COMM_WORLD, rootGroup, &gcomm);

  managerIDgcomm = MPI_PROC_NULL;
  if (gcomm != MPI_COMM_NULL) {
    int gcommSize;
    MPI_Comm_size(gcomm, &gcommSize);
    managerIDgcomm = gcommSize - 1;
    MPI_Comm_rank(gcomm, &globalID);
  }
}


int main(int argc, char** argv) {
  //MPI_Init(&argc, &argv);
  simft::Sim_FT_MPI_Init(&argc, &argv);

  int globalID, size;
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &globalID );

  assert( argc == 2 );
  int rank_fail = atoi(argv[1]);

  if( globalID == rank_fail ){
    std::cout << "globalID " << rank_fail << " will fail" << std::endl;
  }

  MPI_Comm normalComm = MPI_COMM_WORLD;
  simft::Sim_FT_MPI_Comm ftComm = simft::Sim_FT_MPI_COMM_WORLD;

  const int managerIDworld = size - 1;

  size_t ngroup = 1;
  size_t nprocs = 1;
  assert( size  == ngroup*nprocs + 1 );

  int grank = -1;
  int lrank = -1;
  int managerIDgcomm = -1;
  int lglobalID = -1;
  int gglobalID = -1;
  int globalSize = size;

  MPI_Comm gcomm;
  MPI_Comm lcomm;
  createCommunicators(ngroup, nprocs, globalID, globalSize,
                      managerIDgcomm, grank, lrank, gcomm, lcomm);

  for( auto r=0; r<size; ++r ){
    if( globalID == r ){
      std::cout << "globalID = " << globalID << "\t"
                << "gcomm = " << gcomm << "\t"
                << "grank = " << grank << "\t"
                << "lcomm = " << lcomm << "\t"
                << "lrank = " << lrank << "\t"
                << "mgrID = " << managerIDgcomm << "\t"
                << std::endl;
    }
    MPI_Barrier( MPI_COMM_WORLD );
  }

  // create wrapper for gcomm
  simft::Sim_FT_MPI_Comm gcomm_ft = new simft::Sim_FT_MPI_Comm_struct;
  gcomm_ft->c_comm = gcomm;
  Sim_FT_Initialize_new_comm( &gcomm_ft, true );

  // manager code
  if( globalID == managerIDworld ){
    // sleep for some time
    sleep( 1 );

    int signal1(-1);
    simft::Sim_FT_MPI_Request req1;

    simft::Sim_FT_MPI_Irecv( &signal1, 1, MPI_INT, 0, 0, simft::Sim_FT_MPI_COMM_WORLD, &req1 );

    // in case of failure this should be detectable here
    int flag1 = -1;
    simft::Sim_FT_MPI_Status status;
    int err1 = simft::Sim_FT_MPI_Test( &req1, &flag1, (simft::Sim_FT_MPI_Status *)1 );

    std::cout << "signal = " << signal1 << "flag = " << flag1 << " err =" << err1 << std::endl;
  } else{
    if( lrank == 0 ){
      int signal = 100;

      if( globalID == rank_fail )
        simft::Sim_FT_kill_me();

      simft::Sim_FT_MPI_Send( &signal, 1, MPI_INT, 1, 0, simft::Sim_FT_MPI_COMM_WORLD );
    }
  }

  simft::Sim_FT_MPI_Finalize();
}
