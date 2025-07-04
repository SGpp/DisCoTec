#include "discotec/MPI-FT.h"
#include REAL_MPI_INCLUDE
#include <iostream>
#include "discotec/utils/Stats.hpp"

simft::Sim_FT_MPI_Comm simft::Sim_FT_MPI_COMM_WORLD;
simft::Sim_FT_MPI_Comm simft::Sim_FT_MPI_COMM_NULL;

simft::Sim_FT_MPI_Request_struct simft::Sim_FT_MPI_REQUEST_NULL;
std::set<simft::Sim_FT_MPI_Comm>
    simft::Sim_FT_Current_Active_Communicators;  //= { simft::Sim_FT_MPI_COMM_WORLD };
std::map<int, simft::Sim_FT_P2P_Response> simft::Sim_FT_Current_Active_P2P_Requests;

#ifdef USE_NON_BLOCKING_COLLECTIVES
MPI_Op simft::customOp;
#endif

int simft::Sim_FT_MPI_Init(int *argc, char ***argv) {
  int ret = MPI_Init(argc, argv);
  
  simft::Sim_FT_MPI_COMM_WORLD = new simft::Sim_FT_MPI_Comm_struct; //TODO make_unique
  simft::Sim_FT_MPI_COMM_WORLD->c_comm = MPI_COMM_WORLD;

  simft::Sim_FT_MPI_COMM_NULL = new simft::Sim_FT_MPI_Comm_struct;

  simft::Sim_FT_MPI_REQUEST_NULL.c_request = MPI_REQUEST_NULL;

  // int worldSize;
  // MPI_Comm_size(MPI_COMM_WORLD, &worldSize );
  // simft::Sim_FT_MPI_COMM_WORLD->Root_Rank = worldSize - 1;
  // important: initialize MPI_COMM_WORLD in our fault layer
  simft::Sim_FT_Initialize_new_comm(&simft::Sim_FT_MPI_COMM_WORLD, true);

  // combigrid::Stats::initialize();

  return ret;
}

void simft::Sim_FT_MPI_Init_worker() {
  simft::Sim_FT_MPI_COMM_WORLD = new simft::Sim_FT_MPI_Comm_struct;
  simft::Sim_FT_MPI_COMM_WORLD->c_comm = MPI_COMM_WORLD;

  simft::Sim_FT_MPI_COMM_NULL = new simft::Sim_FT_MPI_Comm_struct;

  simft::Sim_FT_MPI_REQUEST_NULL.c_request = MPI_REQUEST_NULL;

  // int ret = MPI_Init(argc, argv); MPI is initialized beforehand
  // int worldSize;
  // MPI_Comm_size(MPI_COMM_WORLD, &worldSize );
  // simft::Sim_FT_MPI_COMM_WORLD->Root_Rank = worldSize - 1;
  // important: initialize MPI_COMM_WORLD in our fault layer
  simft::Sim_FT_Initialize_new_comm(&simft::Sim_FT_MPI_COMM_WORLD, true);
  // combigrid::Stats::initialize();

  // std::cout << "Init FT_MPI worker!";
}
