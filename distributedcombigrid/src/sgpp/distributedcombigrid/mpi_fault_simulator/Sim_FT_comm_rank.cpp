/*
 *  Sim_FT_comm_rank.cpp
 *
 *  Created on: 31.07.2015
 *      Author: Johannes Walter
 */
#include "MPI-FT.h"
#include REAL_MPI_INCLUDE
#include <iostream>

int simft::Sim_FT_MPI_Comm_rank(simft::Sim_FT_MPI_Comm f_comm, int *rank) {
  simft::Sim_FT_decide_kill();
  simft::Sim_FT_Perform_background_operations();

  return MPI_Comm_rank(f_comm->c_comm, rank);
}
