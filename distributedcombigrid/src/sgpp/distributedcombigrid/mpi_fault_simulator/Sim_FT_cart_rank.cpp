/*
 *  Sim_FT_cart_rank.cpp
 *
 *  Created on: 20.10.2015
 *      Author: Johannes Walter
 */
#include "MPI-FT.h"
#include REAL_MPI_INCLUDE

int simft::Sim_FT_MPI_Cart_rank(simft::Sim_FT_MPI_Comm f_comm, int coords[], int *rank) {
  simft::Sim_FT_decide_kill();
  simft::Sim_FT_Perform_background_operations();

  return MPI_Cart_rank(f_comm->c_comm, coords, rank);
}
