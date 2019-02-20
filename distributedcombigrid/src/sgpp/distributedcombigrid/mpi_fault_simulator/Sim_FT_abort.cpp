/*
 *  Sim_FT_abort.cpp
 *
 *  Created on: 31.07.2015
 *      Author: Johannes Walter
 */
#include "MPI-FT.h"
#include REAL_MPI_INCLUDE

int simft::Sim_FT_MPI_Abort(simft::Sim_FT_MPI_Comm f_comm, int errorcode) {
  return MPI_Abort(f_comm->c_comm, errorcode);
}
