#include "../../include/discotec/MPI-FT.h"
#ifndef DISABLE_NBC
#include REAL_MPI_INCLUDE

int simft::Sim_FT_MPI_Ibarrier(simft::Sim_FT_MPI_Comm f_comm, simft::Sim_FT_MPI_Request *request) {
  simft::Sim_FT_decide_kill();
  if (simft::Sim_FT_Check_comm_revoked(f_comm)) {
    return MPI_ERR_REVOKED;
  }

  // create new Sim_FT_Istats object and store necessary information
  simft::Sim_FT_Istats NewCollective;
  NewCollective.Type = 0;  // ibarrier

  (*request) = new simft::Sim_FT_MPI_Request_struct;
  NewCollective.request = *request;
  f_comm->Recent_Icollectives.push_back(NewCollective);
  (*request)->f_comm = f_comm;

  return MPI_Ibarrier(f_comm->c_comm, &(*request)->c_request);
}
#endif
