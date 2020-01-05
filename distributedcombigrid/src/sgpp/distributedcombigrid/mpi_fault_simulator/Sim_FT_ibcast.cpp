#include "MPI-FT.h"
#ifndef DISABLE_NBC
#include REAL_MPI_INCLUDE

int simft::Sim_FT_MPI_Ibcast(void *buffer, int count, MPI_Datatype datatype, int root,
                             simft::Sim_FT_MPI_Comm f_comm, simft::Sim_FT_MPI_Request *request) {
  simft::Sim_FT_decide_kill();
  if (simft::Sim_FT_Check_comm_revoked(f_comm)) {
    return MPI_ERR_REVOKED;
  }

  // create new Sim_FT_Istats object and store necessary information
  simft::Sim_FT_Istats NewCollective;
  NewCollective.Type = 1;  // ibcast
  NewCollective.Datatype = simft::Datatype_to_Dtype(
      datatype);  //"transform" MPI datatype to int, so we can store and send it
  NewCollective.root = root;
  NewCollective.count = count;

  // create nem custom struct object, set the object in NewCollective and attach NewCollective to
  // Recent_Icollectives
  (*request) = new simft::Sim_FT_MPI_Request_struct;
  NewCollective.request = *request;
  f_comm->Recent_Icollectives.push_back(NewCollective);
  (*request)->f_comm = f_comm;

  return MPI_Ibcast(buffer, count, datatype, root, f_comm->c_comm, &(*request)->c_request);
}
#endif
