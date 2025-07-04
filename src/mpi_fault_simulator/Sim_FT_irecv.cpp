#include "../../include/discotec/MPI-FT.h"
#include REAL_MPI_INCLUDE

#include <iostream>

int simft::Sim_FT_MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
                            simft::Sim_FT_MPI_Comm f_comm, simft::Sim_FT_MPI_Request *request) {
  simft::Sim_FT_decide_kill();

  if (simft::Sim_FT_Check_comm_revoked(f_comm)) {
    return MPI_ERR_REVOKED;
  }

  MPI_Request Send_Request;
  char Alive_Status_Send = I_AM_ALIVE_NB;

  // an alive-message can only be sent to a specific process. If source is MPI_ANY_SOURCE, send
  // alive-message in MPI_TEST
  if (source != MPI_ANY_SOURCE) {
    MPI_Isend(&Alive_Status_Send, 1, MPI_CHAR, source, tag + TAGOFFSETALIVE,
              f_comm->c_comm_copy_p2p, &Send_Request);
    MPI_Request_free(&Send_Request);
  }

  (*request) = new simft::Sim_FT_MPI_Request_struct;
  (*request)->P2P_Stats.p2p_rank = source;
  (*request)->f_comm = f_comm;
  (*request)->Request_Type = 1;  // set request type to p2p
  (*request)->P2P_Stats.tag = tag + TAGOFFSETALIVE;
  (*request)->P2P_Stats.Type = 11;

  return MPI_Irecv(buf, count, datatype, source, tag, f_comm->c_comm, &(*request)->c_request);
}
