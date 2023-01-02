#include <stdio.h>
#include <iostream>
#include "MPI-FT.h"
#include REAL_MPI_INCLUDE

int simft::Sim_FT_MPI_Send(void *buf, int count, MPI_Datatype type, int dest, int tag,
                           simft::Sim_FT_MPI_Comm f_comm) {
  simft::Sim_FT_decide_kill();
  simft::Sim_FT_Perform_background_operations();

  if (simft::Sim_FT_Check_comm_revoked(f_comm)) {
    return MPI_ERR_REVOKED;
  }

  int rc = MPI_SUCCESS;

  double Timeout_Begin = MPI_Wtime();

  // wait for request message

  char Alive_Status_Recv = 0;
  char Alive_Status_Send;
  Alive_Status_Send = I_AM_ALIVE;
  MPI_Request Send_Request;
  MPI_Status Recv_Status;
  MPI_Isend(&Alive_Status_Send, 1, MPI_CHAR, dest, tag + TAGOFFSETALIVE, f_comm->c_comm_copy_p2p,
            &Send_Request);
  MPI_Request_free(
      &Send_Request);  // We can't use MPI_Wait or MPI_Test => use MPI_Request_free instead
  bool leaveloop = false;
  while (!leaveloop) {
    int RecFlag = 0;

    while (RecFlag == 0) {
      MPI_Iprobe(dest, tag + TAGOFFSETALIVE, f_comm->c_comm_copy_p2p, &RecFlag, &Recv_Status);
      simft::Sim_FT_Perform_background_operations();
    }

    /*
     * If comm is revoked, a recent "Sim_FT_Perform_background_operations" could already have
     * responded the
     * alive-request from dest with "COMM_IS_REVOKED" => return error and do not try to receive
     * alive-message again
     */
    if (simft::Sim_FT_Check_comm_revoked(f_comm)) {
      return MPI_ERR_REVOKED;  // TODO: wait for incoming SIM_FT_STATUS_TAG, respond with
                               // comm_revoked
    } else {
      MPI_Recv(&Alive_Status_Recv, 1, MPI_CHAR, dest, tag + TAGOFFSETALIVE, f_comm->c_comm_copy_p2p,
               MPI_STATUS_IGNORE);
    }

    if (Alive_Status_Recv == PROBE_REQUEST) {
      MPI_Send(&count, 1, MPI_INT, dest, tag, f_comm->c_comm_copy_probe);
    } else {
      leaveloop = true;
    }
  }

  // check the content of the status message and react to it
  if (Alive_Status_Recv == I_AM_ALIVE || Alive_Status_Recv == I_AM_ALIVE_NB) {
    rc = MPI_Send(buf, count, type, dest, tag, f_comm->c_comm);

  } else if (Alive_Status_Recv == COMM_IS_REVOKED) {
    rc = MPI_ERR_REVOKED;

  } else {
    f_comm->dead_set.insert(dest);
    while (MPI_Wtime() - Timeout_Begin < SIM_FT_TIMEOUT) {
      // simulate timeout
      simft::Sim_FT_Perform_background_operations();
    }

    rc = MPI_ERR_PROC_FAILED;
  }

  return rc;
}
