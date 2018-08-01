#include <stdio.h>
#include "MPI-FT.h"
#include REAL_MPI_INCLUDE
#include <iostream>

int simft::Sim_FT_MPI_Recv(void *buf, int count, MPI_Datatype type, int source, int tag,
                           simft::Sim_FT_MPI_Comm f_comm, simft::Sim_FT_MPI_Status *status) {
  simft::Sim_FT_decide_kill();
  simft::Sim_FT_Perform_background_operations();

  if (simft::Sim_FT_Check_comm_revoked(f_comm)) {
    return MPI_ERR_REVOKED;
  }

  int rc = MPI_SUCCESS;

  double Timeout_Begin = MPI_Wtime();

  char Alive_Status_Recv;
  char Alive_Status_Send;
  Alive_Status_Send = I_AM_ALIVE;

  MPI_Request Send_Request;
  if (source == MPI_ANY_SOURCE) {
    /*
     * if source is MPI_ANY_SOURCE, we cannot send an alive-message,
     * because we don't know the source rank, yet. So we wait for an
     * incoming alive-message and respond it afterwards. If the comm
     * is revoked or any process is dead, return with an error.
     */

    int RFlag = 0;
    MPI_Status RStatus;

    while (1) {
      MPI_Iprobe(source, tag + TAGOFFSETALIVE, f_comm->c_comm_copy_p2p, &RFlag, &RStatus);
      if (f_comm->dead_set.size() > 0 && f_comm->Ack_failed_processes == false) {
        return MPI_ERR_PROC_FAILED_PENDING;
        /* ######## IMPORTANT ########
         * note: in this case, another process possibly already sent an alive-message
         * and waits for a response. It is recommended to immediately revoke the communicator
         * in order to notify a possible sending process.
         * ###########################
         */
      } else if (RFlag == 1) {
        // receive the status message and respond with an alive-message
        MPI_Recv(&Alive_Status_Recv, 1, MPI_CHAR, RStatus.MPI_SOURCE, RStatus.MPI_TAG,
                 f_comm->c_comm_copy_p2p, MPI_STATUS_IGNORE);
        MPI_Isend(&Alive_Status_Send, 1, MPI_CHAR, RStatus.MPI_SOURCE, tag + TAGOFFSETALIVE,
                  f_comm->c_comm_copy_p2p, &Send_Request);
        MPI_Request_free(&Send_Request);
        break;
      } else if (f_comm->comm_revoked) {
        return MPI_ERR_REVOKED;
      }
      simft::Sim_FT_Perform_background_operations();
    }
  } else {
    MPI_Isend(&Alive_Status_Send, 1, MPI_CHAR, source, tag + TAGOFFSETALIVE,
              f_comm->c_comm_copy_p2p, &Send_Request);  // send alive-message
    MPI_Request_free(&Send_Request);

    int RecFlag = 0;

    // wait for a status message
    while (RecFlag == 0) {
      MPI_Iprobe(source, tag + TAGOFFSETALIVE, f_comm->c_comm_copy_p2p, &RecFlag,
                 MPI_STATUS_IGNORE);
      simft::Sim_FT_Perform_background_operations();
    }
    MPI_Recv(&Alive_Status_Recv, 1, MPI_CHAR, source, tag + TAGOFFSETALIVE, f_comm->c_comm_copy_p2p,
             MPI_STATUS_IGNORE);
  }

  // check the content of the status message and react to it
  if (Alive_Status_Recv == I_AM_ALIVE || Alive_Status_Recv == I_AM_ALIVE_NB) {
    if (status != SIM_FT_MPI_STATUS_IGNORE) {  // if the status object does not exist, we aren't
                                               // allowed to access it(!!!)
      rc = MPI_Recv(buf, count, type, source, tag, f_comm->c_comm, &status->c_status);
    } else {
      rc = MPI_Recv(buf, count, type, source, tag, f_comm->c_comm, MPI_STATUS_IGNORE);
    }

  } else if (Alive_Status_Recv == COMM_IS_REVOKED) {
    rc = MPI_ERR_REVOKED;

  } else {
    f_comm->dead_set.insert(source);
    while (MPI_Wtime() - Timeout_Begin < SIM_FT_TIMEOUT) {
      // simulate timeout
      simft::Sim_FT_Perform_background_operations();
    }
    rc = MPI_ERR_PROC_FAILED;
  }

  return rc;
}
