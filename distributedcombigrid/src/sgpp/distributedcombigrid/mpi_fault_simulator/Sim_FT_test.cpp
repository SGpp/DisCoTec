#include "MPI-FT.h"
#include REAL_MPI_INCLUDE
#include <iostream>

int simft::Sim_FT_MPI_Test(simft::Sim_FT_MPI_Request *request, int *flag,
                           simft::Sim_FT_MPI_Status *status) {
  simft::Sim_FT_decide_kill();
  simft::Sim_FT_Perform_background_operations();

  *flag = 0;
  int Ret = 0;

  if ((*request)->Request_Type == 1) {  // P2P

    if ((*request)->f_comm->comm_revoked) {
      if ((*request)->c_request != MPI_REQUEST_NULL) simft::Sim_FT_MPI_Request_free(request);
      *flag = 1;
      if (status != SIM_FT_MPI_STATUS_IGNORE) status->MPI_ERROR = MPI_ERR_REVOKED;
      return MPI_ERR_REVOKED;
    } else if ((*request)->P2P_Stats.p2p_rank == MPI_ANY_SOURCE &&
               (*request)->f_comm->dead_set.size() > 0 &&
               (*request)->f_comm->Ack_failed_processes == false) {
      // ULFM standard behavior for MPI_ANY_SOURCE. it is recommended to revoke the comm after this
      // error occurs!!!
      if ((*request)->c_request != MPI_REQUEST_NULL) simft::Sim_FT_MPI_Request_free(request);
      *flag = 1;
      if (status != SIM_FT_MPI_STATUS_IGNORE) status->MPI_ERROR = MPI_ERR_PROC_FAILED_PENDING;
      return MPI_ERR_PROC_FAILED_PENDING;

    } else if ((*request)->f_comm->dead_set.count((*request)->P2P_Stats.p2p_rank) > 0) {
      // TODO special case?
      if ((*request)->c_request != MPI_REQUEST_NULL) simft::Sim_FT_MPI_Request_free(request);

      *flag = 1;
      if (status != SIM_FT_MPI_STATUS_IGNORE) status->MPI_ERROR = MPI_ERR_PROC_FAILED;
      return MPI_ERR_PROC_FAILED;
    }

    if (!(*request)->completed) {
      int CFlag = 0;
      if (status != SIM_FT_MPI_STATUS_IGNORE) {
        Ret = MPI_Test(&(*request)->c_request, &CFlag, &status->c_status);
      } else {
        Ret = MPI_Test(&(*request)->c_request, &CFlag, MPI_STATUS_IGNORE);
      }
      if (status != SIM_FT_MPI_STATUS_IGNORE) {
        status->MPI_ERROR = status->c_status.MPI_ERROR;
        status->MPI_SOURCE = status->c_status.MPI_SOURCE;
        status->MPI_TAG = status->c_status.MPI_TAG;
        // status->cancelled = status->c_status.cancelled; //this only exists in MPICH, not in Open
        // MPI => better not assignt it
      }
      if (CFlag) {
        (*request)->completed = true;
        (*request)->completition_err_code = Ret;
      }

    } else if ((*request)->completed) {  // operation is completed, continue to check for
                                         // alive/status messages
      int Alive_Flag = 0;
      MPI_Status Alive_Status;
      MPI_Iprobe((*request)->P2P_Stats.p2p_rank, (*request)->P2P_Stats.tag,
                 (*request)->f_comm->c_comm_copy_p2p, &Alive_Flag, &Alive_Status);
      if (Alive_Flag) {
        *flag = 1;
        (*request)->completed = false;
        MPI_Status recv_status;
        char Alive_Status_Recv;
        MPI_Recv(&Alive_Status_Recv, 1, MPI_CHAR, Alive_Status.MPI_SOURCE,
                 (*request)->P2P_Stats.tag, (*request)->f_comm->c_comm_copy_p2p, &recv_status);

        // If source is MPI_ANY_SOURCE, alive-message wasn't sent, yet. Do it now!
        if ((*request)->P2P_Stats.p2p_rank == MPI_ANY_SOURCE) {
          MPI_Request Send_Request;
          char Alive_Status_Send = I_AM_ALIVE_NB;
          MPI_Isend(&Alive_Status_Send, 1, MPI_CHAR, Alive_Status.MPI_SOURCE,
                    (*request)->P2P_Stats.tag, (*request)->f_comm->c_comm_copy_p2p, &Send_Request);
          MPI_Request_free(&Send_Request);
        }
        if (Alive_Status_Recv == I_AM_ALIVE || Alive_Status_Recv == I_AM_ALIVE_NB) {
          if (status != SIM_FT_MPI_STATUS_IGNORE)
            status->MPI_ERROR = (*request)->completition_err_code;
          Ret = (*request)->completition_err_code;
        } else if (Alive_Status_Recv == COMM_IS_REVOKED) {
          if (status != SIM_FT_MPI_STATUS_IGNORE) status->MPI_ERROR = MPI_ERR_REVOKED;
          Ret = MPI_ERR_REVOKED;
        } else {
          //			    	while(MPI_Wtime() - Timeout_Begin < SIM_FT_TIMEOUT){
          //			    		//simulate timeout
          //			    	}
          (*request)->f_comm->dead_set.insert(Alive_Status.MPI_SOURCE);
          if (status != SIM_FT_MPI_STATUS_IGNORE) status->MPI_ERROR = MPI_ERR_PROC_FAILED;
          Ret = MPI_ERR_PROC_FAILED;
        }
        delete (*request);
        *request = &simft::Sim_FT_MPI_REQUEST_NULL;
      }
    }
    return Ret;

  } else if ((*request)->Request_Type == 2) {  // collective

    if ((*request)->f_comm->dead_set.size() > 0) {
      if (status != SIM_FT_MPI_STATUS_IGNORE) status->MPI_ERROR = MPI_ERR_PROC_FAILED;
      return MPI_ERR_PROC_FAILED;
      //		}
    } else {
      // TODO: delete request ? it was created by new...
      int ret;
      if (status != SIM_FT_MPI_STATUSES_IGNORE) {
        ret = MPI_Test(&(*request)->c_request, flag, &status->c_status);
      } else {
        ret = MPI_Test(&(*request)->c_request, flag, MPI_STATUS_IGNORE);
      }
      if (status != SIM_FT_MPI_STATUS_IGNORE) {
        status->MPI_ERROR = status->c_status.MPI_ERROR;
        status->MPI_SOURCE = status->c_status.MPI_SOURCE;
        status->MPI_TAG = status->c_status.MPI_TAG;
        // status->cancelled = status->c_status.cancelled; //this only exists in MPICH, not in Open
        // MPI => better not assignt it
      }
      return ret;
    }
  }

  return -1;
}
