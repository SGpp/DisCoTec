#include "../../include/discotec/MPI-FT.h"
#include REAL_MPI_INCLUDE
#include <iostream>

int simft::Sim_FT_MPI_Comm_free(simft::Sim_FT_MPI_Comm *f_comm) {
  // MPI_Comm_free is a collective function according to the standard, so we use our dead processes
  // functions.

  // first sync to make sure that in case of a revoked comm all processes switched to the revoked
  // tag
  simft::Sim_FT_Check_dead_processes(*f_comm, simft::NextOp::Default);

  if ((*f_comm)->comm_rank == (*f_comm)->Root_Rank) {
    if ((*f_comm)->CurrentMaxNBC_Rank != (*f_comm)->comm_rank) {
      /* ####################
       * request and receive NBC vector
       */
      MPI_Request NBC_Request;
      MPI_Isend(0, 0, MPI_INT, (*f_comm)->CurrentMaxNBC_Rank, SIM_FT_NBC_REQUEST_TAG,
                (*f_comm)->c_comm_copy, &NBC_Request);
      MPI_Request_free(&NBC_Request);

      MPI_Status RStatus;
      int PFlag = 0;
      while (PFlag == 0) {
        MPI_Iprobe((*f_comm)->CurrentMaxNBC_Rank, SIM_FT_NBC_RESPONSE_TAG, (*f_comm)->c_comm_copy,
                   &PFlag, &RStatus);
        simft::Sim_FT_Perform_background_operations();
      }
      int Count = 0;
      MPI_Get_count(&RStatus, MPI_INT, &Count);
      (*f_comm)->NBC_Vector_Send.resize(Count);

      MPI_Recv(&(*f_comm)->NBC_Vector_Send.front(), Count, MPI_INT, (*f_comm)->CurrentMaxNBC_Rank,
               SIM_FT_NBC_RESPONSE_TAG, (*f_comm)->c_comm_copy, MPI_STATUS_IGNORE);
      // ####################
    } else {
      // CurrentMaxNBC_Rank is root => serialize the local recent_icollectives
      (*f_comm)->NBC_Vector_Send = simft::NBC_to_Vector((*f_comm));
    }
  }
  for (int i = 0; i < static_cast<int>((*f_comm)->NBC_Vector_Send.size()); i++) {
    std::cout << "NBCOUT:::" << (*f_comm)->NBC_Vector_Send[i] << "\n";
  }

  // now send NBC vector to all processes

  simft::Sim_FT_Check_dead_processes(*f_comm, simft::NextOp::Comm_free);

  std::vector<simft::Sim_FT_Istats> Deserialized_NBC;
  Deserialized_NBC = simft::Vector_to_NBC((*f_comm)->NBC_Vector_Send);
  simft::Sim_FT_Perform_nb_operations(&Deserialized_NBC, (*f_comm));

  simft::Sim_FT_Complete_nb_operations(&Deserialized_NBC, (*f_comm));

  int Ret = MPI_Comm_free(&(*f_comm)->c_comm);

  if (0 == Ret) {  // if Ret != 0, we have a problem!
    MPI_Comm_free(&(*f_comm)->c_comm_copy);
    MPI_Comm_free(&(*f_comm)->c_comm_copy_p2p);
    MPI_Comm_free(&(*f_comm)->c_comm_copy_coll);
    MPI_Comm_free(&(*f_comm)->c_comm_copy_probe);

    simft::Sim_FT_Current_Active_Communicators.erase(
        *f_comm);    // remove freed communicator from active comm set
    delete *f_comm;  // delete custom comm object
    *f_comm = nullptr;
  }

  return Ret;
}

int simft::Sim_FT_MPI_Comm_free2(simft::Sim_FT_MPI_Comm *f_comm) {  // does not free c_comm
  // MPI_Comm_free is a collective function according to the standard, so we use our dead processes
  // functions.

  // first sync to make sure that in case of a revoked comm all processes switched to the revoked
  // tag
  simft::Sim_FT_Check_dead_processes(*f_comm, simft::NextOp::Default);

  if ((*f_comm)->comm_rank == (*f_comm)->Root_Rank) {
    if ((*f_comm)->CurrentMaxNBC_Rank != (*f_comm)->comm_rank) {
      /* ####################
       * request and receive NBC vector
       */
      MPI_Request NBC_Request;
      MPI_Isend(0, 0, MPI_INT, (*f_comm)->CurrentMaxNBC_Rank, SIM_FT_NBC_REQUEST_TAG,
                (*f_comm)->c_comm_copy, &NBC_Request);
      MPI_Request_free(&NBC_Request);

      MPI_Status RStatus;
      int PFlag = 0;
      while (PFlag == 0) {
        MPI_Iprobe((*f_comm)->CurrentMaxNBC_Rank, SIM_FT_NBC_RESPONSE_TAG, (*f_comm)->c_comm_copy,
                   &PFlag, &RStatus);
        simft::Sim_FT_Perform_background_operations();
      }
      int Count = 0;
      MPI_Get_count(&RStatus, MPI_INT, &Count);
      (*f_comm)->NBC_Vector_Send.resize(Count);

      MPI_Recv(&(*f_comm)->NBC_Vector_Send.front(), Count, MPI_INT, (*f_comm)->CurrentMaxNBC_Rank,
               SIM_FT_NBC_RESPONSE_TAG, (*f_comm)->c_comm_copy, MPI_STATUS_IGNORE);
      // ####################
    } else {
      // CurrentMaxNBC_Rank is root => serialize the local recent_icollectives
      (*f_comm)->NBC_Vector_Send = simft::NBC_to_Vector((*f_comm));
    }
  }
  for (int i = 0; i < static_cast<int>((*f_comm)->NBC_Vector_Send.size()); i++) {
    std::cout << "NBCOUT:::" << (*f_comm)->NBC_Vector_Send[i] << "\n";
  }

  // now send NBC vector to all processes

  simft::Sim_FT_Check_dead_processes(*f_comm, simft::NextOp::Comm_free);

  std::vector<simft::Sim_FT_Istats> Deserialized_NBC;
  Deserialized_NBC = simft::Vector_to_NBC((*f_comm)->NBC_Vector_Send);
  simft::Sim_FT_Perform_nb_operations(&Deserialized_NBC, (*f_comm));

  simft::Sim_FT_Complete_nb_operations(&Deserialized_NBC, (*f_comm));

  // int Ret = MPI_Comm_free(&(*f_comm)->c_comm);
  int Ret = 0;
  if (0 == Ret) {  // if Ret != 0, we have a problem!
    MPI_Comm_free(&(*f_comm)->c_comm_copy);
    MPI_Comm_free(&(*f_comm)->c_comm_copy_p2p);
    MPI_Comm_free(&(*f_comm)->c_comm_copy_coll);
    MPI_Comm_free(&(*f_comm)->c_comm_copy_probe);

    simft::Sim_FT_Current_Active_Communicators.erase(
        *f_comm);    // remove freed communicator from active comm set
    delete *f_comm;  // delete custom comm object
    *f_comm = nullptr;
  }

  return Ret;
}
