#include "discotec/MPI-FT.h"
#include REAL_MPI_INCLUDE
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>

#include <stdlib.h>

int RandKill_cnt = 0;

// check if current process is (simulated) root in any communicator. Implemented but currently not
// used => just return
bool simft::Sim_FT_proc_is_root() {
  return false;
  int rank = -1;
  for (auto Actives : simft::Sim_FT_Current_Active_Communicators) {
    simft::Sim_FT_MPI_Comm_rank(Actives, &rank);
    if (rank == Actives->Root_Rank) {
      return true;
    }
  }
  return false;
}

/* Step1: wait for incoming reduce msg and forward it
 * Step2: wait for incoming bcast msg and forward it
 *
 * This is a collective function. All (alive) ranks of the communicator have to call it in order to
 * continue.
 * The call only gets canceled if NextOp is Deafult and the communicator is revoked! For calls like
 * Shrink or
 * Finalize, the function MUST be called by all alive processes, the application needs to make sure
 * of it,
 * otherwise some processes will be stuck => possible deadlock.
 * */
int simft::Sim_FT_Check_dead_processes(simft::Sim_FT_MPI_Comm f_comm, simft::NextOp option) {
  // if comm is revoked and NextOp is Default, return error (no sync operation is needed in this
  // special case)
  if (simft::Sim_FT_Check_comm_revoked(f_comm) && option == simft::Default) {
    simft::Sim_FT_Perform_background_operations(f_comm, true);
    return MPI_ERR_REVOKED;
  }

  bool revokedAck;

  // If comm is revoked, the sync messages have to be sent via a special revoke tag. Only messages
  // other than Default can be sent this way.
  if (f_comm->comm_revoked && option != simft::Default) {
    revokedAck = true;
  } else {
    revokedAck = false;
  }

  int Reduce_Vector[4];

  int Ret;

  /*
   * Step1: reduce count of dead processes to root
   */
  // std::cout << " size of comm " << f_comm->comm_size << "\n";
  simft::Sim_FT_Custom_Reduce(Reduce_Vector, f_comm, option, revokedAck);

  if (f_comm->comm_rank == f_comm->Root_Rank) {
    f_comm->CurrentMaxNBC = Reduce_Vector[1];
    f_comm->CurrentMaxNBC_Rank = Reduce_Vector[2];
    // std::cout << "Number of dead processes is " << f_comm->Dead_Processes_Root.size() << "\n";
    // std::cout << "Reduce vector value is " << Reduce_Vector[0] << "\n";
    // receive dead messages until the reduced dead count is reached
    while ((int)f_comm->Dead_Processes_Root.size() < Reduce_Vector[0]) {
      simft::Sim_FT_Receive_dead_msgs_root(f_comm);
      simft::Sim_FT_Perform_background_operations(true);
      // std::cout << "Received dead msg from " << f_comm->Dead_Processes_Root[0] << "\n";
    }
    std::vector<int>
        Send_Vector;  //[0]: NextOp; [1]: Dead processes count; [>1]: Dead processes ids

    int addReserveSize = 2;
    if (option == simft::NextOp::Comm_free) {
      addReserveSize +=
         static_cast<int>(f_comm->NBC_Vector_Send.size());  // NBC_Vector_Send will be attached, reserve its size
    }

    Send_Vector.reserve(f_comm->Dead_Processes_Root.size() + addReserveSize);
    // Attach the information, which operation is next (Shrink, Finalize or normal blocking
    // collective)
    Send_Vector.push_back(option);
    //
    Send_Vector.push_back(static_cast<int>(f_comm->Dead_Processes_Root.size()));

    if (f_comm->Dead_Processes_Root.size() > 0) {
      // std::cout << "option = " <<option << "\n";
      Ret = MPI_ERR_PROC_FAILED;

      // attach dead processes
      Send_Vector.insert(Send_Vector.end(), f_comm->Dead_Processes_Root.begin(),
                         f_comm->Dead_Processes_Root.end());

      // if next operation is MPI_Comm_free, attach NBC_Vector_Send
      if (option == simft::NextOp::Comm_free)
        Send_Vector.insert(Send_Vector.end(), f_comm->NBC_Vector_Send.begin(),
                           f_comm->NBC_Vector_Send.end());

      simft::Sim_FT_Custom_Dead_bcast(&Send_Vector, f_comm, revokedAck);
      if (f_comm->Dead_Processes_Root.size() != f_comm->dead_nodes.size()) {
        f_comm->dead_nodes = f_comm->Dead_Processes_Root;
        f_comm->Ack_failed_processes = false;  // new failure detected => reset ack
      }
      // f_comm->Blocking_Collective_Count++; //TODO: not needed, there will be only one cycle per
      // comm (err after cycle complete)
      return Ret;
    } else {
      if (option == simft::NextOp::Shrink) {
        Ret = true;
      } else {
        Ret = false;
      }

      // if next operation is MPI_Comm_free, attach NBC_Vector_Send
      if (option == simft::NextOp::Comm_free)
        Send_Vector.insert(Send_Vector.end(), f_comm->NBC_Vector_Send.begin(),
                           f_comm->NBC_Vector_Send.end());

      simft::Sim_FT_Custom_Dead_bcast(&Send_Vector, f_comm,
                                      revokedAck);  // empty vector (no dead processes)
      // f_comm->Recent_Icollectives.clear(); //after a blocking collective, clear the vector at all
      // ranks
      // f_comm->Recent_Icollectives.reserve(10);

      if (f_comm->dead_nodes.size() != f_comm->dead_set.size()) {
        f_comm->dead_set = std::set<int>(
            f_comm->dead_nodes.begin(),
            f_comm->dead_nodes
                .end());  //"quite" inefficient (O(n logn)) but only occurs if processes failed
      }

      return Ret;
    }

  } else {  // process is not root
    std::vector<int> Bcast_Vector;
    simft::Sim_FT_Custom_Dead_bcast(&Bcast_Vector, f_comm, revokedAck);

    if (Bcast_Vector.size() > 0) {
      if (Bcast_Vector[1] > 0) {
        // at least one process in our comm is dead
        Ret = MPI_ERR_PROC_FAILED;
        if ((unsigned int)Bcast_Vector[1] != f_comm->dead_nodes.size()) {
          f_comm->dead_nodes =
              std::vector<int>(&Bcast_Vector[2], &Bcast_Vector[Bcast_Vector[1] + 2]);  //+1
          f_comm->Ack_failed_processes = false;  // new failure detected => reset ack
        }
        // f_comm->Blocking_Collective_Count++; //TODO: not needed, there will be only one cycle per
        // comm (err after cycle complete)
      } else {
        Ret = false;
      }

      if (Bcast_Vector[0] == 9) {  // comm is revoked, but we started this op while it wasn't
                                   // revoked => set revoked now and return error
        f_comm->comm_revoked = true;
        Ret = MPI_ERR_REVOKED;
      }

      // f_comm->NBC_Vector_Send = std::vector<int>(Bcast_Vector[1]+1,Bcast_Vector.size());
      if (Bcast_Vector.size() > (unsigned int)Bcast_Vector[1] + 2) {
        f_comm->NBC_Vector_Send = std::vector<int>(&Bcast_Vector[Bcast_Vector[1] + 2],
                                                   &Bcast_Vector[Bcast_Vector.size()]);
      }

      if (f_comm->dead_nodes.size() != f_comm->dead_set.size()) {
        f_comm->dead_set = std::set<int>(
            f_comm->dead_nodes.begin(),
            f_comm->dead_nodes
                .end());  //"quite" inefficient (O(n logn)) but only occurs if processes failed
      }
      return Ret;
    } else {
      // this should never happen
      std::cout << "FATAL ERROR in fault layer (check dead processes; non-root)\n";
      exit(1);
    }
  }
}

bool simft::Sim_FT_Receive_dead_msgs_root(simft::Sim_FT_MPI_Comm f_comm) {
  int Dead_Flag = 0;
  MPI_Status Dead_Status;
  MPI_Iprobe(MPI_ANY_SOURCE, SIM_FT_DEAD_NBC_TAG, f_comm->c_comm_copy, &Dead_Flag, &Dead_Status);
  if (Dead_Flag == 1) {
    // new dead message. Receive it.
    int NBC = 0;
    MPI_Recv(&NBC, 1, MPI_INT, Dead_Status.MPI_SOURCE, SIM_FT_DEAD_NBC_TAG, f_comm->c_comm_copy,
             MPI_STATUS_IGNORE);

    // add dead process to dead vector and dead set
    f_comm->Dead_Processes_Root.push_back(Dead_Status.MPI_SOURCE);
    f_comm->dead_set.insert(Dead_Status.MPI_SOURCE);

    // add dead process to recently detected dead processes
    f_comm->root_recent_dead_set.insert(Dead_Status.MPI_SOURCE);

    // if NBC count is lower than CurrentMinNBC, update it
    if (f_comm->CurrentMinNBC > NBC || f_comm->CurrentMinNBC == -1) {
      f_comm->CurrentMinNBC = NBC;
      f_comm->CurrentMinNBC_Rank = Dead_Status.MPI_SOURCE;
    }
    f_comm->CurrentMinNBC_Modified = true;  // tells root to initiate a NBC broadcast (to inform
                                            // comm about new dead processes and NBC count)
    return true;
  } else {
    return false;
  }
}

/* Very important function. Needs to be called frequently
 * so that background messages such as revoke-msg or NBC-msg can be forwarded through the whole
 * communicator.
 * If current rank is root, the function also receives dead-msgs and initiates NBC-broadcst from
 * time to time.
*/
void simft::Sim_FT_Perform_background_operations(simft::Sim_FT_MPI_Comm f_comm,
                                                 bool inside_Check_dead) {
  if (f_comm->Root_Rank == f_comm->comm_rank) {
    // receive all currently incoming dead messages from newly dead processes
    while (simft::Sim_FT_Receive_dead_msgs_root(f_comm))
      ;
    //
    simft::Sim_FT_NBC_bcast_root(f_comm);
  } else {
    // forward possible background dead broadcast
    simft::Sim_FT_NBC_bcast_Nroot(f_comm);

    /* ####################################
     * check for NBC requests from root.
     * if request is received, respond it
     * with the local list of recent icollectives */
    int PFlag = 0;
    MPI_Iprobe(f_comm->Root_Rank, SIM_FT_NBC_REQUEST_TAG, f_comm->c_comm_copy, &PFlag,
               MPI_STATUS_IGNORE);
    if (PFlag == 1) {
      MPI_Recv(0, 0, MPI_INT, f_comm->Root_Rank, SIM_FT_NBC_REQUEST_TAG, f_comm->c_comm_copy,
               MPI_STATUS_IGNORE);

      // root requests the current NBC vector. Send it now.
      MPI_Request NBC_Request;

      std::vector<int> NBCVec;
      // serialize recent NBC operations
      NBCVec = simft::NBC_to_Vector(f_comm);

      // send serialized list to root
      MPI_Isend(&NBCVec[0], static_cast<int>(NBCVec.size()), MPI_INT, f_comm->Root_Rank, SIM_FT_NBC_RESPONSE_TAG,
                f_comm->c_comm_copy, &NBC_Request);
      MPI_Request_free(&NBC_Request);
    }
    /* #################################### */
  }

  /*
   * this check is necessary for consistency. Otherwise a revoke broadcast could overtake a
   * Check_dead_processes call
   * resulting in some processes using the normal tags and some the special revoked tags
  */
  if (!inside_Check_dead) {
    // forward revoke broadcast if one is active (non-root) or initiate if revoke message is
    // received (root)
    simft::Sim_FT_Check_comm_revoked(f_comm);
  }

  if (f_comm->comm_revoked) {
    // Check for incoming alive-messages. Answer them if comm is revoked to free p2p communication
    int Alive_Flag = -1;
    MPI_Status Alive_Status;
    char Alive_Status_Recv = 0;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, f_comm->c_comm_copy_p2p, &Alive_Flag, &Alive_Status);
    if (Alive_Flag == 1) {
      MPI_Recv(&Alive_Status_Recv, 1, MPI_CHAR, Alive_Status.MPI_SOURCE, Alive_Status.MPI_TAG,
               f_comm->c_comm_copy_p2p, MPI_STATUS_IGNORE);
      if (Alive_Status_Recv == I_AM_ALIVE ||
          Alive_Status_Recv ==
              I_AM_ALIVE_NB) {  // should always be true!!! (TODO: maybe remove if???)
        MPI_Request Send_Request;
        char Revoked_Char = COMM_IS_REVOKED;
        MPI_Isend(&Revoked_Char, 1, MPI_CHAR, Alive_Status.MPI_SOURCE, Alive_Status.MPI_TAG,
                  f_comm->c_comm_copy_p2p, &Send_Request);
        MPI_Request_free(&Send_Request);
      }
      std::cout << "sending revoked char to processor: " << Alive_Status.MPI_SOURCE
                << " size: " << f_comm->comm_size << "\n";
    }

    // check for incoming reduce from Check_dead_processes on the non-revoked tag and respond with a
    // bcast-message with NextOp "comm is revoked"
    int Recv_Flag = 0;
    MPI_Status Recv_Status;
    MPI_Iprobe(MPI_ANY_SOURCE, SIM_FT_REDUCE_TAG, f_comm->c_comm_copy_coll, &Recv_Flag,
               &Recv_Status);

    if (Recv_Flag == 1) {
      std::cout << "Answering with revoke comm to rank: " << Recv_Status.MPI_SOURCE
                << " from rank: " << f_comm->comm_rank << " size " << f_comm->comm_size << "\n";
      std::vector<int> Receive_Vector;
      int count;
      MPI_Get_count(&Recv_Status, MPI_INT, &count);
      Receive_Vector.resize(count);
      MPI_Recv(&Receive_Vector[0], count, MPI_INT, Recv_Status.MPI_SOURCE, SIM_FT_REDUCE_TAG,
               f_comm->c_comm_copy_coll, MPI_STATUS_IGNORE);

      // Respond with NextOp RevokedComm
      MPI_Request Send_Request;
      std::vector<int> Send_Vector;
      Send_Vector.resize(3);
      Send_Vector[0] = simft::RevokedComm;
      MPI_Isend(&Send_Vector[0], 1, MPI_INT, Recv_Status.MPI_SOURCE, SIM_FT_DEAD_TAG,
                f_comm->c_comm_copy_coll, &Send_Request);
      MPI_Request_free(&Send_Request);
    }
  }
}

/* Special: perform background operations on all active communicators (possibly slowing down if too
 * many active communicators)
 * */
void simft::Sim_FT_Perform_background_operations(bool inside_Ckeck_dead) {
  // for(int i = 0; i < (int)simft::Sim_FT_Current_Active_Communicators.size(); i++){
  for (auto Actives : simft::Sim_FT_Current_Active_Communicators) {
    simft::Sim_FT_Perform_background_operations(Actives, inside_Ckeck_dead);
  }
}

/* Same function - but for dead processes */
void simft::Sim_FT_Perform_background_operations_dead(simft::Sim_FT_MPI_Comm f_comm) {
  if (f_comm->Root_Rank == f_comm->comm_rank) {
    // std::cout << "ROOT check revoked\n";
    simft::Sim_FT_Check_comm_revoked(f_comm);
    while (simft::Sim_FT_Receive_dead_msgs_root(f_comm))
      ;
    simft::Sim_FT_NBC_bcast_root(f_comm);
    if ((int)simft::Sim_FT_MPI_COMM_WORLD->Dead_Processes_Root.size() ==
        simft::Sim_FT_MPI_COMM_WORLD->comm_size) {
      /* ############################################
       * very special case:
       * if all processes are dead in MPI_COMM_WORLD,
       * initiate finalize protocol and exit
       * ############################################ */

      while (!simft::Sim_FT_Custom_Ireduce(simft::Sim_FT_MPI_COMM_WORLD, true))
        ;

      std::vector<int>
          Send_Vector;  //[0]: NextOp; [1]: Dead processes count; [>1]: Dead processes ids

      Send_Vector.reserve(simft::Sim_FT_MPI_COMM_WORLD->Dead_Processes_Root.size() + 2);
      // Next op is Finalize
      Send_Vector.push_back(2);
      //
      Send_Vector.push_back(static_cast<int>(simft::Sim_FT_MPI_COMM_WORLD->Dead_Processes_Root.size()));

      Send_Vector.insert(Send_Vector.end(),
                         simft::Sim_FT_MPI_COMM_WORLD->Dead_Processes_Root.begin(),
                         simft::Sim_FT_MPI_COMM_WORLD->Dead_Processes_Root.end());

      simft::Sim_FT_Custom_Ibcast_Root(&Send_Vector, &(simft::Sim_FT_MPI_COMM_WORLD->Dead_Buf),
                                       static_cast<int>(Send_Vector.size()), MPI_INT, SIM_FT_DEAD_TAG,
                                       simft::Sim_FT_MPI_COMM_WORLD->c_comm_copy_coll,
                                       simft::Sim_FT_MPI_COMM_WORLD);

      std::cout << ">> Sim_FT_MPI note: all processes are dead. Finalize and exit now! <<\n";
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(0);
    }
  } else {
    simft::Sim_FT_Check_comm_revoked(f_comm);
    simft::Sim_FT_NBC_bcast_Nroot(f_comm);
  }

  /* ####################################
   * check for NBC requests from root.
   * if request is received, respond it
   * with the local list of recent icollectives */
  int PFlag = 0;
  MPI_Iprobe(f_comm->Root_Rank, SIM_FT_NBC_REQUEST_TAG, f_comm->c_comm_copy, &PFlag,
             MPI_STATUS_IGNORE);
  if (PFlag == 1) {
    MPI_Recv(0, 0, MPI_INT, f_comm->Root_Rank, SIM_FT_NBC_REQUEST_TAG, f_comm->c_comm_copy,
             MPI_STATUS_IGNORE);

    // root requests the current NBC vector. Send it now.
    MPI_Request NBC_Request;

    std::vector<int> NBCVec;
    // serialize recent NBC operations
    NBCVec = simft::NBC_to_Vector(f_comm);
    // send serialized list to root
    MPI_Isend(&NBCVec[0], static_cast<int>(NBCVec.size()), MPI_INT, f_comm->Root_Rank, SIM_FT_NBC_RESPONSE_TAG,
              f_comm->c_comm_copy, &NBC_Request);
    MPI_Request_free(&NBC_Request);
  }
  /* #################################### */

  /* ########################
   * check for alive-requests */
  int Status_Flag = 0;
  MPI_Status Status_Status;
  MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, f_comm->c_comm_copy_p2p, &Status_Flag, &Status_Status);
  if (Status_Flag == 1) {
    char Rec_Char;
    MPI_Recv(&Rec_Char, 1, MPI_CHAR, Status_Status.MPI_SOURCE, Status_Status.MPI_TAG,
             f_comm->c_comm_copy_p2p, MPI_STATUS_IGNORE);
    if (Rec_Char == I_AM_ALIVE) {  // I_AM_ALIVE_NB will not be responded! Instead, the alive
                                   // process will know about the dead process after receiving nbc
                                   // bcast
      // Respond alive-message with dead-message
      char DeadMsg = I_AM_DEAD;
      MPI_Send(&DeadMsg, 1, MPI_CHAR, Status_Status.MPI_SOURCE, Status_Status.MPI_TAG,
               f_comm->c_comm_copy_p2p);
    }
  }
  /* ########################*/

  if (f_comm->comm_revoked) {
    // check for incoming reduce from Check_dead_processes on the non-revoked tag and respond with a
    // bcast-message with NextOp "comm is revoked"
    int Recv_Flag = 0;
    MPI_Status Recv_Status;
    MPI_Iprobe(MPI_ANY_SOURCE, SIM_FT_REDUCE_TAG, f_comm->c_comm_copy_coll, &Recv_Flag,
               &Recv_Status);

    if (Recv_Flag == 1) {
      std::vector<int> Receive_Vector;
      int count;
      MPI_Get_count(&Recv_Status, MPI_INT, &count);
      Receive_Vector.resize(count);
      MPI_Recv(&Receive_Vector[0], count, MPI_INT, Recv_Status.MPI_SOURCE, SIM_FT_REDUCE_TAG,
               f_comm->c_comm_copy_coll, MPI_STATUS_IGNORE);

      // Respond with NextOp RevokedComm
      MPI_Request Send_Request;
      std::vector<int> Send_Vector;
      Send_Vector.resize(3);
      Send_Vector[0] = simft::RevokedComm;
      MPI_Isend(&Send_Vector[0], 1, MPI_INT, Recv_Status.MPI_SOURCE, SIM_FT_DEAD_TAG,
                f_comm->c_comm_copy_coll, &Send_Request);
      MPI_Request_free(&Send_Request);
    }
  }
}

void simft::Sim_FT_NBC_bcast_root(simft::Sim_FT_MPI_Comm f_comm) {
  if (f_comm->CurrentMinNBC_Modified) {
    if (f_comm->LastNBC_bcast + SIM_FT_WAIT_NBC_BCAST < MPI_Wtime()) {
      int SendNBC = f_comm->CurrentMaxNBC;
      std::vector<int> Send_Vector;
      std::vector<int> Buffer_Vector;
      Send_Vector.reserve(f_comm->root_recent_dead_set.size() + 1);
      Send_Vector.push_back(SendNBC);
      // attach newly dead processes to the nbc broadcast
      Send_Vector.insert(Send_Vector.end(), f_comm->root_recent_dead_set.begin(),
                         f_comm->root_recent_dead_set.end());
      f_comm->root_recent_dead_set.clear();

      simft::Sim_FT_Custom_Ibcast_Root(&Send_Vector, &Buffer_Vector, static_cast<int>(Send_Vector.size()), MPI_INT,
                                       SIM_FT_BCAST_NBC_TAG, f_comm->c_comm_copy, f_comm);
      f_comm->CurrentMinNBC_Modified = false;
      f_comm->LastNBC_bcast = MPI_Wtime();
      std::cout << "Root start NBC bcast \n";
    }
  }
}

int simft::Sim_FT_NBC_bcast_Nroot(simft::Sim_FT_MPI_Comm f_comm) {
  int Count = 0;
  std::vector<int> NBC_Count;
  std::vector<int> Buffer_Vector;
  if (simft::Sim_FT_Custom_Ibcast_NRoot(&NBC_Count, &Buffer_Vector, &Count, MPI_INT,
                                        SIM_FT_BCAST_NBC_TAG, f_comm->c_comm_copy, f_comm)) {
    for (unsigned int i = 1; i < NBC_Count.size(); i++) {
      f_comm->dead_set.insert(NBC_Count[i]);
    }
    if (NBC_Count[0] < f_comm->CurrentMaxNBC || f_comm->CurrentMaxNBC == -1) {
      f_comm->CurrentMaxNBC = NBC_Count[0];
      return NBC_Count[0];
    } else {
      return f_comm->CurrentMaxNBC;
    }
  } else {
    return f_comm->CurrentMaxNBC;
  }
}

void simft::Sim_FT_Perform_nb_operations(std::vector<simft::Sim_FT_Istats> *NB_ops_to_perform,
                                         simft::Sim_FT_MPI_Comm f_comm) {
#ifndef DISABLE_NBC
  if (f_comm->Recent_Icollectives.size() <
      NB_ops_to_perform->size()) {  // else: no NB Collective to be done
    // std::cout << "ICOLLECTIVE WITH DUMMY DATA NEXT!!!\n";
    // std::vector<MPI_Request> Coll_Req;
    // Coll_Req.reserve(NB_ops_to_perform->size() - (f_comm->Recent_Icollectives.size() + 1));

    for (int i = static_cast<int>(f_comm->Recent_Icollectives.size()); i < (int)NB_ops_to_perform->size(); i++) {
      std::cout << "op type = " << (*NB_ops_to_perform)[i].Type << "\n";
      (*NB_ops_to_perform)[i].request = new simft::Sim_FT_MPI_Request_struct;

      switch ((*NB_ops_to_perform)[i].Type) {
        case 0:  // ibarrier
        {
          // std::cout << "Perform_NB ibarrier\n";
          MPI_Ibarrier(f_comm->c_comm, &(*NB_ops_to_perform)[i].request->c_request);
          break;
        }
        case 1:  // ibcast
        {
          if ((*NB_ops_to_perform)[i].Datatype == 0) {  // MPI_CHAR

            (*NB_ops_to_perform)[i].BuffChar.resize((*NB_ops_to_perform)[i].count);

            // begin Ibcast with dummy Array to free other nodes with active Ibcasts
            MPI_Ibcast(&(*NB_ops_to_perform)[i].BuffChar.front(), (*NB_ops_to_perform)[i].count,
                       MPI_CHAR, (*NB_ops_to_perform)[i].root, f_comm->c_comm,
                       &(*NB_ops_to_perform)[i].request->c_request);

          } else if ((*NB_ops_to_perform)[i].Datatype == 1) {  // MPI_INT

            (*NB_ops_to_perform)[i].BuffInt.resize((*NB_ops_to_perform)[i].count);

            // begin Ibcast with dummy Array to free other nodes with active Ibcasts
            MPI_Ibcast(&(*NB_ops_to_perform)[i].BuffInt.front(), (*NB_ops_to_perform)[i].count,
                       MPI_INT, (*NB_ops_to_perform)[i].root, f_comm->c_comm,
                       &(*NB_ops_to_perform)[i].request->c_request);

          } else if ((*NB_ops_to_perform)[i].Datatype == 2) {  // MPI_DOUBLE

            (*NB_ops_to_perform)[i].BuffDouble.resize((*NB_ops_to_perform)[i].count);

            // begin Ibcast with dummy Array to free other nodes with active Ibcasts
            MPI_Ibcast(&(*NB_ops_to_perform)[i].BuffDouble.front(), (*NB_ops_to_perform)[i].count,
                       MPI_DOUBLE, (*NB_ops_to_perform)[i].root, f_comm->c_comm,
                       &(*NB_ops_to_perform)[i].request->c_request);

          } else {
            /*###################################
             * Add other MPI Datatypes if needed!
             *###################################*/

            // Datatype not yet implemented, cancel now
            std::cout << "SIM_FT FATAL ERROR: In simft::Sim_FT_Perform_nb_operations, ibcast: "
                         "Datatype not yet implemented - "
                      << (*NB_ops_to_perform)[i].Datatype << "\n";
            exit(-1);
          }

          break;
        }
        case 2:  // iallreduce
        {
          if ((*NB_ops_to_perform)[i].Datatype == 0) {  // MPI_CHAR

            (*NB_ops_to_perform)[i].BuffChar.resize((*NB_ops_to_perform)[i].count);
            (*NB_ops_to_perform)[i].BuffChar2.resize((*NB_ops_to_perform)[i].count);

            // begin Iallreduce with dummy Array to free other nodes with active Iallreduces
            MPI_Iallreduce(&(*NB_ops_to_perform)[i].BuffChar.front(),
                           &(*NB_ops_to_perform)[i].BuffChar2.front(),
                           (*NB_ops_to_perform)[i].count, MPI_CHAR,
                           simft::Op_to_MOp((*NB_ops_to_perform)[i].op), f_comm->c_comm,
                           &(*NB_ops_to_perform)[i].request->c_request);

          } else if ((*NB_ops_to_perform)[i].Datatype == 1) {  // MPI_INT

            (*NB_ops_to_perform)[i].BuffInt.resize((*NB_ops_to_perform)[i].count);
            (*NB_ops_to_perform)[i].BuffInt2.resize((*NB_ops_to_perform)[i].count);

            // begin Iallreduce with dummy Array to free other nodes with active Iallreduces
            MPI_Iallreduce(&(*NB_ops_to_perform)[i].BuffInt.front(),
                           &(*NB_ops_to_perform)[i].BuffInt2.front(), (*NB_ops_to_perform)[i].count,
                           MPI_INT, simft::Op_to_MOp((*NB_ops_to_perform)[i].op), f_comm->c_comm,
                           &(*NB_ops_to_perform)[i].request->c_request);

          } else if ((*NB_ops_to_perform)[i].Datatype == 2) {  // MPI_DOUBLE

            (*NB_ops_to_perform)[i].BuffDouble.resize((*NB_ops_to_perform)[i].count);
            (*NB_ops_to_perform)[i].BuffDouble2.resize((*NB_ops_to_perform)[i].count);

            // begin Iallreduce with dummy Array to free other nodes with active Iallreduces
            MPI_Iallreduce(&(*NB_ops_to_perform)[i].BuffDouble.front(),
                           &(*NB_ops_to_perform)[i].BuffDouble2.front(),
                           (*NB_ops_to_perform)[i].count, MPI_DOUBLE,
                           simft::Op_to_MOp((*NB_ops_to_perform)[i].op), f_comm->c_comm,
                           &(*NB_ops_to_perform)[i].request->c_request);

          } else {
            /*###################################
             * Add other MPI Datatypes if needed!
             *###################################*/

            // Datatype not yet implemented, cancel now
            std::cout << simft::Sim_FT_MPI_COMM_WORLD->comm_rank;
            std::cout << "SIM_FT FATAL ERROR: In simft::Sim_FT_Perform_nb_operations, iallreduce: "
                         "Datatype not implemented yet - "
                      << (*NB_ops_to_perform)[i].Datatype << "\n";
            exit(1);
          }
          break;
        }

        default:
          std::cout << "SIM_FT FATAL ERROR: In simft::Sim_FT_Perform_nb_operations: Function Type "
                       "not implemented - "
                    << (*NB_ops_to_perform)[i].Type << "\n";
          exit(-1);
      }
    }
  }
#endif
}

void simft::Sim_FT_Complete_nb_operations(
    std::vector<simft::Sim_FT_Istats> *Deserialized_NBC_Vector, simft::Sim_FT_MPI_Comm f_comm) {
  std::vector<MPI_Request> Requests_To_Finish;
  Requests_To_Finish.reserve(f_comm->NBC_Vector_Send.size());

  // get all unfinished operations from Recent_icollectives and add them to Requests_To_Finish
  for (unsigned int i = 0; i < f_comm->Recent_Icollectives.size(); i++) {
    if (f_comm->Recent_Icollectives[i].request != nullptr) {
      if (f_comm->Recent_Icollectives[i].request->c_request != MPI_REQUEST_NULL) {
        Requests_To_Finish.push_back(f_comm->Recent_Icollectives[i].request->c_request);
        delete f_comm->Recent_Icollectives[i].request;  // created with new => need to delete
        f_comm->Recent_Icollectives[i].request = nullptr;
      }
    }
  }
  /* do the same for all other non-blocking collectives, received from synchronization
   * and NOT known to the process prior to the synchronization protocol
   */
  for (unsigned int i = static_cast<unsigned int>(f_comm->Recent_Icollectives.size()); i < (*Deserialized_NBC_Vector).size();
       i++) {
    if ((*Deserialized_NBC_Vector)[i].request != nullptr) {
      if ((*Deserialized_NBC_Vector)[i].request->c_request != MPI_REQUEST_NULL) {
        Requests_To_Finish.push_back((*Deserialized_NBC_Vector)[i].request->c_request);
        delete (*Deserialized_NBC_Vector)[i].request;  // created with new => need to delete
        (*Deserialized_NBC_Vector)[i].request = nullptr;
      }
    }
  }

  /* complete all requests. This can be called without having to call Perform_background_operations,
   * because all processes in the communicator (including the dead ones) called the same
   * non-blocking collective operations up to this point and so all requests will be completed.
   */
  MPI_Waitall(static_cast<int>(Requests_To_Finish.size()), &Requests_To_Finish.front(), MPI_STATUSES_IGNORE);

  for (unsigned int i = static_cast<unsigned int>(f_comm->Recent_Icollectives.size()); i < (*Deserialized_NBC_Vector).size();
       i++) {
    /* In case any of our buffers were allocated using new, delete them.
     * This is safe because we initialized them with nullptr. (maybe use smart pointers instead?)
     * */
    //		if((*Deserialized_NBC_Vector)[i].BuffChar != nullptr)delete[]
    //(*Deserialized_NBC_Vector)[i].BuffChar;
    //		if((*Deserialized_NBC_Vector)[i].BuffChar2 != nullptr)delete[]
    //(*Deserialized_NBC_Vector)[i].BuffChar2;
    //		if((*Deserialized_NBC_Vector)[i].BuffInt != nullptr)delete[]
    //(*Deserialized_NBC_Vector)[i].BuffInt;
    //		while(1);
    //		delete[] (*Deserialized_NBC_Vector)[i].BuffInt2;
    //		delete[] (*Deserialized_NBC_Vector)[i].BuffDouble;
    //		delete[] (*Deserialized_NBC_Vector)[i].BuffDouble2;

    //		(*Deserialized_NBC_Vector)[i].BuffChar = nullptr;
    //		(*Deserialized_NBC_Vector)[i].BuffChar2 = nullptr;
    //		(*Deserialized_NBC_Vector)[i].BuffInt = nullptr;
    //		(*Deserialized_NBC_Vector)[i].BuffInt2 = nullptr;
    //		(*Deserialized_NBC_Vector)[i].BuffDouble = nullptr;
    //		(*Deserialized_NBC_Vector)[i].BuffDouble2 = nullptr;
  }
}

bool simft::Sim_FT_Check_comm_revoked(simft::Sim_FT_MPI_Comm f_comm) {
  if (f_comm->comm_rank == f_comm->Root_Rank) {
    if (f_comm->comm_size == 1) {
      return f_comm->comm_revoked;
    }
    if (f_comm->comm_revoked == true) {
      // if another rank decides to revoke comm, receive message and "discard" it to free the
      // resource
      int Revoke_Flag = 0;
      MPI_Status Revoke_Status;
      MPI_Iprobe(MPI_ANY_SOURCE, SIM_FT_REVOKE_TAG, f_comm->c_comm_copy, &Revoke_Flag,
                 &Revoke_Status);
      if (Revoke_Flag == 1) {
        MPI_Recv(0, 0, MPI_CHAR, Revoke_Status.MPI_SOURCE, SIM_FT_REVOKE_TAG, f_comm->c_comm_copy,
                 MPI_STATUS_IGNORE);
      }
      return true;
    } else {
      int Revoke_Flag = 0;
      MPI_Status Revoke_Status;
      MPI_Iprobe(MPI_ANY_SOURCE, SIM_FT_REVOKE_TAG, f_comm->c_comm_copy, &Revoke_Flag,
                 &Revoke_Status);
      if (Revoke_Flag == 1) {
        // receive probed revoke message
        MPI_Recv(0, 0, MPI_CHAR, Revoke_Status.MPI_SOURCE, SIM_FT_REVOKE_TAG, f_comm->c_comm_copy,
                 MPI_STATUS_IGNORE);
        std::cout << "received revoke message from: " << Revoke_Status.MPI_SOURCE
                  << " I am rank: " << f_comm->comm_rank << " size: " << f_comm->comm_size << "\n";

        if (f_comm->comm_revoked == false) {
          // Comm not revoked yet, initialize background revoke broadcast
          int BCast_Revoke_Message = 2;
          std::vector<int> Send_Vector;
          std::vector<int> Buffer_Vector;
          Send_Vector.push_back(BCast_Revoke_Message);
          simft::Sim_FT_Custom_Ibcast_Root(&Send_Vector, &Buffer_Vector, 1, MPI_INT,
                                           SIM_FT_BCAST_REVOKE_TAG, f_comm->c_comm_copy, f_comm);

          f_comm->comm_revoked = true;
          return true;
        } else {
          // comm already revoked, no need to send revoke bcast again
          return true;
        }
        return true;
      } else {
        return false;
      }
    }
  } else {
    //		if(f_comm->comm_revoked == true){
    //			return true;
    //		}else{

    int Count;
    std::vector<int> BCast_Revoke_Message;
    std::vector<int> Buffer_Vector;

    // check for revoke broadcast
    bool BCast_Ret =
        simft::Sim_FT_Custom_Ibcast_NRoot(&BCast_Revoke_Message, &Buffer_Vector, &Count, MPI_INT,
                                          SIM_FT_BCAST_REVOKE_TAG, f_comm->c_comm_copy, f_comm);
    if (BCast_Ret) {
      // broadcast finished, set comm to revoked
      if (BCast_Revoke_Message[0] == 2) {
        std::cout << "revoking communicator at rank: " << f_comm->comm_rank
                  << " size: " << f_comm->comm_size << "\n";
        f_comm->comm_revoked = true;
        return true;
      } else {
        std::cout << "++####### ERROR RECV MS Msg =" << BCast_Revoke_Message[0] << "bei Rank "
                  << f_comm->comm_rank << "\n";
        return false;
      }
    } else {
      return f_comm->comm_revoked;
    }
    //		}
  }
}

void simft::Sim_FT_Send_revoke_message(simft::Sim_FT_MPI_Comm f_comm) {
  if (f_comm->comm_size == 1 || f_comm->comm_revoked == true) {
    // if comm size is 1 or comm already revoked, no need to send a revoke message
    f_comm->comm_revoked = true;
  } else {
    MPI_Request Revoke_Req;
    MPI_Isend(0, 0, MPI_CHAR, f_comm->Root_Rank, SIM_FT_REVOKE_TAG, f_comm->c_comm_copy,
              &Revoke_Req);
    // std::cout << "sending revoke message to root: " << f_comm->Root_Rank << " from rank: " <<
    // f_comm->comm_rank << "\n";
    MPI_Request_free(&Revoke_Req);

    // do not set f_comm->comm_revoked to true right now, instead wait to receive a revoke bcast and
    // only then set it to true!
    /*if(f_comm->comm_rank != f_comm->Root_Rank){
      f_comm->comm_revoked = true; //we do it because of dead locks
    }*/
  }
}

void simft::Sim_FT_Check_comm_root(simft::Sim_FT_MPI_Comm f_comm) {
  return;

  // for better performance: if necessary, this function can be extended to change a root of our
  // custom comm object, if the process is already root of another comm (not implemented yet!)
  /*
  int current_rank;
  MPI_Comm_rank(f_comm->c_comm, &current_rank);
  if(current_rank == f_comm->Root_Rank){
          if(simft::Sim_FT_proc_is_root()){
                  //TODO: set new root
                  //int Comm_Msg = SIM_FT_CHANGE_ROOT;
                  //MPI_Bcast(&Comm_Msg, 1, MPI_INT, f_comm->Root_Rank, f_comm->c_comm);
          }else{
                  //int Comm_Msg = SIM_FT_KEEP_ROOT;
                  //MPI_Bcast(&Comm_Msg, 1, MPI_INT, f_comm->Root_Rank, f_comm->c_comm);
          }
  }else{
          //int Comm_Msg;
          //MPI_Bcast(&Comm_Msg, 1, MPI_INT, f_comm->Root_Rank, f_comm->c_comm);
          //if(Comm_Msg == SIM_FT_KEEP_ROOT){
                  //root is ok, go on
          //}
  }
  */
}

// Because of the custom communicator, we have to initialize some variables
void simft::Sim_FT_Initialize_new_comm(simft::Sim_FT_MPI_Comm *f_new_comm, bool firstInit) {
  // if c_comm is MPI_COMM_NULL, we can't (and mustn't) call MPI_Comm_dup => set f_new_comm to NULL
  // and return
  if ((*f_new_comm)->c_comm == MPI_COMM_NULL) {
    delete (*f_new_comm);
    (*f_new_comm) = simft::Sim_FT_MPI_COMM_NULL;
    return;
  }

  // if root == -1 set to last process
  if ((*f_new_comm)->Root_Rank == -1) {
    int size;
    MPI_Comm_size((*f_new_comm)->c_comm, &size);
    std::cout << "Setting root rank to " << size - 1 << "\n";
    (*f_new_comm)->Root_Rank = size - 1;
  }

  if (!firstInit) {
    simft::Sim_FT_Check_comm_root((*f_new_comm));
  }

  MPI_Comm_dup((*f_new_comm)->c_comm, &(*f_new_comm)->c_comm_copy);
  MPI_Comm_dup((*f_new_comm)->c_comm, &(*f_new_comm)->c_comm_copy_coll);
  MPI_Comm_dup((*f_new_comm)->c_comm, &(*f_new_comm)->c_comm_copy_p2p);
  MPI_Comm_dup((*f_new_comm)->c_comm, &(*f_new_comm)->c_comm_copy_probe);

  MPI_Comm_size((*f_new_comm)->c_comm, &(*f_new_comm)->comm_size);
  MPI_Comm_rank((*f_new_comm)->c_comm, &(*f_new_comm)->comm_rank);

  // calculate successors and predecessor for our custom tree topology
  (*f_new_comm)->Bcast_Successors = simft::Bcast_Get_Successors(
      (*f_new_comm)->Root_Rank, (*f_new_comm)->comm_rank, (*f_new_comm)->comm_size);
  (*f_new_comm)->Bcast_Predecessor = simft::Bcast_Get_Predecessor(
      (*f_new_comm)->Root_Rank, (*f_new_comm)->comm_rank, (*f_new_comm)->comm_size);

  (*f_new_comm)->Recv_Request.resize((*f_new_comm)->Bcast_Successors.size());
  (*f_new_comm)->Received.resize((*f_new_comm)->Bcast_Successors.size());
  (*f_new_comm)->Send_Vector.resize(3);

  (*f_new_comm)->LastNBC_bcast = MPI_Wtime();

  // add the new communicator to set of active communicators
  simft::Sim_FT_Current_Active_Communicators.insert(*f_new_comm);
}

void simft::Sim_FT_kill_me() {
  // Send current Non-Blocking Collectives count to all root processes - at the same time notify
  // them about dead status
  // for(unsigned int i = 0; i < simft::Sim_FT_Current_Active_Communicators.size(); i++){
  for (auto Actives : simft::Sim_FT_Current_Active_Communicators) {
    MPI_Request Send_Request;
    int NBC = static_cast<int>(Actives->Recent_Icollectives.size());
    // Send dead-message to root
    MPI_Isend(&NBC, 1, MPI_INT, Actives->Root_Rank, SIM_FT_DEAD_NBC_TAG, Actives->c_comm_copy,
              &Send_Request);
    MPI_Request_free(&Send_Request);

    // Initialize values for next reduce
    Actives->Send_Vector.resize(4);
    Actives->Send_Vector[0] = 0;
    Actives->Send_Vector[1] = static_cast<int>(Actives->Recent_Icollectives.size());
    Actives->Send_Vector[2] = Actives->comm_rank;
    Actives->Send_Vector[3] = 0;

    if (Actives
            ->comm_revoked) {  // comm is revoked, use the alternative tags for reduce and broadcast
      Actives->ReduceTag = SIM_FT_REDUCE_TAG_REVOKED;
      Actives->DeadTag = SIM_FT_DEAD_TAG_REVOKED;
    }
  }

  while (1) {
    // loop through the set of active communicators
    for (auto Actives : simft::Sim_FT_Current_Active_Communicators) {
      if (Actives->comm_size ==
          1) {  // if comm_size is 1, no alive process is in there: free comm and delete it
        MPI_Comm_free(&Actives->c_comm_copy);
        MPI_Comm_free(&Actives->c_comm_copy_p2p);
        MPI_Comm_free(&Actives->c_comm_copy_coll);
        MPI_Comm_free(&Actives->c_comm_copy_probe);

        std::cout << "DEAD DELETE COMM OBJ\n";
        simft::Sim_FT_Current_Active_Communicators.erase(
            Actives);    // remove freed communicator from active comm set
        delete Actives;  // delete custom comm object
        Actives = nullptr;
        break;
      }
      // std::cout << "revoked: " << Actives->comm_revoked <<" size: " << Actives->comm_size << "
      // name: "  << Actives << "rank" << Actives->comm_rank << "\n";
      if (Actives->comm_revoked) {  // change tags if communicator is revoked
        Actives->ReduceTag = SIM_FT_REDUCE_TAG_REVOKED;
        Actives->DeadTag = SIM_FT_DEAD_TAG_REVOKED;
      }
      std::vector<int> buf;
      int count = 0;

      if (Actives->Root_Rank != Actives->comm_rank) {
        if (Actives->Ireduce_Active && simft::Sim_FT_Custom_Ibcast_NRoot(
                                           &buf, &(Actives->Dead_Buf), &count, MPI_INT,
                                           Actives->DeadTag, Actives->c_comm_copy_coll, Actives)) {
          // Ibcast completed, set Ireduce_Active to false, so that a new reduce operation can be
          // started
          Actives->Ireduce_Active = false;

          if (Actives->dead_nodes.size() != Actives->dead_set.size()) {
            Actives->dead_set = std::set<int>(
                Actives->dead_nodes.begin(),
                Actives->dead_nodes
                    .end());  // inefficient (O(n logn)) but only occurs if new processes failed
          }
          // std::cout <<"Command tag: " << buf[0] << "\n";
          if (buf[0] == 1) {  // Shrink next
            MPI_Comm EmptyComm;
            MPI_Comm_split(Actives->c_comm, MPI_UNDEFINED, 0, &EmptyComm);
            break;
          } else if (buf[0] == 2) {  // Finalize next
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Finalize();
            exit(0);
          } else if (buf[0] == 3) {  // comm_free next

            // extract the NBC vector from the last received bcast vector
            if (buf.size() > (unsigned int)buf[1] + 2)
              Actives->NBC_Vector_Send = std::vector<int>(&buf[buf[1] + 2], &buf[buf.size()]);
            std::cout << "DEAD NBC_SIZE: " << Actives->NBC_Vector_Send.size() << "\n";
            // f_comm->NBC_Vector_Send =
            // std::vector<int>(&Bcast_Vector[Bcast_Vector[1]+1],&Bcast_Vector[Bcast_Vector.size()-1]);

            std::vector<simft::Sim_FT_Istats> Deserialized_NBC;
            Deserialized_NBC = simft::Vector_to_NBC(Actives->NBC_Vector_Send);

            // first initiate outstanding non-blocking operations
            simft::Sim_FT_Perform_nb_operations(&Deserialized_NBC, Actives);
            // then free them and the already started local ones using MPI_Waitall
            simft::Sim_FT_Complete_nb_operations(&Deserialized_NBC, Actives);

            // free all "real" communicators belonging to the current custom communicator
            MPI_Comm_free(&Actives->c_comm_copy);
            MPI_Comm_free(&Actives->c_comm_copy_p2p);
            MPI_Comm_free(&Actives->c_comm_copy_coll);
            MPI_Comm_free(&Actives->c_comm_copy_probe);

            simft::Sim_FT_Current_Active_Communicators.erase(
                Actives);    // remove freed communicator from active comm set
            delete Actives;  // delete custom comm object
            Actives = nullptr;
            break;
          } else if (buf[0] == 9 || Actives->comm_revoked) {  // comm is revoked
            // at this point, switch from normal tag to reduced tag, if comm is revoked or revoke
            // message via broadcast is recieved
            Actives->ReduceTag = SIM_FT_REDUCE_TAG_REVOKED;
            Actives->DeadTag = SIM_FT_DEAD_TAG_REVOKED;
          }
        }
        simft::Sim_FT_Custom_Ireduce(Actives, true);
      } else {
        // if process is root, first receive Ireduce. After it is completed, initiate broadcast
        // (immediately finished)
        if (simft::Sim_FT_Custom_Ireduce(Actives, true)) {
          count = 4;
          Actives->CurrentMaxNBC = Actives->Last_Send_Vector[1];
          Actives->CurrentMaxNBC_Rank = Actives->Last_Send_Vector[2];
          while ((int)Actives->Dead_Processes_Root.size() < Actives->Last_Send_Vector[0]) {
            simft::Sim_FT_Receive_dead_msgs_root(Actives);
          }

          std::vector<int>
              Send_Vector;  //[0]: NextOp; [1]: Dead processes count; [>1]: Dead processes id's
          int reserveAddSize = 2;

          /*if next operation is MPI_Comm_free, request NBC vector and receive it  */
          if (Actives->Last_Send_Vector[3] == 3) {
            /* ####################
             * request and receive NBC vector
             */
            MPI_Request NBC_Request;
            MPI_Isend(0, 0, MPI_INT, Actives->CurrentMaxNBC_Rank, SIM_FT_NBC_REQUEST_TAG,
                      Actives->c_comm_copy, &NBC_Request);
            MPI_Request_free(&NBC_Request);

            MPI_Status RStatus;
            int PFlag = 0;
            while (PFlag == 0) {
              MPI_Iprobe(Actives->CurrentMaxNBC_Rank, SIM_FT_NBC_RESPONSE_TAG, Actives->c_comm_copy,
                         &PFlag, &RStatus);
              simft::Sim_FT_Perform_background_operations_dead(Actives);
            }

            int Count = 0;
            MPI_Get_count(&RStatus, MPI_INT, &Count);
            Actives->NBC_Vector_Send.resize(Count);

            MPI_Recv(&Actives->NBC_Vector_Send[0], Count, MPI_INT, Actives->CurrentMaxNBC_Rank,
                     SIM_FT_NBC_RESPONSE_TAG, Actives->c_comm_copy, MPI_STATUS_IGNORE);
            // ####################

            reserveAddSize += static_cast<int>(Actives->NBC_Vector_Send.size());
          }

          Send_Vector.reserve(Actives->Dead_Processes_Root.size() + reserveAddSize);
          // Attach the information, which operation is next (Shrink, Finalize, Comm_free,
          // Comm_agree or normal blocking collective)
          Send_Vector.push_back(Actives->Last_Send_Vector[3]);
          // Attach amount of dead processes and append all dead processes
          Send_Vector.push_back(static_cast<int>(Actives->Dead_Processes_Root.size()));
          Send_Vector.insert(Send_Vector.end(), Actives->Dead_Processes_Root.begin(),
                             Actives->Dead_Processes_Root.end());

          // attach NBC_Vector_Send to our Send_Vector
          if (Send_Vector[0] == 3) {
            Send_Vector.insert(Send_Vector.end(), Actives->NBC_Vector_Send.begin(),
                               Actives->NBC_Vector_Send.end());
          }

          simft::Sim_FT_Custom_Ibcast_Root(&Send_Vector, &(Actives->Dead_Buf), static_cast<int>(Send_Vector.size()),
                                           MPI_INT, Actives->DeadTag, Actives->c_comm_copy_coll,
                                           Actives);

          // set Ireduce_Active to false, so that a new reduce operation can be started
          Actives->Ireduce_Active = false;

          if (Actives->dead_nodes.size() != Actives->dead_set.size()) {
            Actives->dead_set = std::set<int>(
                Actives->dead_nodes.begin(),
                Actives->dead_nodes
                    .end());  // inefficient (O(n logn)) but only occurs if new processes failed
          }

          if (Send_Vector[0] == 1) {  // Shrink next
            MPI_Comm EmptyComm;
            MPI_Comm_split(Actives->c_comm, MPI_UNDEFINED, 0, &EmptyComm);
            break;
          } else if (Send_Vector[0] == 2) {  // Finalize next

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Finalize();
            exit(0);
          } else if (Send_Vector[0] == 3) {  // comm_free next

            std::vector<simft::Sim_FT_Istats> Deserialized_NBC;
            Deserialized_NBC = simft::Vector_to_NBC(Actives->NBC_Vector_Send);

            simft::Sim_FT_Perform_nb_operations(&Deserialized_NBC, Actives);
            // while(1);
            simft::Sim_FT_Complete_nb_operations(&Deserialized_NBC, Actives);

            // free all "real" communicators belonging to the current custom communicator
            MPI_Comm_free(&Actives->c_comm_copy);
            MPI_Comm_free(&Actives->c_comm_copy_p2p);
            MPI_Comm_free(&Actives->c_comm_copy_coll);
            MPI_Comm_free(&Actives->c_comm_copy_probe);

            simft::Sim_FT_Current_Active_Communicators.erase(
                Actives);    // remove freed communicator from active comm set
            delete Actives;  // delete custom comm object
            Actives = nullptr;
            break;
          } else if (Actives->comm_revoked) {  // comm is revoked
            // at this point, switch from normal tag to reduced tag, if comm is revoked or revoke
            // message via broadcast is recieved
            Actives->ReduceTag = SIM_FT_REDUCE_TAG_REVOKED;
            Actives->DeadTag = SIM_FT_DEAD_TAG_REVOKED;
          }
        }
      }
      // printf("Dead background operations \n");
      simft::Sim_FT_Perform_background_operations_dead(Actives);
    }
  }
}

/*
 * Decides (on function call) if a process should be killed right now.
 * Function can (and should) be customized.
 */
void simft::Sim_FT_decide_kill() {
  return;  //(if no simulated process failures should occur, uncomment)

  RandKill_cnt++;
  if ((simft::Sim_FT_MPI_COMM_WORLD->comm_rank == 2) && RandKill_cnt == 1) {
    simft::Sim_FT_kill_me();
  }
}

// if custom IBcast returns true, buf is received and can be read OR current rank is root; in this
// case buf should not be touched!
bool simft::Sim_FT_Custom_Ibcast_Root(std::vector<int> *buf_out, std::vector<int> *buf_copy,
                                      int count, MPI_Datatype datatype, int tag, MPI_Comm comm,
                                      simft::Sim_FT_MPI_Comm f_comm) {
  if (count != 0) {
    *buf_copy = *buf_out;  // set copy-buffer to out-buffer
  }

  for (unsigned int i = 0; i < f_comm->Bcast_Successors.size(); i++) {
    MPI_Request Send_Request;  // will never be read again because we assume that eventually all
                               // Isend operations finish successfully
    if (count == 0) {
      MPI_Isend(0, 0, datatype, f_comm->Bcast_Successors[i], tag, comm, &Send_Request);
    } else {
      MPI_Isend(&buf_copy->front(), count, datatype, f_comm->Bcast_Successors[i], tag, comm,
                &Send_Request);
    }
    MPI_Request_free(&Send_Request);
  }
  return true;
}

bool simft::Sim_FT_Custom_Ibcast_NRoot(std::vector<int> *buf_out, std::vector<int> *buf_copy,
                                       int *count, MPI_Datatype datatype, int tag, MPI_Comm comm,
                                       simft::Sim_FT_MPI_Comm f_comm) {
  int Bcast_Flag = 0;
  MPI_Status Bcast_Status;
  MPI_Iprobe(f_comm->Bcast_Predecessor, tag, comm, &Bcast_Flag, &Bcast_Status);
  if (Bcast_Flag == 1) {
    // our Ibcast can dynamically forward messages of any size, get count from probed message
    MPI_Get_count(&Bcast_Status, MPI_INT, count);
    int sendcount = *count;

    if (sendcount == 0) {
      MPI_Recv(0, 0, datatype, Bcast_Status.MPI_SOURCE, tag, comm, MPI_STATUS_IGNORE);
    } else {
      buf_out->resize(sendcount);
      MPI_Recv(&buf_out->front(), sendcount, datatype, Bcast_Status.MPI_SOURCE, tag, comm,
               MPI_STATUS_IGNORE);
      *buf_copy = *buf_out;  // set copy-buffer to out-buffer
    }

    for (unsigned int i = 0; i < f_comm->Bcast_Successors.size(); i++) {
      MPI_Request Send_Request;  // will never be read again because we assume that eventually all
                                 // Isend operations finish successfully

      if (sendcount == 0) {
        MPI_Isend(0, 0, datatype, f_comm->Bcast_Successors[i], tag, comm, &Send_Request);
      } else {
        MPI_Isend(&buf_copy->front(), sendcount, datatype, f_comm->Bcast_Successors[i], tag, comm,
                  &Send_Request);
      }
      MPI_Request_free(&Send_Request);
    }
    return true;
  } else {
    return false;
  }
}

void simft::Sim_FT_Custom_Dead_bcast(std::vector<int> *buf, simft::Sim_FT_MPI_Comm f_comm,
                                     bool revoked) {
  int BcastTag;

  // If comm is revoked, use a special tag
  if (!revoked) {
    BcastTag = SIM_FT_DEAD_TAG;
  } else {
    BcastTag = SIM_FT_DEAD_TAG_REVOKED;
  }

  if (f_comm->comm_rank == f_comm->Root_Rank) {
    f_comm->Dead_Buf = *buf;  // copy dead processes to buffer

    simft::Sim_FT_Custom_Ibcast_Root(buf, &(f_comm->Dead_Buf), static_cast<int>(buf->size()), MPI_INT, BcastTag,
                                     f_comm->c_comm_copy_coll, f_comm);

  } else {
    int count = 0;

    while (!simft::Sim_FT_Custom_Ibcast_NRoot(buf, &(f_comm->Dead_Buf), &count, MPI_INT, BcastTag,
                                              f_comm->c_comm_copy_coll, f_comm)) {
      simft::Sim_FT_Perform_background_operations(true);
    }
  }
  simft::Sim_FT_Perform_background_operations(true);
}

void simft::customReduceOp(void *inV, void *inoutV, int *len, MPI_Datatype *dptr) {
  int *in = (int *)inV;
  int *inout = (int *)inoutV;
  if (*(inout + 1) < *(in + 1)) {  // check, if NBC count and id has to be updated (we want the max.
                                   // NBC count and its proc id)
    *(inout + 1) = *(in + 1);
    *(inout + 2) = *(in + 2);
  }
  *inout = *inout + *in;  // add dead count from successor to local dead count
  if (*(inout + 3) == 0) {
    *(inout + 3) = *(in + 3);
  } else if (*(inout + 3) == 5) {  // 4 => false in comm_agree; 5 => true in comm_agree
    *(in + 3) = *(inout + 3);  // comm_agree performs "and"-operation => if local value is true, set
                               // it to the remote value (true or false)
  }
}

// Reduce Vector has the following form: 4x int: ( dead count; max. NBC; rank max. NBC; NextOp )
int simft::Sim_FT_Custom_Reduce(int Send_Vector[], simft::Sim_FT_MPI_Comm f_comm,
                                simft::NextOp option, bool Revoked) {
  //#ifndef USE_NON_BLOCKING_COLLECTIVES
  int Recv_Flag = 0;
  unsigned int Received_Cnt = 0;
  int count = 4;

  int ReduceTag;

  if (!Revoked) {
    ReduceTag = SIM_FT_REDUCE_TAG;
  } else {
    ReduceTag = SIM_FT_REDUCE_TAG_REVOKED;
  }
  // std::cout << "Reduce Tag " << ReduceTag <<" of rank: " << f_comm->comm_rank << " revoked: " <<
  // Revoked << "\n";
  MPI_Comm comm = f_comm->c_comm_copy_coll;

  Send_Vector[0] = 0;
  Send_Vector[1] = static_cast<int>(f_comm->Recent_Icollectives.size());
  Send_Vector[2] = f_comm->comm_rank;
  Send_Vector[3] = (int)option;  // if the next operation shall be Shrink or Finalize etc., this
                                 // must be propagated through the comm to inform the dead processes
                                 // about it

  MPI_Status Reduce_Status;

  // if a successor message is received, set the corresponding Received to true
  bool *Received;
  if (f_comm->Bcast_Successors.size() > 0) {
    Received = new bool[f_comm->Bcast_Successors.size()];
  }

  for (unsigned int i = 0; i < f_comm->Bcast_Successors.size(); i++) {
    Received[i] = false;
  }

  std::vector<int> Receive_Vectors;
  Receive_Vectors.resize(count);

  while (Received_Cnt < f_comm->Bcast_Successors.size()) {
    for (unsigned int i = 0; i < f_comm->Bcast_Successors.size(); i++) {
      if (!Received[i]) {
        MPI_Iprobe(f_comm->Bcast_Successors[i], ReduceTag, comm, &Recv_Flag, &Reduce_Status);
        // std::cout << "Bcast Successor " << f_comm->Bcast_Successors[i] <<" of rank " <<
        // f_comm->comm_rank << "\n";
        if (Recv_Flag == 1) {
          if (count ==
              0) {  //"deprecated", should not happen, because size of reduce msg is at least 3 int!
            std::cout << ">>> Sim_FT note: empty reduce message received! <<<\n";
            MPI_Recv(0, count, MPI_INT, f_comm->Bcast_Successors[i], ReduceTag, comm,
                     MPI_STATUS_IGNORE);
          } else {
            MPI_Recv(&Receive_Vectors[0], count, MPI_INT, f_comm->Bcast_Successors[i], ReduceTag,
                     comm, MPI_STATUS_IGNORE);
          }

          simft::customReduceOp(&Receive_Vectors[0], &Send_Vector[0]);

          Received[i] = true;
          Received_Cnt++;
        }
      }
    }
    simft::Sim_FT_Perform_background_operations(true);
  }

  if (f_comm->Bcast_Successors.size() > 0) {
    // clean up after operation is finished
    delete[] Received;
  }

  if (f_comm->comm_rank != f_comm->Root_Rank) {
    // std::cout << "sending from rank " << f_comm->comm_rank << " to: " <<
    // f_comm->Bcast_Predecessor << " size: " << f_comm->comm_size << " revoked: " << Revoked
    // <<"\n";
    MPI_Request Send_Request;
    // use the Bcast-topology to find a predecessor
    MPI_Isend(&Send_Vector[0], 4, MPI_INT, f_comm->Bcast_Predecessor, ReduceTag, comm,
              &Send_Request);
    MPI_Request_free(&Send_Request);
  }

  if (f_comm->comm_revoked) {
    return MPI_ERR_REVOKED;
  } else {
    return MPI_SUCCESS;
  }
}

// If a process is dead, it needs to continuously forward reduce-messages, as well as add the own
// dead-message
bool simft::Sim_FT_Custom_Ireduce(simft::Sim_FT_MPI_Comm f_comm,
                                  bool proc_dead) {  //, std::vector<std::vector<int>>
                                                     //&Receive_Vectors, std::vector<bool>
                                                     //&Received, std::vector<MPI_Request>
                                                     //&Recv_Request){

  // explanation: if a reduce message is already forwarded to the predecessor, a new reduce mustn't
  // be startet until the bcast is recieved and Ireduce_Active is set to false
  if (f_comm->Ireduce_Active == true) return false;
  int Recv_Flag = 0;
  MPI_Status Reduce_Status;
  int ReduceTag = f_comm->ReduceTag;
  /*
          if( !f_comm->comm_revoked){
      ReduceTag = SIM_FT_REDUCE_TAG;
    }else{
      ReduceTag = SIM_FT_REDUCE_TAG_REVOKED;
    }*/
  // std::cout << "Reduce Tag: " << ReduceTag << "\n";

  // check all successors for incoming reduce messages
  for (unsigned int i = 0; i < f_comm->Bcast_Successors.size(); i++) {
    if (!f_comm->Received[i]) {
      MPI_Iprobe(f_comm->Bcast_Successors[i], ReduceTag, f_comm->c_comm_copy_coll, &Recv_Flag,
                 &Reduce_Status);

      if (Recv_Flag == 1) {
        std::vector<int> Receive_Vector;
        Receive_Vector.resize(4);
        /*if(f_comm->comm_rank == 1){
          std::cout << "received from " << f_comm->Bcast_Successors[i] << " size: "<<
        f_comm->comm_size << "\n";
        }*/
        MPI_Recv(&Receive_Vector[0], 4, MPI_INT, f_comm->Bcast_Successors[i], ReduceTag,
                 f_comm->c_comm_copy_coll, MPI_STATUS_IGNORE);

        // call reduce function
        simft::customReduceOp(&Receive_Vector[0], &f_comm->Send_Vector[0]);

        // after a successor message is received, set Received to true
        f_comm->Received[i] = true;
        f_comm->Received_Cnt++;
      }
    }
  }

  simft::Sim_FT_Check_comm_revoked(f_comm);
  // std::cout << "Receive count " << f_comm->Received_Cnt << " Successor count " <<
  // f_comm->Bcast_Successors.size()<< " rank: " << f_comm->comm_rank << " size: " <<
  // f_comm->comm_size <<"\n";

  if (f_comm->Received_Cnt < f_comm->Bcast_Successors.size()) {
    // not all successors have sent reduce messages, wait
    return false;
  } else {
    if (proc_dead) {
      f_comm->Send_Vector[0]++;
    }

    f_comm->Last_Send_Vector = f_comm->Send_Vector;
    if (f_comm->comm_rank != f_comm->Root_Rank) {  // if root, process has to access
                                                   // f_comm->Last_Send_Vector to get all failed
                                                   // processes
      MPI_Request Send_Request;

      // use the Bcast-topology to get the processes' predecessor
      MPI_Isend(&f_comm->Last_Send_Vector[0], static_cast<int>(f_comm->Last_Send_Vector.size()), MPI_INT,
                f_comm->Bcast_Predecessor, ReduceTag, f_comm->c_comm_copy_coll, &Send_Request);
      MPI_Request_free(&Send_Request);
    }

    f_comm->Ireduce_Active = true;
    // std::cout << "Active Ireduce" << f_comm << " " << f_comm->comm_rank <<" rvoked: " <<
    // f_comm->comm_revoked <<" size " << f_comm->comm_size <<"\n";
    // std::cout << "sending to " << f_comm->Bcast_Predecessor <<" size_ " << f_comm->comm_size
    // <<"\n";
    // std::cout << "Send command is " << f_comm->Send_Vector[3] << "\n";
    /*
     * Custom Ireduce is finished, clear variables for a possible next Ireduce
     */
    for (unsigned int i = 0; i < f_comm->Bcast_Successors.size(); i++) {
      f_comm->Received[i] = false;
    }
    f_comm->Received_Cnt = 0;

    // if(f_comm->Root_Rank == f_comm->comm_rank)
    // f_comm->Last_Send_Vector = f_comm->Send_Vector;

    f_comm->Send_Vector[0] = 0;
    f_comm->Send_Vector[1] = static_cast<int>(f_comm->Recent_Icollectives.size());
    f_comm->Send_Vector[2] = f_comm->comm_rank;
    f_comm->Send_Vector[3] = 0;

    return true;
  }
}

std::vector<int> simft::Bcast_Get_Successors(int root_node, int node_id, int tree_size,
                                             int Successor_Count) {
  if (Successor_Count == -1) {
    Successor_Count = TREE_SUCCESSOR_COUNT;
  }
  std::vector<int> ReturnVec;
  for (int i = 0; i < Successor_Count; i++) {
    // normalize id, so that root has the id 0 and calculate the successor
    int Normalized =
        (simft::ConvertId(simft::NormalizeId(node_id, root_node), tree_size) * Successor_Count + i +
         1);
    if (Normalized < tree_size) {
      // revert the normalization, so that root has its original id
      ReturnVec.push_back(simft::ConvertId(simft::DeNormalizeId(Normalized, root_node), tree_size));
    }
  }
  return ReturnVec;
}

int simft::Bcast_Get_Predecessor(int root_node, int node_id, int tree_size, int Successor_Count) {
  if (Successor_Count == -1) {
    Successor_Count = TREE_SUCCESSOR_COUNT;
  }
  int ReturnInt;
  int Normalized = (ConvertId(NormalizeId(node_id, root_node), tree_size));
  ReturnInt = DeNormalizeId(
      (Normalized - 1 - ((Normalized - 1) % Successor_Count)) / Successor_Count, root_node);
  return ReturnInt;
}

int simft::ConvertId(int id, int tree_size) {
  if (id >= tree_size) {
    return id - tree_size;
  } else if (id < 0) {
    return tree_size + id;
  } else {
    return id;
  }
}

int simft::NormalizeId(int id, int root) { return id - root; }

int simft::DeNormalizeId(int id, int root) { return id + root; }

int simft::Datatype_to_Dtype(MPI_Datatype Datatype) {
  if (Datatype == MPI_CHAR) {
    return 0;
  } else if (Datatype == MPI_INT) {
    return 1;
  } else if (Datatype == MPI_DOUBLE) {
    return 2;
  }  // if other Datatypes are used, please extend this function, as well as Perform_nb_operations
  return 1;  // should not happen
}

MPI_Datatype simft::Dtype_to_Datatype(int Datatype) {
  if (Datatype == 0) {
    return MPI_CHAR;
  } else if (Datatype == 1) {
    return MPI_INT;
  } else if (Datatype == 2) {
    return MPI_DOUBLE;
  }  // if other Datatypes are used, please extend this function, as well as Perform_nb_operations
  return MPI_INT;  // should not happen
}

int simft::MOp_to_Op(MPI_Op Operation) {
  if (Operation == MPI_MIN) {
    return 0;
  } else if (Operation == MPI_MAX) {
    return 1;
  } else if (Operation == MPI_SUM) {
    return 2;
  }          // add more MPI operations if needed
  return 0;  // should not happen
}

MPI_Op simft::Op_to_MOp(int Operation) {
  if (Operation == 0) {
    return MPI_MIN;
  } else if (Operation == 1) {
    return MPI_MAX;
  } else if (Operation == 2) {
    return MPI_SUM;
  }                // add more MPI operations if needed
  return MPI_MIN;  // should not happen
}

std::vector<int> simft::NBC_to_Vector(simft::Sim_FT_MPI_Comm f_comm) {
  std::vector<int> returnVec;
  returnVec.reserve(f_comm->Recent_Icollectives.size() * 5);

  for (unsigned int i = 0; i < f_comm->Recent_Icollectives.size(); i++) {
    // serialize the NBC informations
    returnVec.push_back(f_comm->Recent_Icollectives[i].Type);
    returnVec.push_back(f_comm->Recent_Icollectives[i].count);
    returnVec.push_back(f_comm->Recent_Icollectives[i].Datatype);
    returnVec.push_back(f_comm->Recent_Icollectives[i].root);
    returnVec.push_back(f_comm->Recent_Icollectives[i].op);
  }
  return returnVec;
}

std::vector<simft::Sim_FT_Istats> simft::Vector_to_NBC(std::vector<int> NBC_Vec) {
  std::vector<simft::Sim_FT_Istats> returnVec;
  returnVec.reserve(
      NBC_Vec.size() /
      5);  // assume our NBC_Vec has the correct form, so we don't check if the size is correct

  unsigned int i = 0;
  while (i < NBC_Vec.size()) {
    //"deserialize" the NBC informations
    simft::Sim_FT_Istats Store_Stats;
    Store_Stats.Type = NBC_Vec[i++];
    Store_Stats.count = NBC_Vec[i++];
    Store_Stats.Datatype = NBC_Vec[i++];
    Store_Stats.root = NBC_Vec[i++];
    Store_Stats.op = NBC_Vec[i++];
    returnVec.push_back(Store_Stats);
  }

  return returnVec;
}

//################# overload functions #################

bool operator==(const simft::Sim_FT_MPI_Request &r1, const simft::Sim_FT_MPI_Request_struct &r2) {
  if (r1 == nullptr) {
    if (r2.c_request == MPI_REQUEST_NULL) {
      return true;
    } else {
      return false;
    }
  }
  return (r1->c_request == r2.c_request);  // TODO maybe need to compare more?
}

bool operator==(const simft::Sim_FT_MPI_Request_struct &r1, const simft::Sim_FT_MPI_Request &r2) {
  return (r2 == r1);
}

bool operator!=(const simft::Sim_FT_MPI_Request &r1, const simft::Sim_FT_MPI_Request_struct &r2) {
  if (r1 == nullptr) {
    if (r2.c_request != MPI_REQUEST_NULL) {
      return true;
    } else {
      return false;
    }
  }
  return (r1->c_request != r2.c_request);  // TODO maybe need to compare more?
}

bool operator!=(const simft::Sim_FT_MPI_Request_struct &r1, const simft::Sim_FT_MPI_Request &r2) {
  return (r2 != r1);
}

bool operator==(const simft::Sim_FT_MPI_Comm &c1, const simft::Sim_FT_MPI_Comm_struct &c2) {
  if (c1 == nullptr) {
    if (c2.c_comm == MPI_COMM_NULL) {
      return true;
    } else {
      return false;
    }
  }
  return (c1->c_comm == c2.c_comm);  // TODO maybe need to compare more?
}

bool operator==(const simft::Sim_FT_MPI_Comm_struct &c1, const simft::Sim_FT_MPI_Comm &c2) {
  return (c2 == c1);
}

bool operator!=(const simft::Sim_FT_MPI_Comm &c1, const simft::Sim_FT_MPI_Comm_struct &c2) {
  if (c1 == nullptr) {
    if (c2.c_comm != MPI_COMM_NULL) {
      return true;
    } else {
      return false;
    }
  }
  return (c1->c_comm != c2.c_comm);  // TODO maybe need to compare more?
}

bool operator!=(const simft::Sim_FT_MPI_Comm_struct &c1, const simft::Sim_FT_MPI_Comm &c2) {
  return (c2 != c1);  // TODO maybe need to compare more?
}
