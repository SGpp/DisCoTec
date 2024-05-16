#include <map>
#include <set>
#include <vector>

#include <mpi.h>

#ifndef MPI_FT_H
#define MPI_FT_H

// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1

namespace simft {
// how long does a process wait until it detects a simulated process failure (in p2p communication)
#define SIM_FT_TIMEOUT 0
// how long does the root wait until it sends a NBC broadcast to inform comm about a dead process
// (including non-blocking count)
#define SIM_FT_WAIT_NBC_BCAST 5

//#define DISABLE_NBC //if defined, do not use non-blocking collectives. Required for MPI
//implementations not supporting non-blocking collectives.

#define TREE_SUCCESSOR_COUNT 6

enum NextOp {
  Default = 0,
  Shrink = 1,
  Finalize = 2,
  Comm_free = 3,

  // The following values are used for MPI_Comm_agree
  AFalse = 4,
  ATrue = 5,

  /*
   * If the communicator is revoked, a special NextOp instead of Default is used,
   * making sure all processes know about the revoked comm and already initialized
   * blocking collectives are canceled. If a check_dead_processes is then called
   * with using Default instead of RevokedComm, it is immediately responded with
   * a broadcast containing NextOp RevokedComm.
  */
  RevokedComm = 9
};

struct Sim_FT_MPI_Status;
#define SIM_FT_MPI_STATUSES_IGNORE (simft::Sim_FT_MPI_Status *)1
#define SIM_FT_MPI_STATUS_IGNORE (simft::Sim_FT_MPI_Status *)1
// extern Sim_FT_MPI_Status* Sim_FT_MPI_STATUS_IGNORE;
// extern Sim_FT_MPI_Status* Sim_FT_MPI_STATUSES_IGNORE;

struct Sim_FT_MPI_Comm_struct;
typedef Sim_FT_MPI_Comm_struct *Sim_FT_MPI_Comm;

struct Sim_FT_MPI_Request_struct;
typedef Sim_FT_MPI_Request_struct *Sim_FT_MPI_Request;

typedef void(Sim_FT_MPI_Comm_errhandler_function)(Sim_FT_MPI_Comm, int *, ...);

/*####################
 * special tags and message texts for our layer functions */
#define ARE_YOU_ALIVE 1
#define I_AM_ALIVE_NB 7
//#define I_AM_DEAD_NB 8
#define I_AM_ALIVE 9
#define TAGOFFSETALIVE 10000
#define I_AM_DEAD 10
#define PROBE_REQUEST 11
#define PROBE_ANSWER 12
#define COMM_IS_REVOKED 13
#define SIM_FT_BEGIN_BROADCAST 198
#define SIM_FT_BEGIN_BROADCAST_KILL 199
#define SIM_FT_BEGIN_BROADCAST_KILL_DUMMY 200
// used for dead broadcast in Check_dead_processes
#define SIM_FT_DEAD_TAG 990
#define SIM_FT_DEAD_NBC_TAG 991
#define SIM_FT_FINALIZE_TAG 994
#define SIM_FT_COMMAND_TAG 996
#define SIM_FT_KILL_TAG 997
#define SIM_FT_STATUS_TAG 998
#define SIM_FT_STATUS_TAG_NB 999

#define SIM_FT_REVOKE_TAG 1000
#define SIM_FT_BCAST_REVOKE_TAG 1001
#define SIM_FT_BCAST_NBC_TAG 1002
// used for dead reduce in Check_dead_processes
#define SIM_FT_REDUCE_TAG 1009

// used for dead broadcast in Check_dead_processes in a revoked communicator
#define SIM_FT_DEAD_TAG_REVOKED 1012
// used for dead reduce in Check_dead_processes in a revoked communicator
#define SIM_FT_REDUCE_TAG_REVOKED 1013

#define SIM_FT_NBC_REQUEST_TAG 1020
#define SIM_FT_NBC_RESPONSE_TAG 1021

// gene implementation tags

#define FT_NEW_RANK_TAG 1030
#define FT_RECOVERY_STATUS_TAG 1031
#define FT_EXCLUDE_TAG 1032
#define FT_SHRINK_TAG 1033
#define FT_REUSABLE_RANK_TAG 1034
#define FT_FAILED_RANK_TAG 1035

/*not implemented yet:
  #define SIM_FT_CHANGE_ROOT 1
  #define SIM_FT_KEEP_ROOT 0 */
//####################

#ifndef MPI_ERR_PROC_FAILED
// New errors according to ULFM, only define if not defined by MPI yet
#define MPI_ERR_PROC_FAILED 54
#define MPI_ERR_REVOKED 55
#define MPI_ERR_PROC_FAILED_PENDING 56
#endif // MPI_ERR_PROC_FAILED

// use a custom status object, as defined in MPICH MPI.
struct Sim_FT_MPI_Status {
  MPI_Status c_status;
  int MPI_SOURCE;
  int MPI_TAG;
  int MPI_ERROR;
  int count;
  int cancelled;
  int abi_slush_fund[2];
};

/* struct is used to store all important informations about a non-blocking operation
 * that are required to be able to complete them. Has to be extended if more non-blocking ops are to
 * be used! */
struct Sim_FT_Istats {
  int Type = -1;  // 0 => ibarrier; 1 => ibcast 2 => iallreduce; 10 => Isend; 11 => Irecv
  // variables needed for ibcast:
  int count = -1;
  int Datatype = -1;
  int root = -1;
  // + needed for iallreduce:
  int op = -1;

  // + needed for Isend and Irecv
  int tag = -1;

  // special variable for p2p. Will hold the "partner" rank from a p2p-communication.
  int p2p_rank = -1;

  /* Stores the pointer to the request belonging to the current non-blocking operation.
   * If the operations are to be completed by our fault layer due to a failed process,
   * this request is used to complete the operation if not already completed
   */
  simft::Sim_FT_MPI_Request request;

  /* special variables to hold the buffer to complete NBC operations with our completion protocol.
   */
  std::vector<char> BuffChar;
  std::vector<int> BuffInt;
  std::vector<double> BuffDouble;
  // another version for iallreduce
  std::vector<char> BuffChar2;
  std::vector<int> BuffInt2;
  std::vector<double> BuffDouble2;
};

struct Sim_FT_MPI_Comm_struct {
  MPI_Comm c_comm = MPI_COMM_NULL;  //"real" comm object, used by the apllication

  MPI_Comm c_comm_copy = MPI_COMM_NULL;  // copy used send/recv background messages in the layer
  MPI_Comm c_comm_copy_coll =
      MPI_COMM_NULL;  // copy used send/recv messages in the layer (used for collectives)
  MPI_Comm c_comm_copy_p2p =
      MPI_COMM_NULL;  // copy used send/recv messages in the layer (used for p2p)
  MPI_Comm c_comm_copy_probe =
      MPI_COMM_NULL;                   // copy used send/recv messages in the layer (used for probe)
  bool comm_revoked = false;           // set true if the communicator is revoked
  bool manual_revoke_message = false;  // indicates if a revoke message has to be manually forwared
                                       // (by using Sim_FT_Only_propagate_comm_revoked)

  bool Finalize_Next = false;  // set to true if MPI_Finalize is about to be called

  int comm_size = -1;
  int comm_rank = -1;

  std::vector<int> Bcast_Successors;  // in the custom broadcast used to identify the successor
                                      // nodes of the custom tree topology
  int Bcast_Predecessor;  //

  /*one process should only be root in max. one communicator at once;
  switch Root_Rank to another process if this one is already root.
  also -1 can be used to set root to last process in communicator */
  int Root_Rank = 0;

  std::vector<int> dead_nodes;  // used to mark nodes as failed (simulated)
  std::set<int>
      root_recent_dead_set;  // for nbc broadcast. attach recent dead processes not broadcasted yet
  std::set<int> dead_set;    // used to quickly get the information if a given process is already
                           // detected as dead

  bool err_return = false;  // by setting errhandler to MPI_ERRS_RETURN, this will be set to true

  // variable needed for MPI_Comm_failure_ack; set to true if current failed processes are
  // acknowledged
  bool Ack_failed_processes = false;

  /* after MPI_Comm_failure_ack, the current failed processes will be saved here.
   * MPI_Comm_failure_get_acked will read the failed group from this */
  MPI_Group failedgrp = MPI_GROUP_NULL;

  // Variable used to buffer dead processes for dead broadcast
  std::vector<int> Dead_Buf;

  // vector used to save last Non-Blocking collectives
  std::vector<Sim_FT_Istats> Recent_Icollectives;

  // used to store the serialized updated last non-blocking collectives (received from another
  // process), will be sent via broadcast
  std::vector<int> NBC_Vector_Send;

  // at a collective op, root waits for all msgs from dead processes, add dead ranks to vector
  std::vector<int> Dead_Processes_Root;

  int CurrentMaxNBC = -1;       // stores the current maximum NBC count known
  int CurrentMaxNBC_Rank = -1;  // stores the rank with the maximum NBC count

  int CurrentMinNBC = -1;       // stores the current minimum NBC count known
  int CurrentMinNBC_Rank = -1;  // stores the rank with the minimum NBC count
  bool CurrentMinNBC_Modified =
      false;  // set to true if CurrentMinNBC changed after last NBC broadcast (need new broadcast!)
  double LastNBC_bcast = 0;

  // int Blocking_Collective_Count = 0; //Needed to identify if NBC-message is up-to-date

  //    /* ############################################
  //	 * Variables used by the layer if use of non-
  //	 * blocking collectives is enabled (TODO?)
  //	 * ############################################*/
  //    MPI_Request Revoke_Request; //If comm is to be revoked, broadcast is executed by root.
  //    Request is only used once at comm creation or revoke initiation (root).
  ////TODO

  /* ##################################################
   * Variables used by a dead process in Sim_FT_kill_me
   * ################################################## */
  int DeadTag = SIM_FT_DEAD_TAG;
  int ReduceTag = SIM_FT_REDUCE_TAG;
  // variables needed to handle custom_ireduce
  int Recv_Flag = 0;
  unsigned int Received_Cnt = 0;
  int count = 0;
  std::vector<MPI_Request> Recv_Request;
  bool Ireduce_Active = false;
  // if a successor message is received, set to true
  std::vector<bool> Received;
  std::vector<int> Send_Vector;
  std::vector<int> Last_Send_Vector;
};

extern Sim_FT_MPI_Comm Sim_FT_MPI_COMM_WORLD;
extern Sim_FT_MPI_Comm Sim_FT_MPI_COMM_NULL;

// Use a customized Request-Struct be able to get the communicator, the request type and more just
// from the request
struct Sim_FT_MPI_Request_struct {
  MPI_Request c_request;

  char Request_Type = 0;  // 1 -> p2p; 2 -> collective

  bool completed = false;  // set to true if the corresponding operation has completed
  int completition_err_code = -1;

  Sim_FT_MPI_Comm f_comm;
  Sim_FT_Istats P2P_Stats;
  bool P2P_Send = false;  // indicates if a P2P request belongs to a Isend or Irecv
};

extern Sim_FT_MPI_Request_struct Sim_FT_MPI_REQUEST_NULL;

// struct Sim_FT_P2P_Response {
//	int Alive_Status = -1;
//};
typedef int Sim_FT_P2P_Response;

//#define REGULAR_MPI_COMM_WORLD MPI_COMM_WORLD
//#define REGULAR_MPI_COMM_WORLD ((MPI_Comm)0x44000000)

extern std::set<Sim_FT_MPI_Comm> Sim_FT_Current_Active_Communicators;

extern std::map<int, Sim_FT_P2P_Response> Sim_FT_Current_Active_P2P_Requests;
// extern std::map<int,Sim_FT_MPI_Request*> Sim_FT_Current_Active_Alive_Messages;

#ifdef USE_NON_BLOCKING_COLLECTIVES
extern MPI_Op customOp;
#endif
void customReduceOp(void *inV, void *inoutV, int *len = nullptr, MPI_Datatype *dptr = nullptr);

//#endif

/* ################################################
 * Custom functions to support simulated failures
 * ################################################ */

// used to send current dead processes in comm
void Sim_FT_Custom_Dead_bcast(std::vector<int> *buf, Sim_FT_MPI_Comm f_comm, bool revoked);

// used in blocking collectives to synchronize dead processes along the whole communicator
int Sim_FT_Check_dead_processes(Sim_FT_MPI_Comm f_comm, NextOp option = NextOp::Default);

// custom non-blocking broadcast used to send background informations and inside
// Check_dead_processes (root version)
bool Sim_FT_Custom_Ibcast_Root(std::vector<int> *buf_out, std::vector<int> *buf_copy, int count,
                               MPI_Datatype datatype, int tag, MPI_Comm comm,
                               Sim_FT_MPI_Comm f_comm);
// custom non-blocking broadcast used to send background informations and inside
// Check_dead_processes (non-root version)
bool Sim_FT_Custom_Ibcast_NRoot(std::vector<int> *buf_out, std::vector<int> *buf_copy, int *count,
                                MPI_Datatype datatype, int tag, MPI_Comm comm,
                                Sim_FT_MPI_Comm f_comm);

// custom reduce used inside Check_dead_processes. NextOp defines
int Sim_FT_Custom_Reduce(int *Send_Vector, Sim_FT_MPI_Comm f_comm, NextOp option = NextOp::Default,
                         bool Revoked = false);
bool Sim_FT_Custom_Ireduce(Sim_FT_MPI_Comm f_comm, bool proc_dead = true);

/* #########################################
 * functions to support custom tree topology
 * ######################################### */

/* Calculates the successor of node_id in our simple tree topology given root, current node id
 * (0,...,tree_size -1) and max. Successor_Count.
 * Successor_Count 2 would generate a binary tree. */
std::vector<int> Bcast_Get_Successors(int root_node, int node_id, int tree_size,
                                      int Successor_Count = -1);
/* Calculates the predecessor of node_id in our simple tree topology given root, current node id
 * (0,...,tree_size -1) and max. Successor_Count.
 * Successor_Count 2 would generate a binary tree. */
int Bcast_Get_Predecessor(int root_node, int node_id, int tree_size, int Successor_Count = -1);

// helper function create the custom tree topology (makes the id behave like Z/nZ with n =
// tree_size)
int ConvertId(int id, int tree_size);
int NormalizeId(int id, int root);
int DeNormalizeId(int id, int root);
void Sim_FT_Check_comm_root(Sim_FT_MPI_Comm f_comm);

/* #######################################
 * functions to call frequently to support
 * background message transfer etc.
 * ####################################### */

// check if communicator is revoked. Forward revoke message if received.
bool Sim_FT_Check_comm_revoked(Sim_FT_MPI_Comm f_comm);

/*
 * If we have a dead process and LastNBC_bcast is some time ago, broadcast the current lowest count
 * of
 * active non-blocking collectives from all dead processes known to root and a list of all recently
 * new dead processes.
 *
 * If an active process waits for a non-blocking collective to finish, possibly a dead process never
 * contributed
 * to the collective and the process waits infinitely. After receiving a NBC broadcast, the process
 * is
 * knows about it and returns the wait with a failure.
 */
void Sim_FT_NBC_bcast_root(Sim_FT_MPI_Comm f_comm);
// background broadcast containing count of recent non-blocking collectives and list of recently
// killed processes (non-root version)
int Sim_FT_NBC_bcast_Nroot(Sim_FT_MPI_Comm f_comm);
bool Sim_FT_Only_check_comm_revoked(Sim_FT_MPI_Comm f_comm);
void Sim_FT_Only_propagate_comm_revoked(Sim_FT_MPI_Comm f_comm);
// if a process dies, a dead message is sent to root. This function will receive these messages.
bool Sim_FT_Receive_dead_msgs_root(Sim_FT_MPI_Comm f_comm);
// perform all important background operations for a specific communicator
void Sim_FT_Perform_background_operations(Sim_FT_MPI_Comm f_comm, bool inside_Check_dead = false);
// perform all important background operations for all communicators the process is active
void Sim_FT_Perform_background_operations(bool inside_Ckeck_dead = false);
// perform all important background operations for a specific communicator (special vesion for dead
// processes)
void Sim_FT_Perform_background_operations_dead(Sim_FT_MPI_Comm f_comm);

// important function. New communicators have to be initialized by our fault layer
void Sim_FT_Initialize_new_comm(Sim_FT_MPI_Comm *f_new_comm, bool firstInit = false);

// if a process wants the communicator to be revoked, it calls this function and a message is sent
// to root
void Sim_FT_Send_revoke_message(Sim_FT_MPI_Comm f_comm);

// if a node wants to simulated kill itself, this function has to be called
void Sim_FT_kill_me();

// decides, when a process should get killed (customizable, for example randomly or after a certain
// amount of MPI calls)
void Sim_FT_decide_kill();

/* Function initiates some Non-Blocking-Collective operations.
 * It is used to free ranks who already began Non-Blocking Collectives just before a process died.
 * Begin/Finish all outstanding Collective operations with dummy data.
 *
 * Note: This function (as well as the struct simft::Sim_FT_Istats) has to
 *       be extended depending on which Non-Blocking operations and datatypes are used.
 * */
void Sim_FT_Perform_nb_operations(std::vector<Sim_FT_Istats> *NB_ops_to_perform,
                                  Sim_FT_MPI_Comm f_comm);

/* Function completes all outstanding operations in a communicator
 * stored in Recent_Icollectives and NBC_Vector_Send in the custom
 * communicator f_comm.
 * Should only be called after Sim_FT_Perform_nb_operations has been
 * called at all ranks in the communicator with the same MPI
 * operations (including the dead processes in the comm). Otherwise
 * the function will loop indefinitely.
 * */
void Sim_FT_Complete_nb_operations(std::vector<simft::Sim_FT_Istats> *Deserialized_NBC_Vector,
                                   Sim_FT_MPI_Comm f_comm);

/* #############################################
 * functions to override "classic" MPI functions
 * ############################################# */

int Sim_FT_MPI_Abort(Sim_FT_MPI_Comm f_comm, int errorcode);
int Sim_FT_MPI_Allreduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                         Sim_FT_MPI_Comm f_comm);
int Sim_FT_MPI_Barrier(Sim_FT_MPI_Comm f_comm);
int Sim_FT_MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root,
                     Sim_FT_MPI_Comm f_comm);
int Sim_FT_MPI_Cart_coords(Sim_FT_MPI_Comm f_comm, int rank, int maxdims, int coords[]);
int Sim_FT_MPI_Cart_create(Sim_FT_MPI_Comm comm_old, int ndims, int dims[], int periods[],
                           int reorder, Sim_FT_MPI_Comm *comm_cart);
int Sim_FT_MPI_Cart_rank(Sim_FT_MPI_Comm f_comm, int coords[], int *rank);
int Sim_FT_MPI_Comm_create(Sim_FT_MPI_Comm f_comm, MPI_Group group, Sim_FT_MPI_Comm *newcomm);
int Sim_FT_MPI_Comm_free(Sim_FT_MPI_Comm *comm);
int Sim_FT_MPI_Comm_free2(Sim_FT_MPI_Comm *comm);

int Sim_FT_MPI_Comm_rank(Sim_FT_MPI_Comm f_comm, int *rank);
int Sim_FT_MPI_Comm_group(Sim_FT_MPI_Comm f_comm, MPI_Group *group);
int Sim_FT_MPI_Comm_set_errhandler(Sim_FT_MPI_Comm f_comm, MPI_Errhandler errhandler);
int Sim_FT_MPI_Comm_size(Sim_FT_MPI_Comm f_comm, int *size);
int Sim_FT_MPI_Comm_split(Sim_FT_MPI_Comm f_comm, int color, int key, Sim_FT_MPI_Comm *f_newcomm);
int Sim_FT_MPI_Finalize(void);
void Sim_FT_MPI_Finalize_worker(void);
int Sim_FT_MPI_Get_count(const Sim_FT_MPI_Status *status, MPI_Datatype datatype, int *count);
#ifndef DISABLE_NBC
int Sim_FT_MPI_Iallreduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                          Sim_FT_MPI_Comm f_comm, Sim_FT_MPI_Request *request);
int Sim_FT_MPI_Ibarrier(Sim_FT_MPI_Comm f_comm, Sim_FT_MPI_Request *request);
int Sim_FT_MPI_Ibcast(void *buffer, int count, MPI_Datatype datatype, int root,
                      Sim_FT_MPI_Comm f_comm, Sim_FT_MPI_Request *request);
#endif
int Sim_FT_MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
                     Sim_FT_MPI_Comm f_comm, Sim_FT_MPI_Request *request);
int Sim_FT_MPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                     Sim_FT_MPI_Comm f_comm, Sim_FT_MPI_Request *request);
int Sim_FT_MPI_Probe(int source, int tag, Sim_FT_MPI_Comm f_comm, Sim_FT_MPI_Status *status);
int Sim_FT_MPI_Recv(void *buf, int count, MPI_Datatype type, int source, int tag,
                    Sim_FT_MPI_Comm comm, Sim_FT_MPI_Status *status);
int Sim_FT_MPI_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
                      int root, Sim_FT_MPI_Comm f_comm);
int Sim_FT_MPI_Request_free(Sim_FT_MPI_Request *request);
int Sim_FT_MPI_Send(void *buf, int count, MPI_Datatype type, int dest, int tag,
                    Sim_FT_MPI_Comm comm);
int Sim_FT_MPI_Test(Sim_FT_MPI_Request *request, int *flag, Sim_FT_MPI_Status *status);
int Sim_FT_MPI_Wait(Sim_FT_MPI_Request *request, Sim_FT_MPI_Status *status);
int Sim_FT_MPI_Waitall(int count, Sim_FT_MPI_Request array_of_requests[],
                       Sim_FT_MPI_Status array_of_statuses[]);
int Sim_FT_MPI_Init(int *argc, char ***argv);

void Sim_FT_MPI_Init_worker();

// Custom helper-functions
int Sim_FT_Collective_failure_propagation(Sim_FT_MPI_Comm f_comm);
bool Sim_FT_proc_is_root();

/* We need to serialize MPI_Datatype for our perform_NB_operations.
 * These functions transform ints to Datatypes (and the other way around)
 */
int Datatype_to_Dtype(MPI_Datatype Datatype);
MPI_Datatype Dtype_to_Datatype(int Datatype);

int MOp_to_Op(MPI_Op Operation);
MPI_Op Op_to_MOp(int Operation);

std::vector<int> NBC_to_Vector(simft::Sim_FT_MPI_Comm f_comm);
std::vector<Sim_FT_Istats> Vector_to_NBC(std::vector<int> NBC_Vec);

/*
 * Functions needed to produce bit flips
 */

// Function needed to be able to simulate bit errors; special case for integer arrays. p specifies
// the flip probability of a bit
void Sim_FT_Manipulate_bits(int *BitArray, int ArraySize, double p);

// use the c++ internal function to calculate for a binomial distribution the amount of successful
// events (in our case: bit flips)
int getAmountOfFlips(int n, double p);

// flip a specific bit of a given integer array
void flipBit(int *BitArray, int BitPos);

}  // namespace simft

/* ##############
 * ULFM functions
 * ##############
 * */

int MPI_Comm_revoke(simft::Sim_FT_MPI_Comm f_comm);
int MPI_Comm_shrink(simft::Sim_FT_MPI_Comm f_comm, simft::Sim_FT_MPI_Comm *newcomm);
int MPI_Comm_failure_ack(simft::Sim_FT_MPI_Comm f_comm);
int MPI_Comm_failure_get_acked(simft::Sim_FT_MPI_Comm f_comm, MPI_Group *failedgrp);
int MPI_Comm_agree(simft::Sim_FT_MPI_Comm f_comm, int *flag);

/* because we override MPI_Request and MPI_Comm with custom definitions
 * in case one wants to check e.g. if Request == MPI_REQUEST_NULL,
 * we need to provide custom operators == and !=
 */
bool operator==(const simft::Sim_FT_MPI_Request &r1, const simft::Sim_FT_MPI_Request_struct &r2);
bool operator==(const simft::Sim_FT_MPI_Request_struct &r1, const simft::Sim_FT_MPI_Request &r2);
bool operator!=(const simft::Sim_FT_MPI_Request &r1, const simft::Sim_FT_MPI_Request_struct &r2);
bool operator!=(const simft::Sim_FT_MPI_Request_struct &r1, const simft::Sim_FT_MPI_Request &r2);
bool operator==(const simft::Sim_FT_MPI_Comm &c1, const simft::Sim_FT_MPI_Comm_struct &c2);
bool operator==(const simft::Sim_FT_MPI_Comm_struct &c1, const simft::Sim_FT_MPI_Comm &c2);
bool operator!=(const simft::Sim_FT_MPI_Comm &c1, const simft::Sim_FT_MPI_Comm_struct &c2);
bool operator!=(const simft::Sim_FT_MPI_Comm_struct &c1, const simft::Sim_FT_MPI_Comm &c2);

#endif /* MPI_FT_H */
