#include "discotec/MPI-FT.h"
#include REAL_MPI_INCLUDE

int simft::Sim_FT_MPI_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                             MPI_Op op, int root, simft::Sim_FT_MPI_Comm f_comm) {
  simft::Sim_FT_decide_kill();
  simft::Sim_FT_Perform_background_operations();

  if (simft::Sim_FT_Check_dead_processes(f_comm)) {
    return MPI_ERR_PROC_FAILED;
  }

  return MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, f_comm->c_comm);
}
