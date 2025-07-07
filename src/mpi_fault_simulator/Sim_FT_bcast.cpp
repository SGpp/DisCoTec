#include "discotec/MPI-FT.h"
#include REAL_MPI_INCLUDE
#include <iostream>

int simft::Sim_FT_MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root,
                            simft::Sim_FT_MPI_Comm f_comm) {
  simft::Sim_FT_decide_kill();
  simft::Sim_FT_Perform_background_operations();

  int ErrLayer = simft::Sim_FT_Check_dead_processes(f_comm);
  if (MPI_ERR_PROC_FAILED == ErrLayer) {
    return MPI_ERR_PROC_FAILED;
  } else if (MPI_ERR_REVOKED == ErrLayer) {
    return MPI_ERR_REVOKED;
  }

  return MPI_Bcast(buffer, count, datatype, root, f_comm->c_comm);
}
