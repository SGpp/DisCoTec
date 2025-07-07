#include "discotec/MPI-FT.h"
#include REAL_MPI_INCLUDE

int simft::Sim_FT_MPI_Comm_size(simft::Sim_FT_MPI_Comm f_comm, int *size) {
  simft::Sim_FT_decide_kill();
  simft::Sim_FT_Perform_background_operations();

  return MPI_Comm_size(f_comm->c_comm, size);
}
