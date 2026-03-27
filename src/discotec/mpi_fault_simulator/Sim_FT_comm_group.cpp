#include "MPI-FT.h"
#include REAL_MPI_INCLUDE

int simft::Sim_FT_MPI_Comm_group(simft::Sim_FT_MPI_Comm f_comm, MPI_Group *group) {
  simft::Sim_FT_decide_kill();
  simft::Sim_FT_Perform_background_operations();

  return MPI_Comm_group(f_comm->c_comm, group);
}
