#include "MPI-FT.h"
#include REAL_MPI_INCLUDE
#include <iostream>

int simft::Sim_FT_MPI_Comm_create(simft::Sim_FT_MPI_Comm f_comm, MPI_Group group,
                                  simft::Sim_FT_MPI_Comm *newcomm) {
  simft::Sim_FT_decide_kill();
  simft::Sim_FT_Perform_background_operations();

  if (simft::Sim_FT_Check_dead_processes(f_comm)) {
    return MPI_ERR_PROC_FAILED;
  }

  (*newcomm) = new simft::Sim_FT_MPI_Comm_struct;

  int Ret = MPI_Comm_create(f_comm->c_comm, group, &(*newcomm)->c_comm);

  if (MPI_SUCCESS == Ret) {
    simft::Sim_FT_Initialize_new_comm((newcomm),
                                      true);  // TODO: if "reroot" implemented, set to false

  } else {
    // if unsuccessful, immediately delete the comm object
    delete (*newcomm);
    (*newcomm) = nullptr;  //= &simft::Sim_FT_MPI_COMM_NULL;
  }

  return Ret;
}
