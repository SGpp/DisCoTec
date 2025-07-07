#include "discotec/MPI-FT.h"
#include REAL_MPI_INCLUDE

int simft::Sim_FT_MPI_Comm_split(simft::Sim_FT_MPI_Comm f_comm, int color, int key,
                                 simft::Sim_FT_MPI_Comm *f_newcomm) {
  simft::Sim_FT_decide_kill();
  simft::Sim_FT_Perform_background_operations();

  if (simft::Sim_FT_Check_dead_processes(f_comm)) {
    return MPI_ERR_PROC_FAILED;
  }

  (*f_newcomm) = new simft::Sim_FT_MPI_Comm_struct;

  int Ret = MPI_Comm_split(f_comm->c_comm, color, key, &((*f_newcomm)->c_comm));

  if (MPI_SUCCESS == Ret) {
    simft::Sim_FT_Initialize_new_comm((f_newcomm),
                                      true);  // TODO: if "reroot" implemented, set to false
  } else {
    delete (*f_newcomm);
    (*f_newcomm) = simft::Sim_FT_MPI_COMM_NULL;
  }

  return Ret;
}
