#include "MPI-FT.h"
#include REAL_MPI_INCLUDE

int simft::Sim_FT_MPI_Cart_create(simft::Sim_FT_MPI_Comm comm_old, int ndims, int dims[],
                                  int periods[], int reorder, simft::Sim_FT_MPI_Comm *comm_cart) {
  simft::Sim_FT_decide_kill();
  simft::Sim_FT_Perform_background_operations();

  if (simft::Sim_FT_Check_dead_processes(comm_old)) {
    return MPI_ERR_PROC_FAILED;
  }

  (*comm_cart) = new simft::Sim_FT_MPI_Comm_struct;

  int Ret = MPI_Cart_create(comm_old->c_comm, ndims, dims, periods, reorder, &(*comm_cart)->c_comm);

  if (MPI_SUCCESS == Ret) {
    simft::Sim_FT_Initialize_new_comm((comm_cart),
                                      true);  // TODO: if "reroot" implemented, set to false
  } else {
    delete (*comm_cart);
    (*comm_cart) = nullptr;
  }

  return Ret;
}
