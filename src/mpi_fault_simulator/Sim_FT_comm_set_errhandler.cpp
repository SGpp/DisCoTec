#include "../../include/discotec/MPI-FT.h"
#include REAL_MPI_INCLUDE

#include <iostream>
#include <vector>

int simft::Sim_FT_MPI_Comm_set_errhandler(simft::Sim_FT_MPI_Comm f_comm,
                                          MPI_Errhandler errhandler) {
  if (errhandler == MPI_ERRORS_RETURN) {
    f_comm->err_return = true;
    return 0;
  } else {
    int ret = MPI_Comm_set_errhandler(f_comm->c_comm, errhandler);
    simft::Sim_FT_decide_kill();
    simft::Sim_FT_Perform_background_operations();
    return ret;
  }
}
