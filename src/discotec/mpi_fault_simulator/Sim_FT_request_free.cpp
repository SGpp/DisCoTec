#include "MPI-FT.h"
#include REAL_MPI_INCLUDE

int simft::Sim_FT_MPI_Request_free(simft::Sim_FT_MPI_Request *request) {
  int Ret = MPI_Request_free(&(*request)->c_request);
  delete (*request);
  *request = nullptr;  // = &simft::Sim_FT_MPI_REQUEST_NULL;
  return Ret;
}
