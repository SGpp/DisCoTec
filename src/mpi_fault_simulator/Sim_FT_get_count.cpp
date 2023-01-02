#include "MPI-FT.h"
#include REAL_MPI_INCLUDE

//
int simft::Sim_FT_MPI_Get_count(const simft::Sim_FT_MPI_Status *status, MPI_Datatype datatype,
                                int *count) {
  *count = status->count;
  return 0;
}
