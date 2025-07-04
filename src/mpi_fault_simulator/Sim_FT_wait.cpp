#include "discotec/MPI-FT.h"
#include REAL_MPI_INCLUDE
#include <iostream>

int simft::Sim_FT_MPI_Wait(simft::Sim_FT_MPI_Request *request, simft::Sim_FT_MPI_Status *status) {
  simft::Sim_FT_decide_kill();
  simft::Sim_FT_Perform_background_operations();

  int Flag = 0;
  int Ret = 0;
  while (Flag == 0) {
    // std::cout << "Test ret = " << Ret << " - req complete = " << (*request)->completed << " -
    // Flag = " << Flag << "\n";
    Ret = simft::Sim_FT_MPI_Test(request, &Flag, status);
    // std::cout << "Test ret = " << Ret << " - req complete = " << (*request)->completed << " -
    // Flag = " << Flag << "\n";
  }
  return Ret;
}
