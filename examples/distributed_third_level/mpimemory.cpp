/*
 * mpimemory.cpp
 */
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include "discotec/mpi/MPIMemory.hpp"

#include <mpi.h>

#include <iostream>
#include <string>
#include <vector>

#include "discotec/io/ParameterIO.hpp"
#include "discotec/mpi/MPISystem.hpp"
#include "discotec/utils/Stats.hpp"

using namespace combigrid;

int main(int argc, char** argv) {
  [[maybe_unused]] auto mpiOnOff = MpiOnOff(&argc, &argv);

  mpimemory::print_memory_usage_world();
  /* when using timers (TIMING is defined in Stats), the Stats class must be
   * initialized at the beginning of the program. (and finalized in the end)
   */
  Stats::initialize();
  mpimemory::print_memory_usage_world();
  // read ctparam: rank 0 reads and broadcasts
  std::string paramfile = "ctparam";
  if (argc > 1) paramfile = argv[1];
  auto [ngroup, nprocs, cfg] = combigrid::readParameterFile(paramfile, MPI_COMM_WORLD);

  // divide the MPI processes into process group and initialize the
  // corresponding communicators
  theMPISystem()->init(ngroup, nprocs);
  WORLD_MANAGER_EXCLUSIVE_SECTION { std::cout << "after mpi system init " << std::flush; }

  mpimemory::print_memory_usage_world();

  Stats::finalize();

  /* write stats to json file for postprocessing */
  Stats::write("timers_memory.json");

  WORLD_MANAGER_EXCLUSIVE_SECTION { std::cout << "after everything" << std::flush; }
  mpimemory::print_memory_usage_world();

  return 0;
}
