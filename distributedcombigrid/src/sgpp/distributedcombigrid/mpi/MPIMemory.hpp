#ifndef SRC_SGPP_COMBIGRID_MPI_MEMORY_HPP_
#define SRC_SGPP_COMBIGRID_MPI_MEMORY_HPP_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"

/**
 * @brief utility to get memory usage in the MPI processes. slow!
 * combined from the examples at
 * https://hpcf.umbc.edu/general-productivity/checking-memory-usage/
 */

namespace combigrid {
/*
 * Look for lines in the procfile contents like:
 * VmRSS:         5560 kB
 * VmSize:         5560 kB
 *
 * Grab the number between the whitespace and the "kB"
 * If 1 is returned in the end, there was a serious problem
 * (we could not find one of the memory usages)
 */
int get_memory_usage_kb(unsigned long* vmrss_kb, unsigned long* vmsize_kb) {
  /* Get the the current process' status file from the proc filesystem */
  FILE* procfile = fopen("/proc/self/status", "r");

  long to_read = 8192;
  char buffer[to_read];
  int read = fread(buffer, sizeof(char), to_read, procfile);
  fclose(procfile);

  short found_vmrss = 0;
  short found_vmsize = 0;
  char* search_result;

  /* Look through proc status contents line by line */
  char delims[] = "\n";
  char* line = strtok(buffer, delims);

  while (line != NULL && (found_vmrss == 0 || found_vmsize == 0)) {
    search_result = strstr(line, "VmRSS:");
    if (search_result != NULL) {
      sscanf(line, "%*s %lu", vmrss_kb);
      found_vmrss = 1;
    }

    search_result = strstr(line, "VmSize:");
    if (search_result != NULL) {
      sscanf(line, "%*s %lu", vmsize_kb);
      found_vmsize = 1;
    }

    line = strtok(NULL, delims);
  }

  return (found_vmrss == 1 && found_vmsize == 1) ? 0 : 1;
}

int get_cluster_memory_usage_kb(unsigned long* vmrss_per_process, unsigned long* vmsize_per_process, int root, int np,
                                CommunicatorType comm = MPI_COMM_WORLD) {
  unsigned long vmrss_kb;
  unsigned long vmsize_kb;
  int ret_code = get_memory_usage_kb(&vmrss_kb, &vmsize_kb);

  if (ret_code != 0) {
    printf("Could not gather memory usage!\n");
    return ret_code;
  }

  MPI_Gather(&vmrss_kb, 1, MPI_UNSIGNED_LONG, vmrss_per_process, 1, MPI_UNSIGNED_LONG, root, comm);

  MPI_Gather(&vmsize_kb, 1, MPI_UNSIGNED_LONG, vmsize_per_process, 1, MPI_UNSIGNED_LONG, root,
             comm);

  return 0;
}

int get_all_memory_usage_kb(unsigned long* vmrss, unsigned long* vmsize, int np,
                            CommunicatorType comm = MPI_COMM_WORLD) {
  unsigned long vmrss_per_process[np];
  unsigned long vmsize_per_process[np];
  int ret_code = get_cluster_memory_usage_kb(vmrss_per_process, vmsize_per_process, 0, np, comm);

  if (ret_code != 0) {
    return ret_code;
  }

  *vmrss = 0;
  *vmsize = 0;
  for (int i = 0; i < np; i++) {
    *vmrss += vmrss_per_process[i];
    *vmsize += vmsize_per_process[i];
  }

  return 0;
}

int get_memory_usage_local_kb(unsigned long* local_vmrss, unsigned long* local_vmsize) {
  return get_all_memory_usage_kb(local_vmrss, local_vmsize, theMPISystem()->getNumProcs(),
                                 theMPISystem()->getLocalComm());
}

void print_memory_usage_local() {
  unsigned long local_vmrss, local_vmsize;
  get_memory_usage_local_kb(&local_vmrss, &local_vmsize);
  if (theMPISystem()->getLocalRank() == 0) {
    printf("\n local memory usage: VmRSS = %6ld KB, VmSize = %6ld KB\n", local_vmrss, local_vmsize);
  }
}

void print_memory_usage_world() {
  unsigned long global_vmrss, global_vmsize;
  // get_all_memory_usage_kb(&global_vmrss, &global_vmsize, getCommSize(theMPISystem()->getWorldComm()), theMPISystem()->getWorldComm());
  get_all_memory_usage_kb(&global_vmrss, &global_vmsize, getCommSize(MPI_COMM_WORLD), MPI_COMM_WORLD);
  WORLD_MANAGER_EXCLUSIVE_SECTION {
    printf("\n world memory usage: VmRSS = %6ld KB, VmSize = %6ld KB\n", global_vmrss, global_vmsize);
  }
}

}  // namespace combigrid
#endif /* SRC_SGPP_COMBIGRID_MPI_MEMORY_HPP_ */
