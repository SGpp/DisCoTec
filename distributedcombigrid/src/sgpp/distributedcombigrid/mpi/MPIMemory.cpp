#include "sgpp/distributedcombigrid/mpi/MPIMemory.hpp"

#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"

namespace combigrid {
namespace mpimemory {
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

int get_all_memory_usage_kb(unsigned long* vmrss, unsigned long* vmsize, int np,
                            CommunicatorType comm) {

  unsigned long vmrss_kb, vmsize_kb;
  int ret_code = get_memory_usage_kb(&vmrss_kb, &vmsize_kb);
  if (ret_code != 0) {
    printf("Could not gather memory usage!\n");
    return ret_code;
  }

  MPI_Allreduce(&vmrss_kb, vmrss, 1, MPI_UNSIGNED_LONG, MPI_SUM, comm);

  MPI_Allreduce(&vmsize_kb, vmsize, 1, MPI_UNSIGNED_LONG, MPI_SUM, comm);
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
    printf("\n ,local memory usage (KB): %3ld, %6ld, %6ld\n", theMPISystem()->getNumProcs(), local_vmrss, local_vmsize);
  }
}

void print_memory_usage_world() {
  unsigned long global_vmrss, global_vmsize;
  auto commSize = getCommSize(MPI_COMM_WORLD);
  // get_all_memory_usage_kb(&global_vmrss, &global_vmsize,
  // getCommSize(theMPISystem()->getWorldComm()), theMPISystem()->getWorldComm());
  get_all_memory_usage_kb(&global_vmrss, &global_vmsize, commSize,
                          MPI_COMM_WORLD);
  if (getCommRank(MPI_COMM_WORLD) == 0) {
    printf("\n ,world memory usage (KB): %3ld, %6ld, %6ld\n", commSize, global_vmrss,
           global_vmsize);
  }
}

}  // namespace mpimemory
}  // namespace combigrid
