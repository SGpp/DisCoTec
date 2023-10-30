#include "mpi/MPIMemory.hpp"

#include "mpi/MPISystem.hpp"

#include <vector>

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
  std::vector<char> buffer(to_read);
  fread(buffer.data(), sizeof(char), to_read, procfile);
  fclose(procfile);

  short found_vmrss = 0;
  short found_vmsize = 0;
  char* search_result;

  /* Look through proc status contents line by line */
  char delims[] = "\n";
  char* line = strtok(buffer.data(), delims);

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

int get_all_memory_usage_kb(unsigned long* vmrss, unsigned long* vmsize, CommunicatorType comm) {
  MPI_Barrier(comm);
  unsigned long vmbuf[2]; /* vmrss_kb, vmsize_kb */
  int ret_code = get_memory_usage_kb(&vmbuf[0], &vmbuf[1]);
  if (ret_code != 0) {
    printf("Could not gather memory usage!\n");
    return ret_code;
  }

  MPI_Allreduce(MPI_IN_PLACE, vmbuf, 2, MPI_UNSIGNED_LONG, MPI_SUM, comm);

  *vmrss = vmbuf[0];
  *vmsize = vmbuf[1];
  return 0;
}

int get_memory_usage_local_kb(unsigned long* local_vmrss, unsigned long* local_vmsize) {
  return get_all_memory_usage_kb(local_vmrss, local_vmsize, theMPISystem()->getLocalComm());
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
  // get_all_memory_usage_kb(&global_vmrss, &global_vmsize, theMPISystem()->getWorldComm());
  get_all_memory_usage_kb(&global_vmrss, &global_vmsize, MPI_COMM_WORLD);
  if (getCommRank(MPI_COMM_WORLD) == 0) {
    printf("\n ,world memory usage (KB): %3d, %6ld, %6ld\n", commSize, global_vmrss,
           global_vmsize);
  }
}

void print_memory_usage_comm(CommunicatorType comm) {
  unsigned long comm_vmrss, comm_vmsize;
  get_all_memory_usage_kb(&comm_vmrss, &comm_vmsize, comm);
  auto commSize = getCommSize(comm);
  if (getCommRank(comm) == 0) {
    printf(", comm memory usage (KB): %3d, %6ld, %6ld\n", commSize, comm_vmrss,
           comm_vmsize);
  }
}

void print_memory_usage_comm_minus_ref(unsigned long ref_vmrss, unsigned long ref_vmsize, CommunicatorType comm) {
  unsigned long comm_vmrss, comm_vmsize;
  get_all_memory_usage_kb(&comm_vmrss, &comm_vmsize, comm);
  auto commSize = getCommSize(comm);
  if (getCommRank(comm) == 0) {
    printf(", comm memory usage (KB): %3d, %6ld, %6ld\n", commSize, comm_vmrss-ref_vmrss, comm_vmsize-ref_vmsize);
  }
}

}  // namespace mpimemory
}  // namespace combigrid
