#ifndef SRC_SGPP_COMBIGRID_MPI_MEMORY_HPP_
#define SRC_SGPP_COMBIGRID_MPI_MEMORY_HPP_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include "sgpp/distributedcombigrid/utils/Types.hpp"

/**
 * @brief utility to get memory usage in the MPI processes. slow!
 * combined from the examples at
 * https://hpcf.umbc.edu/general-productivity/checking-memory-usage/
 */

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
int get_memory_usage_kb(unsigned long* vmrss_kb, unsigned long* vmsize_kb);

int get_all_memory_usage_kb(unsigned long* vmrss, unsigned long* vmsize, int np,
                            CommunicatorType comm);

int get_memory_usage_local_kb(unsigned long* local_vmrss, unsigned long* local_vmsize);

void print_memory_usage_local();

void print_memory_usage_world();
}  // namespace mpimemory
}  // namespace combigrid
#endif /* SRC_SGPP_COMBIGRID_MPI_MEMORY_HPP_ */
