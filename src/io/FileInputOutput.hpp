#pragma once

// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <string>

namespace combigrid {
void writeConcatenatedFileRootOnly(const char* data, int sizeOfData, const std::string& path,
                                   MPI_Comm comm, bool replaceExistingFile = false);

}  // namespace combigrid