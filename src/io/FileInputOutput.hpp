#pragma once

// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <string>

#include "utils/Types.hpp"

namespace combigrid {
void writeConcatenatedFileRootOnly(const char* data, size_t sizeOfData, const std::string& path,
                                   MPI_Comm comm, bool replaceExistingFile = false);

bool getFileExistsRootOnly(const std::string& fileName, CommunicatorType comm, RankType rootRank = 0);

}  // namespace combigrid
