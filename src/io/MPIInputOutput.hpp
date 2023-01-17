#pragma once

// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <numeric>

#include "utils/Types.hpp"

namespace combigrid {
namespace mpiio {

template <typename T>
bool writeValuesConsecutive(const T* valuesStart, MPI_Offset numValues, const std::string& fileName,
                            combigrid::CommunicatorType comm) {
  // get offset in file
  MPI_Offset pos = 0;
  MPI_Exscan(&numValues, &pos, 1, getMPIDatatype(abstraction::getabstractionDataType<MPI_Offset>()),
             MPI_SUM, comm);

  // see: https://wickie.hlrs.de/platforms/index.php/MPI-IO
  // better if set externally
  MPI_Info info = MPI_INFO_NULL;

  // open file
  MPI_File fh;
  int err = MPI_File_open(comm, fileName.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &fh);
  if (err != MPI_SUCCESS) {
    // file already existed, delete it and create new file
    int mpi_rank;
    MPI_Comm_rank(comm, &mpi_rank);
    if (mpi_rank == 0) {
      MPI_File_delete(fileName.c_str(), MPI_INFO_NULL);
    }
    err = MPI_File_open(comm, fileName.c_str(), MPI_MODE_CREATE | MPI_MODE_EXCL | MPI_MODE_WRONLY,
                        info, &fh);
  }

  if (err == MPI_SUCCESS) {
    // write to single file with MPI-IO
    MPI_Datatype dataType = getMPIDatatype(abstraction::getabstractionDataType<T>());
    MPI_Status status;
    err = MPI_File_write_at_all(fh, pos * sizeof(T), valuesStart, static_cast<int>(numValues),
                                dataType, &status);
    if (err != MPI_SUCCESS) {
      std::cerr << err << " in MPI_File_write_at_all" << std::endl;
    }
#ifndef NDEBUG
    int numWritten = 0;
    MPI_Get_count(&status, dataType, &numWritten);
    if (numWritten != numValues) {
      std::cout << "not written enough: " << numWritten << " instead of " << numValues << std::endl;
      err = ~MPI_SUCCESS;
    }
#endif  // !NDEBUG
  }

  MPI_File_close(&fh);
  return err == MPI_SUCCESS;
}
}  // namespace mpiio
}  // namespace combigrid