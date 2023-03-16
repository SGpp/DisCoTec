#pragma once

// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <algorithm>
#include <numeric>

#include "utils/Types.hpp"

namespace combigrid {
namespace mpiio {

static MPI_Info getNewConsecutiveMpiInfo() {
  // see: https://wickie.hlrs.de/platforms/index.php/MPI-IO
  // to be modified externally e.g. via romio hints
  MPI_Info info = MPI_INFO_NULL;
  MPI_Info_create(&info);
  // will not overlap read and write operations
  MPI_Info_set(info, "romio_no_indep_rw", "true");
  // disable caching of file contents in kernel
  MPI_Info_set(info, "direct_io", "true");
  MPI_Info_set(info, "direct_read", "true");
  MPI_Info_set(info, "direct_write", "true");
  // disable ROMIO's data-sieving
  MPI_Info_set(info, "romio_ds_read", "disable");
  MPI_Info_set(info, "romio_ds_write", "disable");
  // disable ROMIO's collective buffering //TODO test
  MPI_Info_set(info, "romio_cb_write", "disable");
  MPI_Info_set(info, "romio_cb_read", "disable");
  return info;
}

template <typename T>
bool writeValuesConsecutive(const T* valuesStart, MPI_Offset numValues, const std::string& fileName,
                            combigrid::CommunicatorType comm) {
  // get offset in file
  MPI_Offset pos = 0;
  MPI_Exscan(&numValues, &pos, 1, MPI_OFFSET, MPI_SUM, comm);

  MPI_Info info = getNewConsecutiveMpiInfo();

  // open file
  MPI_File fh;
  int err = MPI_File_open(comm, fileName.c_str(), MPI_MODE_CREATE | MPI_MODE_EXCL | MPI_MODE_WRONLY,
                          info, &fh);
  if (err == MPI_ERR_FILE_EXISTS) {
    // file already existed, delete it and create new file
    int mpi_rank;
    MPI_Comm_rank(comm, &mpi_rank);
    if (mpi_rank == 0) {
      MPI_File_delete(fileName.c_str(), MPI_INFO_NULL);
    }
    MPI_Barrier(comm);
    err = MPI_File_open(comm, fileName.c_str(), MPI_MODE_CREATE | MPI_MODE_EXCL | MPI_MODE_WRONLY,
                        info, &fh);
  }
  assert(err == MPI_SUCCESS);

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
  MPI_Info_free(&info);
  return err == MPI_SUCCESS;
}

template <typename T>
bool readValuesConsecutive(T* valuesStart, MPI_Offset numValues, const std::string& fileName,
                           combigrid::CommunicatorType comm) {
  MPI_Offset pos = 0;
  MPI_Exscan(&numValues, &pos, 1, MPI_OFFSET, MPI_SUM, comm);

  MPI_Info info = getNewConsecutiveMpiInfo();

  // open file
  MPI_File fh;
  int err = MPI_File_open(comm, fileName.c_str(), MPI_MODE_RDONLY, info, &fh);
  if (err != MPI_SUCCESS) {
    std::cerr << err << " while reading OneFileFromDisk " << fileName << std::endl;
    throw std::runtime_error("read: could not open!");
  }
#ifndef NDEBUG
  MPI_Offset fileSize = 0;
  MPI_File_get_size(fh, &fileSize);
  if (fileSize < numValues * sizeof(T)) {
    // loud failure if file is too small
    std::cerr << fileSize << " and not " << numValues << std::endl;
    throw std::runtime_error("read: file size too small!");
  }
#endif

  // read from single file with MPI-IO
  MPI_Datatype dataType = getMPIDatatype(abstraction::getabstractionDataType<T>());
  MPI_Status status;
  err = MPI_File_read_at_all(fh, pos * sizeof(T), valuesStart, static_cast<int>(numValues),
                             dataType, &status);
  MPI_File_close(&fh);
  MPI_Info_free(&info);
  if (err != MPI_SUCCESS) {
    // non-failure
    std::cerr << err << " in MPI_File_read_at_all" << std::endl;
    return false;
  }

#ifndef NDEBUG
  int readcount = 0;
  MPI_Get_count(&status, dataType, &readcount);
  if (readcount < numValues) {
    // loud failure
    std::cerr << "read " << readcount << " and not " << numValues << std::endl;
    throw std::runtime_error("read: not enough data read!");
  }
#endif

  return true;
}

template <typename T, typename ReduceFunctionType>
bool readReduceValuesConsecutive(T* valuesStart, MPI_Offset numValues, const std::string& fileName,
                                 combigrid::CommunicatorType comm, int numElementsToBuffer,
                                 ReduceFunctionType reduceFunction) {
  MPI_Offset pos = 0;
  MPI_Exscan(&numValues, &pos, 1, MPI_OFFSET, MPI_SUM, comm);
  MPI_Info info = getNewConsecutiveMpiInfo();

  // open file
  MPI_File fh;
  int err = MPI_File_open(comm, fileName.c_str(), MPI_MODE_RDONLY, info, &fh);
  if (err != MPI_SUCCESS) {
    // silent failure
    std::cerr << err << " while reducing OneFileFromDisk " << fileName << std::endl;
    throw std::runtime_error("read: could not open!");
  }

  // read from single file with MPI-IO
  MPI_Datatype dataType = getMPIDatatype(abstraction::getabstractionDataType<T>());
  MPI_Status status;
  int readcount = 0;
  std::vector<T> buffer(numElementsToBuffer);
  auto writePointer = valuesStart;
  while (readcount < numValues) {
    // if there is less to read than the buffer length, overwrite numElementsToBuffer
    if (numValues - readcount < numElementsToBuffer) {
      numElementsToBuffer = numValues - readcount;
      buffer.resize(numElementsToBuffer);
    }
    err = MPI_File_read_at(fh, pos * sizeof(T), buffer.data(),
                           static_cast<int>(numElementsToBuffer), dataType, &status);
    int readcountIncrement = 0;
    MPI_Get_count(&status, dataType, &readcountIncrement);
    assert(readcountIncrement > 0);
    // reduce with present sparse grid data
    std::transform(buffer.cbegin(), buffer.cend(), writePointer, writePointer, reduceFunction);
    readcount += readcountIncrement;
    pos += readcountIncrement;
    std::advance(writePointer, readcountIncrement);
  }
  MPI_File_close(&fh);
  MPI_Info_free(&info);
  return true;
}
}  // namespace mpiio
}  // namespace combigrid