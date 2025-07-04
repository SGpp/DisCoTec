#pragma once

// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <algorithm>
#include <numeric>

#include "../utils/Types.hpp"

namespace combigrid {
static inline std::string getMpiErrorString(int err) {
  int len;
  char errString[MPI_MAX_ERROR_STRING];
  MPI_Error_string(err, errString, &len);
  return std::string(errString);
}
namespace mpiio {

static MPI_Info getNewConsecutiveMpiInfo(bool withCollectiveBuffering) {
  // see: https://wickie.hlrs.de/platforms/index.php/MPI-IO
  // to be further modified externally e.g. via romio hints
  MPI_Info info = MPI_INFO_NULL;
  MPI_Info_create(&info);

  // always disable ROMIO's data-sieving
  MPI_Info_set(info, "romio_ds_read", "disable");
  MPI_Info_set(info, "romio_ds_write", "disable");

  if (withCollectiveBuffering) {
    // enable ROMIO's collective buffering
    MPI_Info_set(info, "collective_buffering", "true");
    MPI_Info_set(info, "romio_no_indep_rw", "true");
    MPI_Info_set(info, "romio_cb_write", "enable");
    MPI_Info_set(info, "romio_cb_read", "enable");
  } else {
    // disable ROMIO's collective buffering
    MPI_Info_set(info, "collective_buffering", "false");
    MPI_Info_set(info, "romio_no_indep_rw", "false");
    MPI_Info_set(info, "romio_cb_write", "disable");
    MPI_Info_set(info, "romio_cb_read", "disable");
  }
  return info;
}

inline MPI_Offset getPositionFromNumValues(MPI_Offset numValues, combigrid::CommunicatorType comm) {
  MPI_Offset position = 0;
  MPI_Exscan(&numValues, &position, 1, MPI_OFFSET, MPI_SUM, comm);
  return position;
}

template <typename T>
class [[nodiscard]] MPIFileConsecutive {
 public:
  // move but no copy
  MPIFileConsecutive(const MPIFileConsecutive&) = delete;
  MPIFileConsecutive& operator=(const MPIFileConsecutive&) = delete;
  MPIFileConsecutive(MPIFileConsecutive&& other) { *this = std::move(other); }
  MPIFileConsecutive& operator=(MPIFileConsecutive&& other) {
    file_ = other.file_;
    other.file_ = MPI_FILE_NULL;
    return *this;
  };

  ~MPIFileConsecutive() { MPI_File_close(&file_); }
  static MPIFileConsecutive getFileToWrite(const std::string& fileName,
                                           combigrid::CommunicatorType comm,
                                           bool replaceExistingFile, bool withCollectiveBuffering,
                                           bool writeOnce = true) {
    MPI_Info info = getNewConsecutiveMpiInfo(withCollectiveBuffering);
    if (writeOnce) {
      MPI_Info_set(info, "access_style", "write_once");
    }
    int commSize;
    MPI_Comm_size(comm, &commSize);
    std::string commSizeStr = std::to_string(commSize);
    MPI_Info_set(info, "nb_procs", commSizeStr.c_str());
    MPI_File file;

    int err = MPI_File_open(comm, fileName.c_str(),
                            MPI_MODE_CREATE | MPI_MODE_EXCL | MPI_MODE_WRONLY, info, &file);
    if (err != MPI_SUCCESS) {
      auto openErrorString = getMpiErrorString(err);
      if (err != MPI_ERR_FILE_EXISTS && openErrorString.find("File exists") == std::string::npos) {
        // there are some weird MPI-IO implementations that define a new error code for "file
        // exists"...
        std::cerr << "potential write/open error: " << std::to_string(err) << " " << openErrorString
                  << std::endl;
      }
      if (replaceExistingFile) {
        // file already existed, delete it and create new file
        int mpi_rank;
        MPI_Comm_rank(comm, &mpi_rank);
        if (mpi_rank == 0) {
          MPI_File_delete(fileName.c_str(), MPI_INFO_NULL);
        }
        MPI_Barrier(comm);
        err = MPI_File_open(comm, fileName.c_str(),
                            MPI_MODE_CREATE | MPI_MODE_EXCL | MPI_MODE_WRONLY, info, &file);
      } else {
        // open file without creating it
        err = MPI_File_open(comm, fileName.c_str(), MPI_MODE_WRONLY, info, &file);
      }
    }
    if (err != MPI_SUCCESS) {
      std::cerr << "Open error " << fileName << " :" << std::to_string(err) << " "
                << getMpiErrorString(err) << std::endl;
    }
    MPI_Info_free(&info);
    return MPIFileConsecutive(file);
  }

  static MPIFileConsecutive getFileToRead(const std::string& fileName,
                                          combigrid::CommunicatorType comm,
                                          bool withCollectiveBuffering, bool readOnce) {
    MPI_Info info = getNewConsecutiveMpiInfo(withCollectiveBuffering);
    if (readOnce) {
      MPI_Info_set(info, "access_style", "read_once");
    }

    MPI_File file;
    int err = MPI_File_open(comm, fileName.c_str(), MPI_MODE_RDONLY, info, &file);
    if (err != MPI_SUCCESS) {
      std::cerr << err << " while reading OneFileFromDisk " << fileName << std::endl;
      throw std::runtime_error("read: could not open! " + fileName + ": " + getMpiErrorString(err));
    }
    return MPIFileConsecutive(file);
  }

  int writeValuesToFileAtPosition(const T* valuesStart, MPI_Offset numValues, MPI_Offset position) {
    // write to single file with MPI-IO
    MPI_Datatype dataType = getMPIDatatype(abstraction::getabstractionDataType<T>());
    MPI_Status status;
    int err = MPI_File_write_at_all(file_, position * sizeof(T), valuesStart,
                                    static_cast<int>(numValues), dataType, &status);
    if (err != MPI_SUCCESS) {
      std::cerr << err << " in MPI_File_write_at_all" << std::endl;
      std::cerr << getMpiErrorString(err) << std::endl;
    }
#ifndef NDEBUG
    int numWritten = 0;
    MPI_Get_count(&status, dataType, &numWritten);
    if (numWritten != numValues) {
      std::cout << "not written enough: " << numWritten << " instead of " << numValues << std::endl;
      err = ~MPI_SUCCESS;
    }
#endif  // !NDEBUG
    return (err == MPI_SUCCESS) ? static_cast<int>(numValues) : 0;
  }

  int readValuesFromFileAtPosition(T* valuesStart, MPI_Offset numValues, MPI_Offset position) {
    // read from single file with MPI-IO
    MPI_Datatype dataType = getMPIDatatype(abstraction::getabstractionDataType<T>());
    MPI_Status status;
    auto err = MPI_File_read_at_all(file_, position * sizeof(T), valuesStart,
                                    static_cast<int>(numValues), dataType, &status);
    if (err != MPI_SUCCESS) {
      // non-failure
      std::cerr << err << " in MPI_File_read_at_all" << std::endl;
      return 0;
    }
    int readcountIncrement = 0;
    MPI_Get_count(&status, dataType, &readcountIncrement);
    if (readcountIncrement < numValues) {
      // non-failure
      std::cerr << "read " << readcountIncrement << " and not " << numValues << std::endl;
    }
    return readcountIncrement;
  }

  int readValuesFromFileAtPositionSingleRank(T* valuesStart, MPI_Offset numValues,
                                             MPI_Offset position) const {
    MPI_Datatype dataType = getMPIDatatype(abstraction::getabstractionDataType<T>());
    MPI_Status status;
    auto err = MPI_File_read_at(file_, position * sizeof(T), valuesStart,
                                static_cast<int>(numValues), dataType, &status);
    if (err != MPI_SUCCESS) {
      // non-failure
      std::cerr << err << " in MPI_File_read_at_all" << std::endl;
      return 0;
    }
    int readcountIncrement = 0;
    MPI_Get_count(&status, dataType, &readcountIncrement);
    if (readcountIncrement < numValues) {
      // non-failure
      std::cerr << "read " << readcountIncrement << " and not " << numValues << std::endl;
    }
    return readcountIncrement;
  }

  bool checkFileSizeConsecutive(MPI_Offset myNumValues, CommunicatorType comm) const {
    // get total number of values
    MPI_Offset totalNumValues;
    MPI_Allreduce(&myNumValues, &totalNumValues, 1, MPI_OFFSET, MPI_SUM, comm);

    // if (rankInComm == 0) { // TODO is it important that only rank 0 checks the file size?
    // get file size
    MPI_Offset fileSize;
    MPI_File_get_size(file_, &fileSize);

    // check whether file size is correct
    if (fileSize != static_cast<MPI_Offset>(totalNumValues * sizeof(T))) {
      throw std::runtime_error("file size does not match number of values; should be " +
                               std::to_string(totalNumValues * sizeof(T)) + " but is " +
                               std::to_string(fileSize) + " bytes");
    }
    return true;
  }

 private:
  explicit MPIFileConsecutive(MPI_File file) : file_(file) {}
  MPI_File file_;
};

template <typename T>
int writeValuesConsecutive(const T* valuesStart, MPI_Offset numValues, const std::string& fileName,
                           combigrid::CommunicatorType comm, bool replaceExistingFile = false,
                           bool withCollectiveBuffering = false) {
  auto file = MPIFileConsecutive<T>::getFileToWrite(fileName, comm, replaceExistingFile,
                                                    withCollectiveBuffering);
  auto position = getPositionFromNumValues(numValues, comm);
  return file.writeValuesToFileAtPosition(valuesStart, numValues, position);
}

template <typename T>
int readValuesConsecutive(T* valuesStart, MPI_Offset numValues, const std::string& fileName,
                          combigrid::CommunicatorType comm, bool withCollectiveBuffering = false) {
  auto file = MPIFileConsecutive<T>::getFileToRead(fileName, comm, withCollectiveBuffering, true);
  file.checkFileSizeConsecutive(numValues, comm);

  auto position = getPositionFromNumValues(numValues, comm);
  return file.readValuesFromFileAtPosition(valuesStart, numValues, position);
}

template <typename T, typename ReduceFunctionType>
int readReduceValuesConsecutive(T* valuesStart, MPI_Offset numValues, const std::string& fileName,
                                combigrid::CommunicatorType comm, int numElementsToBuffer,
                                ReduceFunctionType reduceFunction,
                                bool withCollectiveBuffering = false) {
  auto file = MPIFileConsecutive<T>::getFileToRead(fileName, comm, withCollectiveBuffering, false);
  file.checkFileSizeConsecutive(numValues, comm);

  MPI_Offset position = getPositionFromNumValues(numValues, comm);
  int readcount = 0;
  std::vector<T> buffer(numElementsToBuffer);
  auto writePointer = valuesStart;
  while (readcount < numValues) {
    // if there is less to read than the buffer length, overwrite numElementsToBuffer
    if (numValues - readcount < numElementsToBuffer) {
      numElementsToBuffer = static_cast<int>(numValues - readcount);
      buffer.resize(numElementsToBuffer);
    }
    auto readcountIncrement =
        file.readValuesFromFileAtPosition(buffer.data(), numElementsToBuffer, position);
    assert(readcountIncrement > 0);
    // reduce with present sparse grid data
    std::transform(buffer.cbegin(), buffer.cend(), writePointer, writePointer, reduceFunction);
    readcount += readcountIncrement;
    position += readcountIncrement;
    std::advance(writePointer, readcountIncrement);
  }
  return readcount;
}

// n-ary std::transform https://stackoverflow.com/a/21768697
template <class Functor, class OutputIterator, class Input1, class... Inputs>
OutputIterator transformN(Functor f, OutputIterator out, Input1 first1, Input1 last1,
                          Inputs... inputs) {
  while (first1 != last1) *out++ = f(*first1++, *inputs++...);
  return out;
}

template <typename T, typename ReduceFunctionType>
int readMultipleReduceValuesConsecutive(T* valuesStart, MPI_Offset numValues,
                                        const std::vector<std::string>& fileNames,
                                        combigrid::CommunicatorType comm, int numElementsToBuffer,
                                        ReduceFunctionType reduceFunction,
                                        bool withCollectiveBuffering = false) {
  // open multiple files
  std::vector<MPIFileConsecutive<T>> files;
  files.reserve(fileNames.size());
  for (size_t i = 0; i < fileNames.size(); ++i) {
    files.push_back(std::move(
        MPIFileConsecutive<T>::getFileToRead(fileNames[i], comm, withCollectiveBuffering, false)));
    files.back().checkFileSizeConsecutive(numValues, comm);
  }

  // read from files with MPI-IO
  int readcount = 0;
  MPI_Offset position = getPositionFromNumValues(numValues, comm);
  std::vector<std::vector<T>> buffers(fileNames.size(), std::vector<T>(numElementsToBuffer));
  auto writePointer = valuesStart;
  while (readcount < numValues) {
    // if there is less to read than the buffer length, overwrite numElementsToBuffer
    if (numValues - readcount < numElementsToBuffer) {
      numElementsToBuffer = static_cast<int>(numValues - readcount);
      for (auto& buffer : buffers) {
        buffer.resize(numElementsToBuffer);
      }
    }
    int readcountIncrement = 0;
    for (size_t i = 0; i < buffers.size(); ++i) {
      readcountIncrement =
          files[i].readValuesFromFileAtPosition(buffers[i].data(), numElementsToBuffer, position);
      assert(readcountIncrement > 0);
    }
    // reduce with present sparse grid data
    // todo is there a variadic way of doing this?
    transformN(reduceFunction, writePointer, buffers[0].cbegin(), buffers[0].cend(), writePointer,
               buffers[1].cbegin());
    readcount += readcountIncrement;
    position += readcountIncrement;
    std::advance(writePointer, readcountIncrement);
  }

  return readcount;
}

}  // namespace mpiio
}  // namespace combigrid
