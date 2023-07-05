#pragma once

#include <string>

#include "io/MPIInputOutput.hpp"
#include "sparsegrid/DistributedSparseGridUniform.hpp"
#include "utils/Types.hpp"

namespace combigrid {

namespace DistributedSparseGridIO {
template <typename SparseGridType>
inline void writeMinMaxCoefficents(const SparseGridType& dsg, const std::string& filename,
                                   size_t outputIndex) {
  assert(dsg.isSubspaceDataCreated());
  bool writerProcess = false;
  std::ofstream ofs;
  if (dsg.getRank() == 0) {
    writerProcess = true;
    ofs = std::ofstream(filename + "_" + std::to_string(outputIndex) + ".txt");
  }
  // iterate subspaces
  assert(dsg.getAllLevelVectors().size() > 0);
  MPI_Datatype dataType = getMPIDatatype(abstraction::getabstractionDataType<combigrid::real>());
  auto realmax = std::numeric_limits<combigrid::real>::max();
  auto realmin = std::numeric_limits<combigrid::real>::min();
  auto smaller_real = [](const typename SparseGridType::ElementType& one,
                         const typename SparseGridType::ElementType& two) {
    return std::real(one) < std::real(two);
  };

  for (typename SparseGridType::SubspaceIndexType i = 0;
       i < static_cast<typename SparseGridType::SubspaceIndexType>(dsg.getAllLevelVectors().size());
       ++i) {
    auto minimumValue = realmax;
    auto maximumValue = realmin;
    if (dsg.getSubspaceDataSizes()[i] > 0) {
      // auto first = subspacesData_.begin();
      auto first = dsg.getData(i);
      auto last = first + dsg.getSubspaceDataSizes()[i];
      auto it = std::min_element(first, last, smaller_real);
      minimumValue = std::real(*it);
      first = dsg.getData(i);
      it = std::max_element(first, last, smaller_real);
      maximumValue = std::real(*it);
    }
    // allreduce the minimum and maximum values
    MPI_Allreduce(MPI_IN_PLACE, &minimumValue, 1, dataType, MPI_MIN, dsg.getCommunicator());
    MPI_Allreduce(MPI_IN_PLACE, &maximumValue, 1, dataType, MPI_MAX, dsg.getCommunicator());

    // if on zero process, write them out to file
    if (writerProcess) {
      const auto& level = dsg.getLevelVector(i);
      if (minimumValue < realmax)
        ofs << level << " : " << minimumValue << ", " << maximumValue << std::endl;
    }
  }
}

template <typename SparseGridType>
void writeToDiskChunked(const SparseGridType& dsg, std::string filePrefix) {
  std::string myFilename = filePrefix + std::to_string(dsg.getRank());
  std::ofstream ofp(myFilename, std::ios::out | std::ios::binary);
  ofp.write(reinterpret_cast<const char*>(dsg.getRawData()),
            dsg.getRawDataSize() * sizeof(typename SparseGridType::ElementType));
  ofp.close();
}

template <typename SparseGridType>
void readFromDiskChunked(SparseGridType& dsg, std::string filePrefix) {
  std::string myFilename = filePrefix + std::to_string(dsg.getRank());
  std::ifstream ifp(myFilename, std::ios::in | std::ios::binary);
  ifp.read(reinterpret_cast<char*>(dsg.getRawData()),
           dsg.getRawDataSize() * sizeof(typename SparseGridType::ElementType));
  ifp.close();
}

template <typename SparseGridType>
int writeOneFile(const SparseGridType& dsg, const std::string& fileName,
                 bool deleteExistingFile = false) {
  auto comm = dsg.getCommunicator();

  MPI_Offset len = dsg.getRawDataSize();
  auto data = dsg.getRawData();
  int numWritten = mpiio::writeValuesConsecutive<typename SparseGridType::ElementType>(
      data, len, fileName, comm, deleteExistingFile);
  return numWritten;
}

template <typename SparseGridType>
int readOneFile(SparseGridType& dsg, const std::string& fileName) {
  auto comm = dsg.getCommunicator();

  // get offset in file
  MPI_Offset len = dsg.getRawDataSize();
  auto data = dsg.getRawData();
  int numRead =
      mpiio::readValuesConsecutive<typename SparseGridType::ElementType>(data, len, fileName, comm);
  return numRead;
}

template <typename SparseGridType>
int readOneFileAndReduce(SparseGridType& dsg, const std::string& fileName,
                         uint32_t maxMiBToReadPerThread) {
  auto comm = dsg.getCommunicator();

  const int numElementsToBuffer =
      CombiCom::getGlobalReduceChunkSize<typename SparseGridType::ElementType>(
          maxMiBToReadPerThread);

  // get offset in file
  const MPI_Offset len = dsg.getRawDataSize();
  auto data = dsg.getRawData();
  int numReduced = mpiio::readReduceValuesConsecutive<typename SparseGridType::ElementType>(
      data, len, fileName, comm, numElementsToBuffer,
      std::plus<typename SparseGridType::ElementType>{});

  return numReduced;
}

template <typename SparseGridType>
int writeSubspaceSizesToFile(const SparseGridType& dsg, const std::string& fileName) {
  auto comm = dsg.getCommunicator();
  MPI_Offset len = dsg.getNumSubspaces();
  int numWritten = mpiio::writeValuesConsecutive<SubspaceSizeType>(
      dsg.getSubspaceDataSizes().data(), len, fileName, comm);
  return numWritten;
}

template <typename SparseGridType>
int readSubspaceSizesFromFile(SparseGridType& dsg, const std::string& fileName,
                              bool withCollectiveBuffering = false) {
  auto comm = dsg.getCommunicator();
  MPI_Offset len = dsg.getNumSubspaces();
  int numRead = mpiio::readValuesConsecutive<SubspaceSizeType>(
      dsg.getSubspaceDataSizes().data(), len, fileName, comm, withCollectiveBuffering);
  return numRead;
}

template <typename SparseGridType, typename ReduceFunctionType>
int readReduceSubspaceSizesFromFile(SparseGridType& dsg, const std::string& fileName,
                                    ReduceFunctionType reduceFunction, int numElementsToBuffer = 0,
                                    bool withCollectiveBuffering = false) {
  auto comm = dsg.getCommunicator();
  MPI_Offset len = dsg.getNumSubspaces();
  if (numElementsToBuffer == 0) {
    numElementsToBuffer = len;
  }

  int numReduced = mpiio::readReduceValuesConsecutive<SubspaceSizeType>(
      dsg.getSubspaceDataSizes().data(), len, fileName, comm, numElementsToBuffer, reduceFunction,
      withCollectiveBuffering);

  return numReduced;
}
}  // namespace DistributedSparseGridIO
}  // namespace combigrid