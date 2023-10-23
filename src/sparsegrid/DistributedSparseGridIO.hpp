#pragma once

#include <filesystem>
#include <fstream>
#include <string>

#include "io/MPIInputOutput.hpp"
#include "sparsegrid/DistributedSparseGridUniform.hpp"
#include "utils/Types.hpp"

namespace combigrid {

namespace DistributedSparseGridIO {

// helper functions to reduce a vector across a communicator
template <typename T>
void reduceVectorTowardsMe(std::vector<T>& vectorToReduce, MPI_Comm comm, MPI_Op operation) {
  MPI_Datatype dataType = getMPIDatatype(abstraction::getabstractionDataType<T>());
  MPI_Reduce(MPI_IN_PLACE, vectorToReduce.data(), vectorToReduce.size(), dataType, operation, 0,
             comm);
}

template <typename T>
void reduceVectorTowardsThem(std::vector<T>& vectorToReduce, MPI_Comm comm, MPI_Op operation) {
  MPI_Datatype dataType = getMPIDatatype(abstraction::getabstractionDataType<T>());
  MPI_Reduce(vectorToReduce.data(), MPI_IN_PLACE, vectorToReduce.size(), dataType, operation, 0,
             comm);
}

template <typename SparseGridType>
inline void writeMinMaxCoefficents(SparseGridType& dsg, const std::string& filename,
                                   size_t outputIndex, MPI_Comm firstCommToReduceAcross = MPI_COMM_SELF,
                                   MPI_Comm secondCommToReduceAcross = MPI_COMM_SELF) {
  if (dsg.getMaxCoefficientsPerSubspace().empty()) {
    dsg.accumulateMinMaxCoefficients();
  }

  // reduce the minimum and maximum values towards zero process
  // first, reduce within process groups
  if (getCommRank(firstCommToReduceAcross) == 0) {
    reduceVectorTowardsMe(dsg.getMinCoefficientsPerSubspace(), firstCommToReduceAcross,
                          MPI_MIN);
    reduceVectorTowardsMe(dsg.getMaxCoefficientsPerSubspace(), firstCommToReduceAcross,
                          MPI_MAX);
    //  then, reduce across process groups (only required on zeroth rank of each process group)
    if (getCommRank(secondCommToReduceAcross) == 0) {
      reduceVectorTowardsMe(dsg.getMinCoefficientsPerSubspace(),
                            secondCommToReduceAcross, MPI_MIN);
      reduceVectorTowardsMe(dsg.getMaxCoefficientsPerSubspace(),
                            secondCommToReduceAcross, MPI_MAX);

      // if on zero process, write them out to file
      std::ofstream ofs = std::ofstream(filename + "_" + std::to_string(outputIndex) + ".txt");
      for (typename SparseGridType::SubspaceIndexType i = 0;
          i <
          static_cast<typename SparseGridType::SubspaceIndexType>(dsg.getAllLevelVectors().size());
          ++i) {
        const auto& level = dsg.getLevelVector(i);
        if (dsg.getMinCoefficientsPerSubspace()[i] < std::numeric_limits<combigrid::real>::max())
          ofs << level << " : " << dsg.getMinCoefficientsPerSubspace()[i] << ", "
              << dsg.getMaxCoefficientsPerSubspace()[i] << std::endl;
      }
    } else {
      reduceVectorTowardsThem(dsg.getMinCoefficientsPerSubspace(),
                              secondCommToReduceAcross, MPI_MIN);
      reduceVectorTowardsThem(dsg.getMaxCoefficientsPerSubspace(),
                              secondCommToReduceAcross, MPI_MAX);
    }
  } else {
    reduceVectorTowardsThem(dsg.getMinCoefficientsPerSubspace(), firstCommToReduceAcross,
                            MPI_MIN);
    reduceVectorTowardsThem(dsg.getMaxCoefficientsPerSubspace(), firstCommToReduceAcross,
                            MPI_MAX);
  }
  dsg.clearMinMaxCoefficientsPerSubspace();
}

template <typename SparseGridType>
void writeToDiskChunked(const SparseGridType& dsg, const std::string& filePrefix) {
  std::string myFilename = filePrefix + std::to_string(dsg.getRank());
  std::ofstream ofp(myFilename, std::ios::out | std::ios::binary);
  ofp.write(reinterpret_cast<const char*>(dsg.getRawData()),
            dsg.getRawDataSize() * sizeof(typename SparseGridType::ElementType));
  ofp.close();
}

template <typename SparseGridType>
void readFromDiskChunked(SparseGridType& dsg, const std::string& filePrefix) {
  std::string myFilename = filePrefix + std::to_string(dsg.getRank());
  // assert that file is large enough
  [[maybe_unused]] size_t fileSize = std::filesystem::file_size(myFilename);
  assert(fileSize == dsg.getRawDataSize() * sizeof(typename SparseGridType::ElementType));
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
int writeSomeFiles(const SparseGridType& dsg, const std::string& fileName,
                   bool deleteExistingFile = false) {
  auto comm = theMPISystem()->getOutputComm();
  auto outputGroupSize = getCommSize(comm);
  auto filePart = theMPISystem()->getOutputGroupRank() / outputGroupSize;
  std::string filePartName = fileName + ".part" + std::to_string(filePart);

  MPI_Offset len = dsg.getRawDataSize();
  auto data = dsg.getRawData();
  int numWritten = mpiio::writeValuesConsecutive<typename SparseGridType::ElementType>(
      data, len, filePartName, comm, deleteExistingFile);
  return numWritten;
}

template <typename SparseGridType>
int readSomeFiles(SparseGridType& dsg, const std::string& fileName) {
  auto comm = theMPISystem()->getOutputComm();
  auto outputGroupSize = getCommSize(comm);
  auto filePart = theMPISystem()->getOutputGroupRank() / outputGroupSize;
  std::string filePartName = fileName + ".part" + std::to_string(filePart);

  // get offset in file
  MPI_Offset len = dsg.getRawDataSize();
  auto data = dsg.getRawData();
  int numRead = mpiio::readValuesConsecutive<typename SparseGridType::ElementType>(
      data, len, filePartName, comm);
  return numRead;
}

template <typename SparseGridType>
int readSomeFilesAndReduce(SparseGridType& dsg, const std::string& fileName,
                           uint32_t maxMiBToReadPerThread) {
  auto comm = theMPISystem()->getOutputComm();
  auto outputGroupSize = getCommSize(comm);
  auto filePart = theMPISystem()->getOutputGroupRank() / outputGroupSize;
  std::string filePartName = fileName + ".part" + std::to_string(filePart);

  const int numElementsToBuffer =
      CombiCom::getGlobalReduceChunkSize<typename SparseGridType::ElementType>(
          maxMiBToReadPerThread);

  // get offset in file
  const MPI_Offset len = dsg.getRawDataSize();
  auto data = dsg.getRawData();
  int numReduced = mpiio::readReduceValuesConsecutive<typename SparseGridType::ElementType>(
      data, len, filePartName, comm, numElementsToBuffer,
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

template <typename SparseGridType, typename ReduceFunctionType>
int readReduceSubspaceSizesFromFiles(SparseGridType& dsg, const std::vector<std::string>& fileNames,
                                     ReduceFunctionType reduceFunction, int numElementsToBuffer = 0,
                                     bool withCollectiveBuffering = false) {
  auto comm = dsg.getCommunicator();
  MPI_Offset len = dsg.getNumSubspaces();
  if (numElementsToBuffer == 0) {
    numElementsToBuffer = len;
  }

  int numReduced = mpiio::readMultipleReduceValuesConsecutive<SubspaceSizeType>(
      dsg.getSubspaceDataSizes().data(), len, fileNames, comm, numElementsToBuffer, reduceFunction,
      withCollectiveBuffering);

  return numReduced;
}
}  // namespace DistributedSparseGridIO
}  // namespace combigrid