#ifndef COMBICOM_HPP_
#define COMBICOM_HPP_

#include "fullgrid/DistributedFullGrid.hpp"
#include "fullgrid/FullGrid.hpp"
#include "mpi/MPISystem.hpp"
#include "mpi/MPITags.hpp"
#include "mpi/OpenMPUtils.hpp"
#include "sparsegrid/DistributedSparseGridUniform.hpp"

namespace combigrid {

namespace CombiCom {

/**
 * @brief check if the allocated subspace sizes of a DistributedSparseGrid are correct
 *
 * They must be either 0 or the same as the declared size of the subspace.
 * Only fails if the sizes are different and the program was compiled in Debug mode;
 * otherwise only prints to std::cout.
 *
 * @tparam SparseGridType (derived from) DistributedSparseGridUniform
 * @param dsg the DistributedSparseGrid
 * @param subspaceSizes the declared sizes of the allocated subspaces
 * @return the sum of all elements in \p subspaceSizes
 */
template <typename SparseGridType>
size_t sumAndCheckSubspaceSizes(const SparseGridType& dsg,
                                const std::vector<SubspaceSizeType>& subspaceSizes) {
  size_t bsize = 0;
  for (size_t i = 0; i < subspaceSizes.size(); ++i) {
    // check for implementation errors, the reduced subspace size should not be
    // different from the size of already initialized subspaces
    bool check = (subspaceSizes[i] == 0 || dsg.getDataSize(i) == 0 ||
                  subspaceSizes[i] == dsg.getDataSize(i));

    if (!check) {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      std::cout << "sumAndCheckSubspaceSizes; l = " << dsg.getLevelVector(i) << " "
                << "rank = " << rank << " "
                << "ssize = " << subspaceSizes[i] << " "
                << "dsize = " << dsg.getDataSize(i) << std::endl;
      assert(false);
    }

    bsize += subspaceSizes[i];
  }
  return bsize;
}

/**
 * @brief check if all allocated subspace sizes of a DistributedSparseGrid are correct
 *
 * like sumAndCheckSubspaceSizes for all subspace data sizes of the DistributedSparseGrid
 *
 * @tparam SparseGridType (derived from) DistributedSparseGridUniform
 * @param dsg the DistributedSparseGrid
 * @return true if the sum of all elements in the subspace sizes is equal to the raw data size
 */
template <typename SparseGridType>
bool sumAndCheckSubspaceSizes(const SparseGridType& dsg) {
  auto bsize = sumAndCheckSubspaceSizes(dsg, dsg.getSubspaceDataSizes());
  return bsize == dsg.getRawDataSize();
}

/**
 * @brief convert Mebibytes to Bytes
 *
 * @param numberOfMiB the number of Mebibytes
 * @return the number of Bytes
 */
static size_t MiBtoBytes(uint32_t numberOfMiB) {
  return powerOfTwoByBitshift(20) * static_cast<size_t>(numberOfMiB);
}

/**
 * @brief get the chunk size for global reduction, in number of elements, depending on the number of
 * threads
 *
 * @tparam SG_ELEMENT the type of the elements in the sparse grid
 * @param maxMiBToSendPerThread the maximum number of MiB to send per OpenMP thread (~= per core)
 * @return the chunk size in number of elements
 */
template <typename SG_ELEMENT>
static int getGlobalReduceChunkSize(uint32_t maxMiBToSendPerThread) {
  if (maxMiBToSendPerThread > static_cast<uint32_t>(std::numeric_limits<int>::max()) ||
      maxMiBToSendPerThread == 0) {
    return std::numeric_limits<int>::max();
  }
  assert(maxMiBToSendPerThread > 0);
  auto maxBytesToSend = CombiCom::MiBtoBytes(maxMiBToSendPerThread);
  // allreduce up to 16MiB at a time (when using double precision)
  auto numOMPThreads = OpenMPUtils::getNumThreads();
  int chunkSize = static_cast<int>(maxBytesToSend / sizeof(SG_ELEMENT) * numOMPThreads);
  assert(chunkSize > 0);
  return chunkSize;
}

/**
 * @brief perform distributed global sparse grid reduce
 *
 * sparse grid reduce: entire sparse grid is used for reduction
 * global = between process groups; this assumes that each process group has already
 * performed its local reduction INTO the sparse grid
 *
 * @tparam SparseGridType (derived from) DistributedSparseGridUniform
 * @param dsg the DistributedSparseGrid
 * @param maxMiBToSendPerThread the maximum number of MiB to send per OpenMP thread (~= per core) at
 * once
 * @param globalReduceRankThatCollects the rank that collects the data, or MPI_PROC_NULL if
 * allreduce
 * @param globalComm the communicator for the global reduction
 */
template <typename SparseGridType>
void distributedGlobalSparseGridReduce(
    SparseGridType& dsg, uint32_t maxMiBToSendPerThread,
    RankType globalReduceRankThatCollects = MPI_PROC_NULL,
    MPI_Comm globalComm = theMPISystem()->getGlobalReduceComm()) {
  assert(globalComm != MPI_COMM_NULL);
  assert(dsg.isSubspaceDataCreated() && "Only perform reduce with allocated data");

  typename SparseGridType::ElementType* subspacesData = dsg.getRawData();
  size_t subspacesDataSize = dsg.getRawDataSize();
  assert(subspacesDataSize == dsg.getAccumulatedDataSize());

  // global reduce
  MPI_Datatype dtype = abstraction::getMPIDatatype(
      abstraction::getabstractionDataType<typename SparseGridType::ElementType>());

  auto chunkSize =
      getGlobalReduceChunkSize<typename SparseGridType::ElementType>(maxMiBToSendPerThread);
  size_t sentRecvd = 0;
  if (globalReduceRankThatCollects == MPI_PROC_NULL) {
    while ((subspacesDataSize - sentRecvd) / chunkSize > 0) {
      MPI_Allreduce(MPI_IN_PLACE, subspacesData + sentRecvd, static_cast<int>(chunkSize), dtype,
                    MPI_SUM, globalComm);
      sentRecvd += chunkSize;
    }
    MPI_Allreduce(MPI_IN_PLACE, subspacesData + sentRecvd,
                  static_cast<int>(subspacesDataSize - sentRecvd), dtype, MPI_SUM, globalComm);
  } else {
    if (theMPISystem()->getGlobalReduceRank() == globalReduceRankThatCollects) {
      // I am the reduce rank that collects the data
      while ((subspacesDataSize - sentRecvd) / chunkSize > 0) {
        MPI_Reduce(MPI_IN_PLACE, subspacesData + sentRecvd, static_cast<int>(chunkSize), dtype,
                   MPI_SUM, globalReduceRankThatCollects, globalComm);
        sentRecvd += chunkSize;
      }
      MPI_Reduce(MPI_IN_PLACE, subspacesData + sentRecvd,
                 static_cast<int>(subspacesDataSize - sentRecvd), dtype, MPI_SUM,
                 globalReduceRankThatCollects, globalComm);
    } else {
      // I only need to send
      while ((subspacesDataSize - sentRecvd) / chunkSize > 0) {
        MPI_Reduce(subspacesData + sentRecvd, MPI_IN_PLACE, static_cast<int>(chunkSize), dtype,
                   MPI_SUM, globalReduceRankThatCollects, globalComm);
        sentRecvd += chunkSize;
      }
      MPI_Reduce(subspacesData + sentRecvd, MPI_IN_PLACE,
                 static_cast<int>(subspacesDataSize - sentRecvd), dtype, MPI_SUM,
                 globalReduceRankThatCollects, globalComm);
    }
  }
}

/**
 * @brief custom MPI_Op for adding indexed elements, to avoid intermediate copies
 *
 * @tparam FG_ELEMENT the type of the elements in the sparse grid
 * @param invec input vector
 * @param inoutvec input/output vector
 * @param len unused, must be 1
 * @param dtype MPI datatype
 * @throw std::runtime_error if len>1, datatype not understood, num_addresses != 0, or
 * num_integers%2
 */
template <typename FG_ELEMENT>
void addIndexedElements(void* invec, void* inoutvec, int* len, MPI_Datatype* dtype) {
  // cf. https://stackoverflow.com/a/29286769
  if (*len != 1) {
    throw std::runtime_error("addIndexedElements: len>1 not implemented.");
  }
  int num_integers, num_addresses, num_datatypes, combiner;
  MPI_Type_get_envelope(*dtype, &num_integers, &num_addresses, &num_datatypes, &combiner);
  if (combiner != MPI_COMBINER_INDEXED || num_datatypes != 1) {
    throw std::runtime_error("addIndexedElements: do not understand datatype.");
  }
  if (num_addresses != 0 || num_integers % 2 != 1) {
    throw std::runtime_error("addIndexedElements: num_addresses != 0 or num_integers%2 != 1.");
  }
  int arrayOfInts[num_integers];
  MPI_Aint addresses[1]; // actually, the requried size is num_addresses=0 but forbidden by ISO C++
  MPI_Datatype types[1];
  MPI_Type_get_contents(*dtype, num_integers, num_addresses, num_datatypes, arrayOfInts, addresses,
                        types);
  if (types[0] != getMPIDatatype(abstraction::getabstractionDataType<FG_ELEMENT>())) {
    throw std::runtime_error("addIndexedElements: datatype not as expected.");
  }

  int numBlocks = (num_integers - 1) / 2;
  int* arrayOfBlocklengths = arrayOfInts + 1;
  int* arrayOfDisplacements = arrayOfBlocklengths + numBlocks;
#pragma omp parallel for default(none) firstprivate( \
        numBlocks, arrayOfDisplacements, arrayOfBlocklengths, invec, inoutvec) schedule(guided)
  for (int i = 0; i < numBlocks; ++i) {
    FG_ELEMENT* inoutElements = reinterpret_cast<FG_ELEMENT*>(inoutvec) + arrayOfDisplacements[i];
    FG_ELEMENT* inElements = reinterpret_cast<FG_ELEMENT*>(invec) + arrayOfDisplacements[i];
#pragma omp simd linear(inoutElements, inElements : 1)
    for (int j = 0; j < arrayOfBlocklengths[i]; ++j) {
      *inoutElements += *inElements;
      ++inoutElements;
      ++inElements;
    }
  }
}

/**
 * @brief get chunked subspaces for a DistributedSparseGrid
 *
 * partitions the subspaces of \p dsg into chunks of not more than \p maxChunkSize size in total,
 * but at least one subspace
 *
 * @tparam FG_ELEMENT the type of the elements in the sparse grid
 * @tparam SubspaceIndexContainer a container of subspace indices
 * @param dsg the DistributedSparseGrid
 * @param siContainer the container of subspace indices
 * @param maxChunkSize the maximum size of a chunk, in number of elements
 * @return a vector of sets of subspace indices, each set representing a chunk
 */
template <typename FG_ELEMENT, typename SubspaceIndexContainer>
std::vector<std::set<typename AnyDistributedSparseGrid::SubspaceIndexType>>& getChunkedSubspaces(
    const DistributedSparseGridUniform<FG_ELEMENT>& dsg, const SubspaceIndexContainer& siContainer,
    int maxChunkSize) {
  static thread_local std::vector<std::set<typename AnyDistributedSparseGrid::SubspaceIndexType>>
      subspaceIndicesChunks;
  subspaceIndicesChunks.clear();

  auto subspaceIt = siContainer.cbegin();
  while (subspaceIt != siContainer.cend()) {
    static thread_local std::set<typename AnyDistributedSparseGrid::SubspaceIndexType>
        chunkSubspaces;
    chunkSubspaces.clear();
    SubspaceSizeType chunkDataSize = 0;
    // select chunk of not more than chunkSize, but at least one index
    auto nextAddedDataSize = dsg.getDataSize(*subspaceIt);
    do {
      chunkDataSize += nextAddedDataSize;
      chunkSubspaces.insert(*subspaceIt);
      ++subspaceIt;
    } while (subspaceIt != siContainer.cend() &&
             (nextAddedDataSize = dsg.getDataSize(*subspaceIt)) &&
             (chunkDataSize + nextAddedDataSize) < static_cast<size_t>(maxChunkSize));
    subspaceIndicesChunks.push_back(std::move(chunkSubspaces));
  }
  return subspaceIndicesChunks;
}

/**
 * @brief get the reduction datatypes for a DistributedSparseGridUniform
 *
 * the datatypes will cover all sparse grid subspaces in \p subspaces
 *
 * @tparam FG_ELEMENT the type of the elements in the sparse grid
 * @tparam SubspaceIndexContainer a container of subspace indices
 * @param dsg the DistributedSparseGrid
 * @param subspaces the container of subspace indices
 * @param maxMiBToSendPerThread the maximum number of MiB to send per OpenMP thread (~= per core) at
 * once
 * @return a vector of pairs of the start index of the subspace and the MPI datatype for the
 * reduction
 */
template <typename FG_ELEMENT, typename SubspaceIndexContainer>
std::vector<std::pair<typename AnyDistributedSparseGrid::SubspaceIndexType, MPI_Datatype>>
getReductionDatatypes(const DistributedSparseGridUniform<FG_ELEMENT>& dsg,
                      const SubspaceIndexContainer& subspaces, uint32_t maxMiBToSendPerThread) {
  // like for sparse grid reduce, allow only up to 16MiB per reduction
  auto chunkSize = getGlobalReduceChunkSize<FG_ELEMENT>(maxMiBToSendPerThread);
  std::vector<std::pair<typename AnyDistributedSparseGrid::SubspaceIndexType, MPI_Datatype>>
      datatypesByStartIndex;

  // iterate subspaces and create MPI datatypes
  // get chunked subspaces for this data type
  {
    auto& chunkedSubspaces = getChunkedSubspaces(dsg, subspaces, chunkSize);
    for (auto& subspacesChunk : chunkedSubspaces) {
      const FG_ELEMENT* rawDataStartFirst = dsg.getData(*subspacesChunk.cbegin());

      // create datatype for this chunk
      static thread_local std::vector<int> arrayOfBlocklengths, arrayOfDisplacements;
      arrayOfBlocklengths.clear();
      arrayOfDisplacements.clear();
      arrayOfBlocklengths.reserve(subspacesChunk.size());
      arrayOfDisplacements.reserve(subspacesChunk.size());
      for (const auto& ss : subspacesChunk) {
        arrayOfBlocklengths.push_back(dsg.getDataSize(ss));
        arrayOfDisplacements.push_back(static_cast<int>(dsg.getData(ss) - rawDataStartFirst));
      }

      MPI_Datatype myIndexedDatatype;
      MPI_Type_indexed(static_cast<int>(subspacesChunk.size()), arrayOfBlocklengths.data(),
                       arrayOfDisplacements.data(),
                       getMPIDatatype(abstraction::getabstractionDataType<FG_ELEMENT>()),
                       &myIndexedDatatype);
      MPI_Type_commit(&myIndexedDatatype);
      datatypesByStartIndex.push_back(std::make_pair(*subspacesChunk.cbegin(), myIndexedDatatype));
    }
  }
  assert(!datatypesByStartIndex.empty() && "No datatypes created");
  return datatypesByStartIndex;
}

/**
 * @brief perform distributed global subspace reduce or outgroup sparse grid reduce
 *
 * sets of subspaces are reduced between the groups that have space allocated for them, with a
 * communicator for each set (this is why subspace reduce cannot be used for really large
 * combination schemes / numbers of subspaces and process groups)
 *
 * The differences between sparse grid reduce, subspace reduce, and outgroup reduce are discussed in
 * more detail in http://elib.uni-stuttgart.de/handle/11682/14229 , section 4.4.
 *
 * @tparam SparseGridType (derived from) DistributedSparseGridUniform
 * @tparam communicateAllAllocated if false, all subspaces stored by communicator in the distributed
 * sparse grid are communicated (different sets -> subspace reduce); if true, all allocated subspaces are
 * communicated (a single set of subspaces -> outgroup reduce).
 * @param dsg sparse grid of type SparseGridType
 * @param maxMiBToSendPerThread the maximum number of MiB to send per OpenMP thread (~= per core) at
 * once
 * @param globalReduceRankThatCollects the rank that collects the data, or MPI_PROC_NULL if
 * allreduce
 */
template <typename SparseGridType, bool communicateAllAllocated = false>
void distributedGlobalSubspaceReduce(SparseGridType& dsg, uint32_t maxMiBToSendPerThread,
                                     RankType globalReduceRankThatCollects = MPI_PROC_NULL) {
  assert(dsg.isSubspaceDataCreated() && "Only perform reduce with allocated data");

  MPI_Op indexedAdd;
  MPI_Op_create(addIndexedElements<typename SparseGridType::ElementType>, true, &indexedAdd);

#pragma omp parallel if (dsg.getSubspacesByCommunicator().size() > 1) default(none) \
    shared(dsg, indexedAdd) firstprivate(globalReduceRankThatCollects, maxMiBToSendPerThread)
#pragma omp for schedule(dynamic)
  for (size_t commIndex = 0; commIndex < dsg.getSubspacesByCommunicator().size(); ++commIndex) {
    const std::pair<CommunicatorType,
                    std::vector<typename AnyDistributedSparseGrid::SubspaceIndexType>>&
        commAndItsSubspaces = dsg.getSubspacesByCommunicator()[commIndex];

    // get reduction datatypes
    std::vector<std::pair<typename AnyDistributedSparseGrid::SubspaceIndexType, MPI_Datatype>>
        datatypesByStartIndex;
    if constexpr (communicateAllAllocated) {
      datatypesByStartIndex =
          getReductionDatatypes(dsg, dsg.getCurrentlyAllocatedSubspaces(), maxMiBToSendPerThread);
      assert(datatypesByStartIndex.size() == 1);
    } else {
      datatypesByStartIndex =
          getReductionDatatypes(dsg, commAndItsSubspaces.second, maxMiBToSendPerThread);
    }
    /*
    // // this would be best for outgroup reduce, but leads to MPI truncation
    // // errors if not ordered (desynchronization between MPI ranks on the same communicators
    // // probably)
    // #pragma omp parallel if (dsg.getSubspacesByCommunicator().size() == 1) default(none) \
    // shared(dsg, indexedAdd, datatypesByStartIndex, commAndItsSubspaces)
    // #pragma omp for ordered schedule(static)
    */
    for (size_t datatypeIndex = 0; datatypeIndex < datatypesByStartIndex.size(); ++datatypeIndex) {
      // reduce for each datatype
      auto& subspaceStartIndex = datatypesByStartIndex[datatypeIndex].first;
      auto& comm = commAndItsSubspaces.first;
      auto& datatype = datatypesByStartIndex[datatypeIndex].second;
      if (globalReduceRankThatCollects == MPI_PROC_NULL) {
        // #pragma omp ordered
        [[maybe_unused]] auto success = MPI_Allreduce(MPI_IN_PLACE, dsg.getData(subspaceStartIndex),
                                                      1, datatype, indexedAdd, comm);
        assert(success == MPI_SUCCESS);

      } else {  // reduce towards only one rank
        assert(dsg.getSubspacesByCommunicator().size() == 1);
        if (theMPISystem()->getGlobalReduceRank() == globalReduceRankThatCollects) {
          // I am the reduce rank that collects the data
          MPI_Reduce(MPI_IN_PLACE, dsg.getData(subspaceStartIndex), 1, datatype, indexedAdd,
                     globalReduceRankThatCollects, comm);
        } else {
          // I only need to send
          MPI_Reduce(dsg.getData(subspaceStartIndex), MPI_IN_PLACE, 1, datatype, indexedAdd,
                     globalReduceRankThatCollects, comm);
        }
      }
      // free datatype -- MPI standard says it will be kept until operation is finished
      MPI_Type_free(&(datatype));
    }
  }

  MPI_Op_free(&indexedAdd);
}

/**
 * @brief Sends all subspace data sizes of \p dsg to rank \p dest in communicator \p comm
 */
template <typename SparseGridType>
static void sendSubspaceDataSizes(SparseGridType& dsg, RankType dest, CommunicatorType comm) {
  assert(dsg.getNumSubspaces() > 0);

  const std::vector<int>& subspacesDataSizes = dsg.getSubspaceDataSizes();
  MPI_Send(subspacesDataSizes.data(), subspacesDataSizes.size(), MPI_INT, dest,
           TRANSFER_SUBSPACE_DATA_SIZES_TAG, comm);
}

/** 
 * @brief Performs a max-allreduce in communicator comm with subspace sizes of the sparse grids
 *
 * Typically used in the GlobalReduceComm to ensure that all workers have the same subspace sizes for sparse grid reduce.
 * 
 * @tparam SparseGridType (derived from) DistributedSparseGridUniform
 * @param dsg the sparse grid of type SparseGridType
 * @param comm the communicator
 */
template <typename SparseGridType>
void reduceSubspaceSizes(SparseGridType& dsg, CommunicatorType comm) {
  assert(dsg.getNumSubspaces() > 0);

  // prepare for MPI call in globalReduceComm
  MPI_Datatype dtype = getMPIDatatype(abstraction::getabstractionDataType<SubspaceSizeType>());

  // perform allreduce
  assert(dsg.getNumSubspaces() <
         static_cast<AnyDistributedSparseGrid::SubspaceIndexType>(std::numeric_limits<int>::max()));
  MPI_Allreduce(MPI_IN_PLACE, dsg.getSubspaceDataSizes().data(),
                static_cast<int>(dsg.getNumSubspaces()), dtype, MPI_MAX, comm);
  // assume that the sizes changed, the buffer might be the wrong size now
  dsg.deleteSubspaceData();
}

/**
 * @brief Performs a max-reduce in communicator globalReduceComm with subspace sizes of the sparse grids
 * 
 * this together with broadcastSubspaceSizes can be used as a two-step replacement for reduceSubspaceSizes.
 * 
 * @tparam SparseGridType (derived from) DistributedSparseGridUniform
 * @param dsg the sparse grid of type SparseGridType
 * @param globalReduceRankThatCollects the rank that collects the data
 * @param globalReduceComm the communicator
 */
template <typename SparseGridType>
void maxReduceSubspaceSizesAcrossGroups(
    SparseGridType& dsg, RankType globalReduceRankThatCollects,
    CommunicatorType globalReduceComm = theMPISystem()->getGlobalReduceComm()) {
  auto numSubspaces = static_cast<int>(dsg.getNumSubspaces());
  MPI_Datatype dtype = getMPIDatatype(abstraction::getabstractionDataType<SubspaceSizeType>());
  if (theMPISystem()->getGlobalReduceRank() == globalReduceRankThatCollects) {
    MPI_Reduce(MPI_IN_PLACE, dsg.getSubspaceDataSizes().data(), numSubspaces, dtype, MPI_MAX,
               globalReduceRankThatCollects, globalReduceComm);
  } else {
    MPI_Reduce(dsg.getSubspaceDataSizes().data(), MPI_IN_PLACE, numSubspaces, dtype, MPI_MAX,
               globalReduceRankThatCollects, globalReduceComm);
  }
  dsg.deleteSubspaceData();
}

/**
 * @brief Broadcasts the subspace sizes of the sparse grid from rank \p sendingRank in communicator \p comm
 *
 * this together with maxReduceSubspaceSizesAcrossGroups can be used as a two-step replacement for reduceSubspaceSizes.
 * 
 * @tparam SparseGridType (derived from) DistributedSparseGridUniform
 * @param dsg the sparse grid of type SparseGridType
 * @param comm the communicator
 * @param sendingRank the rank that sends the data
 */
template <typename SparseGridType>
void broadcastSubspaceSizes(SparseGridType& dsg, CommunicatorType comm, RankType sendingRank) {
  assert(dsg.getNumSubspaces() > 0);
  MPI_Datatype dtype = getMPIDatatype(abstraction::getabstractionDataType<SubspaceSizeType>());

  // perform broadcast
  assert(dsg.getNumSubspaces() <
         static_cast<AnyDistributedSparseGrid::SubspaceIndexType>(std::numeric_limits<int>::max()));
  MPI_Bcast(dsg.getSubspaceDataSizes().data(), static_cast<int>(dsg.getNumSubspaces()), dtype,
            sendingRank, comm);
  // assume that the sizes changed, the buffer might be the wrong size now
  dsg.deleteSubspaceData();
}

/**
 * @brief gather all subspace sizes on the different ranks of \p comm to the rank \p collectorRank
 * 
 * can be used for the widely-distributed combination technique
 */
template <typename SparseGridType>
void sendSubspaceSizesWithGather(SparseGridType& dsg, CommunicatorType comm,
                                 RankType collectorRank) {
  auto numSubspaces = static_cast<int>(dsg.getNumSubspaces());
  assert(numSubspaces > 0);
  MPI_Datatype dtype = getMPIDatatype(abstraction::getabstractionDataType<SubspaceSizeType>());

  // perform gather (towards third level manager)
  // send size of buffer to collectorRank
  MPI_Gather(&numSubspaces, 1, MPI_INT, nullptr, 0, MPI_INT, collectorRank, comm);

  // send subspace sizes to manager
  MPI_Gatherv(dsg.getSubspaceDataSizes().data(), numSubspaces, dtype, nullptr, nullptr, nullptr,
              dtype, collectorRank, comm);
}

/**
 * @brief receive all subspace sizes on the different ranks of \p comm from the rank \p collectorRank
 * 
 * can be used for the widely-distributed combination technique
 */
template <typename SparseGridType>
void receiveSubspaceSizesWithScatter(SparseGridType& dsg, CommunicatorType comm,
                                     RankType collectorRank) {
  auto numSubspaces = static_cast<int>(dsg.getNumSubspaces());
  assert(numSubspaces > 0);
  assert(static_cast<size_t>(numSubspaces) == dsg.getSubspaceDataSizes().size());
  MPI_Datatype dtype = getMPIDatatype(abstraction::getabstractionDataType<SubspaceSizeType>());

  // receive updated sizes from manager
  MPI_Scatterv(nullptr, 0, nullptr, dtype, dsg.getSubspaceDataSizes().data(), numSubspaces, dtype,
               collectorRank, comm);
  // assume that the sizes changed, the buffer might be the wrong size now
  dsg.deleteSubspaceData();
}

/**
 * @brief send the raw sparse grid data to rank \p dest in communicator \p comm
 */
template <typename SparseGridType>
static void sendDsgData(SparseGridType& dsg, RankType dest, CommunicatorType comm) {
  typename SparseGridType::ElementType* data = dsg.getRawData();
  auto dataSize = dsg.getRawDataSize();
  MPI_Datatype dataType =
      getMPIDatatype(abstraction::getabstractionDataType<typename SparseGridType::ElementType>());

  size_t sentRecvd = 0;
  while ((dataSize - sentRecvd) / INT_MAX > 0) {
    MPI_Send(data + sentRecvd, (int)INT_MAX, dataType, dest, TRANSFER_DSGU_DATA_TAG, comm);
    sentRecvd += INT_MAX;
  }
  MPI_Send(data + sentRecvd, (int)(dataSize - sentRecvd), dataType, dest, TRANSFER_DSGU_DATA_TAG,
           comm);
}

/**
 * @brief receive the raw sparse grid data from rank \p source in communicator \p comm
 */
template <typename SparseGridType>
static void recvDsgData(SparseGridType& dsg, RankType source, CommunicatorType comm) {
  typename SparseGridType::ElementType* data = dsg.getRawData();
  auto dataSize = dsg.getRawDataSize();
  MPI_Datatype dataType =
      getMPIDatatype(abstraction::getabstractionDataType<typename SparseGridType::ElementType>());

  size_t sentRecvd = 0;
  while ((dataSize - sentRecvd) / INT_MAX > 0) {
    MPI_Recv(data + sentRecvd, (int)INT_MAX, dataType, source, TRANSFER_DSGU_DATA_TAG, comm,
             MPI_STATUS_IGNORE);
    sentRecvd += INT_MAX;
  }
  MPI_Recv(data + sentRecvd, (int)(dataSize - sentRecvd), dataType, source, TRANSFER_DSGU_DATA_TAG,
           comm, MPI_STATUS_IGNORE);
}

/**
 * @brief asynchronous Bcast of the raw sparse grid data in the communicator \p comm
 *
 * from rank \p root to all others
 *
 * @tparam SparseGridType the type of the sparse grid
 * @param dsg the sparse grid
 * @param root the rank that broadcasts the data
 * @param comm the communicator
 * @param request pointer to an already-allocated MPI_Request
 */
template <typename SparseGridType>
static void asyncBcastDsgData(SparseGridType& dsg, RankType root, CommunicatorType comm,
                              MPI_Request* request) {
  if (dsg.getRawDataSize() >= INT_MAX) {
    throw std::runtime_error(
        "asyncBcastDsgData: Dsg is too large and can not be "
        "transferred in a single MPI Call (not "
        "supported yet) try a more refined"
        "decomposition");
  }

  typename SparseGridType::ElementType* data = dsg.getRawData();
  int dataSize = static_cast<int>(dsg.getRawDataSize());
  MPI_Datatype dataType =
      getMPIDatatype(abstraction::getabstractionDataType<typename SparseGridType::ElementType>());

  [[maybe_unused]] auto success = MPI_Ibcast(data, dataSize, dataType, root, comm, request);
  assert(success == MPI_SUCCESS);
}

/**
 * @brief asynchronous Bcast of the raw sparse grid data in the communicator \p comm
 * 
 * but only for the subspaces allocated by all process groups (and used by at least two groups -> outgroup sparse grid reduce)
 * 
 * @tparam SparseGridType the type of the sparse grid
 * @param dsg the sparse grid of type SparseGridType
 * @param root the rank that broadcasts the data
 * @param comm the communicator
 * @param request pointer to an already-allocated MPI_Request
 */
template <typename SparseGridType, bool communicateAllAllocated = false>
static void asyncBcastOutgroupDsgData(SparseGridType& dsg, RankType root, CommunicatorType comm,
                                      MPI_Request* request) {
  // assert that the dsg is set up for outgroup sparse grid reduce
  assert(dsg.getSubspacesByCommunicator().size() < 2);
  if (!dsg.getSubspacesByCommunicator().empty()) {
    std::vector<std::pair<typename AnyDistributedSparseGrid::SubspaceIndexType, MPI_Datatype>>
        datatypesByStartIndex;
    if constexpr (communicateAllAllocated) {  // assuming byte limit is met by selection of
                                              // allocated spaces
      datatypesByStartIndex = getReductionDatatypes(dsg, dsg.getCurrentlyAllocatedSubspaces(),
                                                    std::numeric_limits<uint32_t>::max());
    } else {
      const auto& commAndItsSubspaces = dsg.getSubspacesByCommunicator()[0];
      datatypesByStartIndex = getReductionDatatypes(dsg, commAndItsSubspaces.second,
                                                    std::numeric_limits<uint32_t>::max());
    }
    assert(datatypesByStartIndex.size() == 1);

    auto& subspaceStartIndex = datatypesByStartIndex[0].first;
    auto& datatype = datatypesByStartIndex[0].second;
    [[maybe_unused]] auto success =
        MPI_Ibcast(dsg.getData(subspaceStartIndex), 1, datatype, root, comm, request);
    assert(success == MPI_SUCCESS);
  }
}

}  // namespace CombiCom

} /* namespace combigrid */

#endif /* COMBICOM_HPP_ */
