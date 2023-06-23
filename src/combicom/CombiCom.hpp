#ifndef COMBICOM_HPP_
#define COMBICOM_HPP_

#include "fullgrid/DistributedFullGrid.hpp"
#include "fullgrid/FullGrid.hpp"
#include "mpi/MPISystem.hpp"
#include "mpi/MPITags.hpp"
#include "sparsegrid/DistributedSparseGridUniform.hpp"

namespace combigrid {

namespace CombiCom {

template <typename FG_ELEMENT>
void FGReduce(FullGrid<FG_ELEMENT>& fg, RankType r, MPI_Comm comm) {
  if (!fg.isGridCreated()) fg.createFullGrid();

  auto& sendbuf = fg.getElementVector();
  std::vector<FG_ELEMENT> recvbuf(sendbuf.size());

  int myrank;
  MPI_Comm_rank(comm, &myrank);

  MPI_Datatype dtype =
      abstraction::getMPIDatatype(abstraction::getabstractionDataType<FG_ELEMENT>());

  if (myrank == r) {
    MPI_Reduce(MPI_IN_PLACE, &sendbuf[0], static_cast<int>(sendbuf.size()), dtype, MPI_SUM, r,
               comm);
  } else {
    MPI_Reduce(&sendbuf[0], &recvbuf[0], static_cast<int>(sendbuf.size()), dtype, MPI_SUM, r, comm);
  }
}

template <typename FG_ELEMENT>
void FGAllreduce(FullGrid<FG_ELEMENT>& fg, MPI_Comm comm) {
  if (!fg.isGridCreated()) fg.createFullGrid();

  auto& sendbuf = fg.getElementVector();
  std::vector<FG_ELEMENT> recvbuf(sendbuf.size());

  MPI_Datatype dtype =
      abstraction::getMPIDatatype(abstraction::getabstractionDataType<FG_ELEMENT>());

  int myrank;
  MPI_Comm_rank(comm, &myrank);

  MPI_Allreduce(MPI_IN_PLACE, &sendbuf[0], static_cast<int>(sendbuf.size()), dtype, MPI_SUM, comm);
}

template <typename SparseGridType>
int sumAndCheckSubspaceSizes(const SparseGridType& dsg,
                             const std::vector<SubspaceSizeType>& subspaceSizes) {
  int bsize = 0;

  for (size_t i = 0; i < subspaceSizes.size(); ++i) {
    // check for implementation errors, the reduced subspace size should not be
    // different from the size of already initialized subspaces
    bool check = (subspaceSizes[i] == 0 || dsg.getDataSize(i) == 0 ||
                  subspaceSizes[i] == int(dsg.getDataSize(i)));

    if (!check) {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      std::cout << "l = " << dsg.getLevelVector(i) << " "
                << "rank = " << rank << " "
                << "ssize = " << subspaceSizes[i] << " "
                << "dsize = " << dsg.getDataSize(i) << std::endl;
      assert(false);
    }

    bsize += subspaceSizes[i];
  }
  return bsize;
}

template <typename SparseGridType>
bool sumAndCheckSubspaceSizes(const SparseGridType& dsg) {
  auto bsize = sumAndCheckSubspaceSizes(dsg, dsg.getSubspaceDataSizes());
  return bsize == dsg.getRawDataSize();
}

template <typename SG_ELEMENT>
static int getGlobalReduceChunkSize(size_t maxBytesToSend = 16777216) {
  // auto chunkSize = std::numeric_limits<int>::max();
  // allreduce up to 16MiB at a time (when using double precision)
  // constexpr size_t sixteenMiBinBytes = 16777216;
  auto numOMPThreads = theMPISystem()->getNumOpenMPThreads();
  int chunkSize = static_cast<int>(maxBytesToSend / sizeof(SG_ELEMENT) / numOMPThreads);
  assert(chunkSize == 2097152 || (numOMPThreads > 1) || (!std::is_same_v<SG_ELEMENT, double>));
  return chunkSize;
}

/*** In this method the global reduction of the distributed sparse grid is
 * performed. The global sparse grid is decomposed geometrically according to
 * the domain decomposition (see Variant 2 in chapter 3.3.3 in marios diss). We
 * first collect all subspace parts of our domain part and fill all missing
 * subspaces with 0. We then perform an MPI_Allreduce with all processes from
 * other process groups that own the same geometrical area. Here we follow the
 * Sparse Grid Reduce strategy from chapter 2.7.2 in marios diss.
 */
template <typename SparseGridType>
void distributedGlobalSparseGridReduce(
    SparseGridType& dsg, RankType globalReduceRankThatCollects = MPI_PROC_NULL,
    MPI_Comm globalComm = theMPISystem()->getGlobalReduceComm()) {
  assert(globalComm != MPI_COMM_NULL);
  assert(dsg.isSubspaceDataCreated() && "Only perform reduce with allocated data");

  typename SparseGridType::ElementType* subspacesData = dsg.getRawData();
  size_t subspacesDataSize = dsg.getRawDataSize();
  assert(subspacesDataSize == dsg.getAccumulatedDataSize());

  // global reduce
  MPI_Datatype dtype = abstraction::getMPIDatatype(
      abstraction::getabstractionDataType<typename SparseGridType::ElementType>());

  auto chunkSize = getGlobalReduceChunkSize<typename SparseGridType::ElementType>();
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

// cf. https://stackoverflow.com/a/29286769
template <typename FG_ELEMENT>
void addIndexedElements(void* invec, void* inoutvec, int* len, MPI_Datatype* dtype) {
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
  MPI_Aint addresses[num_addresses];
  MPI_Datatype types[num_datatypes];
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
    const auto subspaceItBefore = subspaceIt;
    SubspaceSizeType chunkDataSize = 0;
    // select chunk of not more than chunkSize, but at least one index
    auto nextAddedDataSize = dsg.getDataSize(*subspaceIt);
    do {
      chunkDataSize += nextAddedDataSize;
      chunkSubspaces.insert(*subspaceIt);
      ++subspaceIt;
    } while (subspaceIt != siContainer.cend() &&
             (nextAddedDataSize = dsg.getDataSize(*subspaceIt)) &&
             (chunkDataSize + nextAddedDataSize) < maxChunkSize);
    subspaceIndicesChunks.push_back(std::move(chunkSubspaces));
  }
  return subspaceIndicesChunks;
}

template <typename FG_ELEMENT, typename SubspaceIndexContainer>
std::vector<std::pair<typename AnyDistributedSparseGrid::SubspaceIndexType, MPI_Datatype>>
getReductionDatatypes(const DistributedSparseGridUniform<FG_ELEMENT>& dsg,
                      const SubspaceIndexContainer& subspaces, size_t maxBytesToSend = 16777216) {
  // like for sparse grid reduce, allow only up to 16MiB per reduction
  auto chunkSize = getGlobalReduceChunkSize<FG_ELEMENT>(maxBytesToSend);
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
        arrayOfDisplacements.push_back(dsg.getData(ss) - rawDataStartFirst);
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

template <typename SparseGridType, bool communicateAllAllocated = false>
void distributedGlobalSubspaceReduce(SparseGridType& dsg,
                                     RankType globalReduceRankThatCollects = MPI_PROC_NULL) {
  assert(dsg.isSubspaceDataCreated() && "Only perform reduce with allocated data");

  MPI_Op indexedAdd;
  MPI_Op_create(addIndexedElements<typename SparseGridType::ElementType>, true, &indexedAdd);

#pragma omp parallel if (dsg.getSubspacesByCommunicator().size() > 1) default(none) \
    shared(dsg, indexedAdd) firstprivate(globalReduceRankThatCollects)
#pragma omp for schedule(dynamic)
  for (size_t commIndex = 0; commIndex < dsg.getSubspacesByCommunicator().size(); ++commIndex) {
    const std::pair<CommunicatorType,
                    std::vector<typename AnyDistributedSparseGrid::SubspaceIndexType>>&
        commAndItsSubspaces = dsg.getSubspacesByCommunicator()[commIndex];

    // get reduction datatypes
    std::vector<std::pair<typename AnyDistributedSparseGrid::SubspaceIndexType, MPI_Datatype>>
        datatypesByStartIndex;
    if constexpr (communicateAllAllocated) {
      datatypesByStartIndex = getReductionDatatypes(dsg, dsg.getCurrentlyAllocatedSubspaces());
    } else {
      datatypesByStartIndex = getReductionDatatypes(dsg, commAndItsSubspaces.second);
    }

    // // this would be best for outgroup reduce, but leads to MPI truncation
    // // errors if not ordered (desynchronization between MPI ranks on the same communicators
    // // probably)
    // #pragma omp parallel if (dsg.getSubspacesByCommunicator().size() == 1) default(none) \
    // shared(dsg, indexedAdd, datatypesByStartIndex, commAndItsSubspaces)
    // #pragma omp for ordered schedule(static)
    for (size_t datatypeIndex = 0; datatypeIndex < datatypesByStartIndex.size(); ++datatypeIndex) {
      // reduce for each datatype
      auto& subspaceStartIndex = datatypesByStartIndex[datatypeIndex].first;
      auto& comm = commAndItsSubspaces.first;
      auto& datatype = datatypesByStartIndex[datatypeIndex].second;
      if (globalReduceRankThatCollects == MPI_PROC_NULL) {
        // #pragma omp ordered
        auto success = MPI_Allreduce(MPI_IN_PLACE, dsg.getData(subspaceStartIndex), 1, datatype,
                                     indexedAdd, comm);
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
 * Sends the raw dsg data to the destination process in communicator comm.
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
 * Recvs the raw dsg data from the source process in communicator comm.
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
 * Asynchronous Bcast of the raw dsg data in the communicator comm.
 */
template <typename SparseGridType>
static MPI_Request asyncBcastDsgData(SparseGridType& dsg, RankType root, CommunicatorType comm) {
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
  MPI_Request request = MPI_REQUEST_NULL;

  auto success = MPI_Ibcast(data, dataSize, dataType, root, comm, &request);
  assert(success == MPI_SUCCESS);
  return request;
}

/**
 * Asynchronous Bcast of the raw dsg data in the communicator comm.
 */
template <typename SparseGridType>
[[nodiscard]] static MPI_Request asyncBcastOutgroupDsgData(SparseGridType& dsg, RankType root,
                                                           CommunicatorType comm) {
  assert(dsg.getSubspacesByCommunicator().size() < 2);
  MPI_Request request = MPI_REQUEST_NULL;
  if (!dsg.getSubspacesByCommunicator().empty()) {
    const auto& commAndItsSubspaces = dsg.getSubspacesByCommunicator()[0];
    auto datatypesByStartIndex =
        getReductionDatatypes(dsg, commAndItsSubspaces.second, std::numeric_limits<size_t>::max());
    assert(datatypesByStartIndex.size() == 1);

    auto& subspaceStartIndex = datatypesByStartIndex[0].first;
    auto& datatype = datatypesByStartIndex[0].second;
    auto success = MPI_Ibcast(dsg.getData(subspaceStartIndex), 1, datatype, root, comm, &request);
    assert(success == MPI_SUCCESS);
  }
  return request;
}

/**
 * Sends all subspace data sizes to the receiver in communicator comm.
 */
template <typename SparseGridType>
static void sendSubspaceDataSizes(SparseGridType& dsg, RankType dest, CommunicatorType comm) {
  assert(dsg.getNumSubspaces() > 0);

  const std::vector<int>& subspacesDataSizes = dsg.getSubspaceDataSizes();
  MPI_Send(subspacesDataSizes.data(), subspacesDataSizes.size(), MPI_INT, dest,
           TRANSFER_SUBSPACE_DATA_SIZES_TAG, comm);
}

/** Performs a max allreduce in comm with subspace sizes of each dsg
 *
 * After calling, all workers which share the same spatial decomposition will
 * have the same subspace sizes and therefor. in the end have equally sized dsgs.
 */
template <typename SparseGridType>
void reduceSubspaceSizes(SparseGridType& dsg, CommunicatorType comm) {
  assert(dsg.getNumSubspaces() > 0);

  // prepare for MPI call in globalReduceComm
  MPI_Datatype dtype = getMPIDatatype(abstraction::getabstractionDataType<SubspaceSizeType>());

  // perform allreduce
  assert(dsg.getNumSubspaces() < static_cast<SubspaceSizeType>(std::numeric_limits<int>::max()));
  MPI_Allreduce(MPI_IN_PLACE, dsg.getSubspaceDataSizes().data(),
                static_cast<int>(dsg.getNumSubspaces()), dtype, MPI_MAX, comm);
  // assume that the sizes changed, the buffer might be the wrong size now
  dsg.deleteSubspaceData();
}

template <typename SparseGridType>
void broadcastSubspaceSizes(SparseGridType& dsg, CommunicatorType comm, RankType sendingRank) {
  assert(dsg.getNumSubspaces() > 0);
  MPI_Datatype dtype = getMPIDatatype(abstraction::getabstractionDataType<SubspaceSizeType>());

  // perform broadcast
  assert(dsg.getNumSubspaces() < static_cast<SubspaceSizeType>(std::numeric_limits<int>::max()));
  MPI_Bcast(dsg.getSubspaceDataSizes().data(), static_cast<int>(dsg.getNumSubspaces()), dtype,
            sendingRank, comm);
  // assume that the sizes changed, the buffer might be the wrong size now
  dsg.deleteSubspaceData();
}

template <typename SparseGridType>
void localMaxReduceSubspaceSizes(SparseGridType& dsg, CommunicatorType comm, RankType sendingRank) {
  assert(dsg.getNumSubspaces() > 0);
  MPI_Datatype dtype = getMPIDatatype(abstraction::getabstractionDataType<SubspaceSizeType>());

  // make a copy of the subspace sizes
  std::vector<SubspaceSizeType> subspaceSizes = dsg.getSubspaceDataSizes();

  // perform broadcast
  assert(dsg.getNumSubspaces() < static_cast<SubspaceSizeType>(std::numeric_limits<int>::max()));
  MPI_Bcast(subspaceSizes.data(), static_cast<int>(dsg.getNumSubspaces()), dtype, sendingRank,
            comm);
  assert(sumAndCheckSubspaceSizes(dsg, subspaceSizes));

  // perform max reduction
  std::transform(dsg.getSubspaceDataSizes().cbegin(), dsg.getSubspaceDataSizes().cend(),
                 subspaceSizes.cbegin(), dsg.getSubspaceDataSizes().begin(),
                 [](SubspaceSizeType a, SubspaceSizeType b) { return std::max(a, b); });

  // assume that the sizes changed, the buffer might be the wrong size now
  dsg.deleteSubspaceData();
}

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
}

template <typename SparseGridType>
void sendSubspaceSizesWithGather(SparseGridType& dsg, CommunicatorType comm,
                                 RankType collectorRank) {
  auto numSubspaces = static_cast<int>(dsg.getNumSubspaces());
  assert(numSubspaces > 0);
  MPI_Datatype dtype = getMPIDatatype(abstraction::getabstractionDataType<SubspaceSizeType>());

  // perform gather (towards tl-manager)
  // send size of buffer to manager
  MPI_Gather(&numSubspaces, 1, MPI_INT, nullptr, 0, MPI_INT, collectorRank, comm);

  // send subspace sizes to manager
  MPI_Gatherv(dsg.getSubspaceDataSizes().data(), numSubspaces, dtype, nullptr, nullptr, nullptr,
              dtype, collectorRank, comm);
}

template <typename SparseGridType>
void receiveSubspaceSizesWithScatter(SparseGridType& dsg, CommunicatorType comm,
                                     RankType collectorRank) {
  auto numSubspaces = static_cast<int>(dsg.getNumSubspaces());
  assert(numSubspaces > 0);
  assert(numSubspaces == dsg.getSubspaceDataSizes().size());
  MPI_Datatype dtype = getMPIDatatype(abstraction::getabstractionDataType<SubspaceSizeType>());

  // receive updated sizes from manager
  MPI_Scatterv(nullptr, 0, nullptr, dtype, dsg.getSubspaceDataSizes().data(), numSubspaces, dtype,
               collectorRank, comm);
  // assume that the sizes changed, the buffer might be the wrong size now
  dsg.deleteSubspaceData();
}

}  // namespace CombiCom

} /* namespace combigrid */

#endif /* COMBICOM_HPP_ */
