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

/*** In this method the global reduction of the distributed sparse grid is
 * performed. The global sparse grid is decomposed geometrically according to
 * the domain decomposition (see Variant 2 in chapter 3.3.3 in marios diss). We
 * first collect all subspace parts of our domain part and fill all missing
 * subspaces with 0. We then perform an MPI_Allreduce with all processes from
 * other process groups that own the same geometrical area. Here we follow the
 * Sparse Grid Reduce strategy from chapter 2.7.2 in marios diss.
 */
template <typename SparseGridType>
void distributedGlobalSparseGridReduce(SparseGridType& dsg,
                                       RankType globalReduceRankThatCollects = MPI_PROC_NULL,
                                       MPI_Comm globalComm = theMPISystem()->getGlobalReduceComm()) {
  assert(globalComm != MPI_COMM_NULL);
  assert(dsg.isSubspaceDataCreated() && "Only perform reduce with allocated data");

  typename SparseGridType::ElementType* subspacesData = dsg.getRawData();
  size_t subspacesDataSize = dsg.getRawDataSize();

  // global reduce
  MPI_Datatype dtype = abstraction::getMPIDatatype(
      abstraction::getabstractionDataType<typename SparseGridType::ElementType>());

  // auto chunkSize = std::numeric_limits<int>::max();
  // allreduce up to 16MiB at a time (when using double precision)
  auto chunkSize = 2097152;
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
