#ifndef COMBICOM_HPP_
#define COMBICOM_HPP_

#include "fullgrid/DistributedFullGrid.hpp"
#include "fullgrid/FullGrid.hpp"
#include "mpi/MPISystem.hpp"
#include "utils/Stats.hpp"

namespace combigrid {

class CombiCom {
 private:
  template <typename FG_ELEMENT>
  static int sumAndCheckSubspaceSizes(const DistributedSparseGridUniform<FG_ELEMENT>& dsg,
                                      const std::vector<SubspaceSizeType>& subspaceSizes);

 public:
  // reduced fg will ONLY be available on r
  template <typename FG_ELEMENT>
  static void FGReduce(FullGrid<FG_ELEMENT>& fg, RankType r, MPI_Comm comm);

  // reduced fg will be available on all member of comm
  template <typename FG_ELEMENT>
  static void FGAllreduce(FullGrid<FG_ELEMENT>& fg, MPI_Comm comm);

  template <typename FG_ELEMENT>
  static void distributedGlobalReduce(DistributedSparseGridUniform<FG_ELEMENT>& dsg,
                                      RankType globalReduceRankThatCollects = MPI_PROC_NULL);

  template <typename FG_ELEMENT>
  static bool sumAndCheckSubspaceSizes(const DistributedSparseGridUniform<FG_ELEMENT>& dsg);
};

template <typename FG_ELEMENT>
void CombiCom::FGReduce(FullGrid<FG_ELEMENT>& fg, RankType r, MPI_Comm comm) {
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
    MPI_Reduce(&sendbuf[0], &recvbuf[0], static_cast<int>(sendbuf.size()), dtype, MPI_SUM, r,
               comm);
  }
}

template <typename FG_ELEMENT>
void CombiCom::FGAllreduce(FullGrid<FG_ELEMENT>& fg, MPI_Comm comm) {
  if (!fg.isGridCreated()) fg.createFullGrid();

  auto& sendbuf = fg.getElementVector();
  std::vector<FG_ELEMENT> recvbuf(sendbuf.size());

  MPI_Datatype dtype =
      abstraction::getMPIDatatype(abstraction::getabstractionDataType<FG_ELEMENT>());

  int myrank;
  MPI_Comm_rank(comm, &myrank);

  MPI_Allreduce(MPI_IN_PLACE, &sendbuf[0], static_cast<int>(sendbuf.size()), dtype, MPI_SUM, comm);
}

template <typename FG_ELEMENT>
int CombiCom::sumAndCheckSubspaceSizes(const DistributedSparseGridUniform<FG_ELEMENT>& dsg,
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

template <typename FG_ELEMENT>
bool CombiCom::sumAndCheckSubspaceSizes(const DistributedSparseGridUniform<FG_ELEMENT>& dsg) {
  auto bsize = CombiCom::sumAndCheckSubspaceSizes(dsg, dsg.getSubspaceDataSizes());
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
template <typename FG_ELEMENT>
void CombiCom::distributedGlobalReduce(DistributedSparseGridUniform<FG_ELEMENT>& dsg,
                                       RankType globalReduceRankThatCollects) {
  // get global communicator for this operation
  MPI_Comm mycomm = theMPISystem()->getGlobalReduceComm();

  assert(mycomm != MPI_COMM_NULL);

  assert(dsg.isSubspaceDataCreated() && "Only perform reduce with allocated data");

  FG_ELEMENT* subspacesData = dsg.getRawData();
  size_t subspacesDataSize = dsg.getRawDataSize();

  // global reduce
  MPI_Datatype dtype =
      abstraction::getMPIDatatype(abstraction::getabstractionDataType<FG_ELEMENT>());

  // auto chunkSize = std::numeric_limits<int>::max();
  // allreduce up to 16MiB at a time (when using double precision)
  auto chunkSize = 2097152;
  size_t sentRecvd = 0;
  if (globalReduceRankThatCollects == MPI_PROC_NULL) {
    while ((subspacesDataSize - sentRecvd) / chunkSize > 0) {
      MPI_Allreduce(MPI_IN_PLACE, subspacesData + sentRecvd, static_cast<int>(chunkSize), dtype,
                    MPI_SUM, mycomm);
      sentRecvd += chunkSize;
    }
    MPI_Allreduce(MPI_IN_PLACE, subspacesData + sentRecvd,
                  static_cast<int>(subspacesDataSize - sentRecvd), dtype, MPI_SUM, mycomm);
  } else {
    if (theMPISystem()->getGlobalReduceRank() == globalReduceRankThatCollects) {
      // I am the reduce rank that collects the data
      while ((subspacesDataSize - sentRecvd) / chunkSize > 0) {
        MPI_Reduce(MPI_IN_PLACE, subspacesData + sentRecvd, static_cast<int>(chunkSize), dtype,
                   MPI_SUM, globalReduceRankThatCollects, mycomm);
        sentRecvd += chunkSize;
      }
      MPI_Reduce(MPI_IN_PLACE, subspacesData + sentRecvd,
                 static_cast<int>(subspacesDataSize - sentRecvd), dtype, MPI_SUM,
                 globalReduceRankThatCollects, mycomm);
    } else {
      // I only need to send
      while ((subspacesDataSize - sentRecvd) / chunkSize > 0) {
        MPI_Reduce(subspacesData + sentRecvd, MPI_IN_PLACE, static_cast<int>(chunkSize), dtype,
                   MPI_SUM, globalReduceRankThatCollects, mycomm);
        sentRecvd += chunkSize;
      }
      MPI_Reduce(subspacesData + sentRecvd, MPI_IN_PLACE,
                 static_cast<int>(subspacesDataSize - sentRecvd), dtype, MPI_SUM,
                 globalReduceRankThatCollects, mycomm);
    }
  }
}

} /* namespace combigrid */

#endif /* COMBICOM_HPP_ */
