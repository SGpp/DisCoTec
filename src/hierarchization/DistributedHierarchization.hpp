#ifndef DISTRIBUTEDHIERARCHIZATION_HPP_
#define DISTRIBUTEDHIERARCHIZATION_HPP_

#include "boost/lexical_cast.hpp"
#include "fullgrid/DistributedFullGrid.hpp"
#include "utils/IndexVector.hpp"
#include "utils/PowerOfTwo.hpp"
#include "utils/Stats.hpp"

namespace combigrid {

/* The RemoteDataSlice is meant to store a (d-1)-dimensional block of a
 * d-dimensional DistributedFullGrid. The RemoteDataSlice is d-dimensional,
 * but has exactly one point in (at least) one dimension.
 */
template <typename FG_ELEMENT>
class RemoteDataSlice {
 public:
  /** Constructor
   *
   * \param[in] sizes d-dimensional vector specifying the extent in each dim
   * \param[in] dim1d reduced dimension in which the block has only one point
   * \param[in] keyIndex the index of the (d-1)-dimensional subgrid in the
   * d-dimensional grid
   * \param[in] lowerBounds lower bounds of the subdomain where the remote data
   *            comes from. this is required for address calculations
   */
  RemoteDataSlice(IndexType size, IndexType keyIndex) : index1d_(keyIndex) {
    data_.resize(size);
  }

  inline const FG_ELEMENT* getData(size_t idx) const { return &data_[idx]; }

  inline std::vector<FG_ELEMENT>& getElementVector() { return data_; }

  // return index of (d-1)-dimensional subgrid in the d-dimensional grid
  inline IndexType getKeyIndex() const { return index1d_; }

 private:
  // global index of (d-1)-dimensional subgrid in the d-dimensional grid
  IndexType index1d_;

  // data vector
  std::vector<FG_ELEMENT> data_;
};

template <typename FG_ELEMENT>
class RemoteDataCollector : public std::vector<RemoteDataSlice<FG_ELEMENT>> {
 public:
  RemoteDataSlice<FG_ELEMENT>* addDataSlice(IndexType size, IndexType keyIndex) {
    auto& newest = this->emplace_back(size, keyIndex);
    return &newest;
  }
};

template <typename FG_ELEMENT>
static void checkLeftSuccesors(IndexType checkIdx, const IndexType& rootIdx, const DimType& dim,
                               const DistributedFullGrid<FG_ELEMENT>& dfg,
                               std::map<RankType, std::set<IndexType>>& OneDIndices,
                               const LevelType& lmin);

template <typename FG_ELEMENT>
static void checkRightSuccesors(IndexType checkIdx, const IndexType& rootIdx, const DimType& dim,
                                const DistributedFullGrid<FG_ELEMENT>& dfg,
                                std::map<RankType, std::set<IndexType>>& OneDIndices,
                                const LevelType& lmin);

template <typename FG_ELEMENT>
static IndexType checkPredecessors(IndexType idx, const DimType& dim,
                                   const DistributedFullGrid<FG_ELEMENT>& dfg,
                                   std::map<RankType, std::set<IndexType>>& OneDIndices,
                                   const LevelType& lmin);

template <typename FG_ELEMENT>
static IndexType getNextIndex1d(const DistributedFullGrid<FG_ELEMENT>& dfg, DimType d,
                                IndexType idx1d);

template <typename FG_ELEMENT>
static IndexType getFirstIndexOfLevel1d(const DistributedFullGrid<FG_ELEMENT>& dfg, DimType d,
                                        LevelType l);

/**
 * @brief helper function for data exchange, have only one MPI_Isend/Irecv per neighboring rank
 *
 * @param send1dIndices : a vector which holds (for each other rank) a set of own 1D indices to send
 * to
 * @param recv1dIndices : a vector which holds (for each other rank) a set of 1D indices to receive
 * from
 * @param dfg : the DistributedFullGrid where the own values are stored
 * @param dim : the dimension in which we want to exchange
 * @param remoteData : the data structure into which the received data will be stored
 */
template <typename FG_ELEMENT>
static void sendAndReceiveIndicesBlock(const std::map<RankType, std::set<IndexType>>& send1dIndices,
                                       const std::map<RankType, std::set<IndexType>>& recv1dIndices,
                                       const DistributedFullGrid<FG_ELEMENT>& dfg, DimType dim,
                                       RemoteDataCollector<FG_ELEMENT>& remoteData) {
  assert(remoteData.empty());
  // count elements of input indices
  auto numSend = send1dIndices.size();
  auto numRecv = recv1dIndices.size();
  // buffer for requests
  std::vector<MPI_Request> sendRequests(numSend);
  std::vector<MPI_Request> recvRequests(numRecv);

  // create general subarray pointing to the first d-1 dimensional slice
  MPI_Datatype mysubarray;
  {
    // sizes of local grid
    std::vector<int> sizes(dfg.getLocalSizes().begin(), dfg.getLocalSizes().end());
    // sizes of subarray ( full size except dimension d )
    std::vector<int> subsizes = sizes;
    subsizes[dim] = 1;
    // start
    std::vector<int> starts(dfg.getDimension(), 0);
    // create subarray view on data
    MPI_Type_create_subarray(static_cast<int>(dfg.getDimension()), &sizes[0], &subsizes[0],
                             &starts[0], MPI_ORDER_FORTRAN, dfg.getMPIDatatype(), &mysubarray);
    MPI_Type_commit(&mysubarray);
  }
  MPI_Aint dfgStartAddr;
  MPI_Get_address(dfg.getData(), &dfgStartAddr);

  numSend = 0;
  numRecv = 0;
  /*
  // #pragma omp parallel shared(sendRequests, numSend, send1dIndices, recvRequests, numRecv, \
//                                 remoteData, recv1dIndices, dfg, mysubarray)              \
//     firstprivate(dfgStartAddr, dim) default(none)
  // no benefit from parallelization here
  */
  {
#pragma omp for schedule(static) nowait
    for (size_t x = 0; x < send1dIndices.size(); ++x) {
      auto mapIt = send1dIndices.cbegin();
      std::advance(mapIt, x);
      const auto& r = mapIt->first;
      const auto& indices = mapIt->second;
      assert(!indices.empty());
      if (!indices.empty()) {
        // make datatype hblock for all indices
        MPI_Datatype myHBlock;
        std::vector<MPI_Aint> displacements;
        displacements.reserve(indices.size());
        for (const auto& index : indices) {
          // convert global 1d index i to local 1d index
          IndexType localLinearIndex =
              (index - dfg.getLowerBounds()[dim]) * dfg.getLocalOffsets()[dim];
#ifndef NDEBUG
          {
            static thread_local IndexVector lidxvec(dfg.getDimension(), 0);
            lidxvec.resize(dfg.getDimension());
            static thread_local IndexVector gidxvec;
            gidxvec = dfg.getLowerBounds();
            gidxvec[dim] = index;
            [[maybe_unused]] bool tmp = dfg.getLocalVectorIndex(gidxvec, lidxvec);
            assert(tmp && "index to be send not in local domain");
            assert(localLinearIndex == dfg.getLocalLinearIndex(lidxvec));
          }
#endif
          MPI_Aint addr;
          MPI_Get_address(&(dfg.getData()[localLinearIndex]), &addr);
          auto d = MPI_Aint_diff(addr, dfgStartAddr);
          displacements.push_back(d);
        }
        // cannot use MPI_Type_create_indexed_block as subarrays may overlap
        MPI_Type_create_hindexed_block(static_cast<int>(indices.size()), 1, displacements.data(),
                                       mysubarray, &myHBlock);
        MPI_Type_commit(&myHBlock);

        // send to rank r, use first global index as tag
        {
          int dest = static_cast<int>(r);
          int tag = static_cast<int>(*(indices.begin()));
          size_t sendIndex;
#pragma omp atomic capture
          sendIndex = numSend++;
          MPI_Isend(dfg.getData(), 1, myHBlock, dest, tag, dfg.getCommunicator(),
                    &sendRequests[sendIndex]);
        }
        MPI_Type_free(&myHBlock);
      }
    }

#pragma omp for schedule(dynamic) nowait
    for (size_t x = 0; x < recv1dIndices.size(); ++x) {
      auto mapIt = recv1dIndices.cbegin();
      std::advance(mapIt, x);
      const auto& r = mapIt->first;
      const auto& indices = mapIt->second;
      assert(!indices.empty());
      const IndexVector& lowerBoundsNeighbor = dfg.getLowerBounds(static_cast<int>(r));

      std::vector<FG_ELEMENT*> bufs;
      bufs.reserve(indices.size());
      IndexType size = dfg.getNrLocalElements() / dfg.getLocalSizes()[dim];
      int bsize = static_cast<int>(size);

      for (const auto& index : indices) {
        // create RemoteDataSlice to store the subarray
        RemoteDataSlice<FG_ELEMENT>* backData = nullptr;
#pragma omp critical
        backData = remoteData.addDataSlice(size, index);
        auto& buf = backData->getElementVector();
        bufs.push_back(buf.data());
        assert(bsize == static_cast<int>(buf.size()));
      }
      {
        // make datatype hblock for all indices
        MPI_Datatype myHBlock;
        MPI_Aint firstBufAddr;
        MPI_Get_address(bufs[0], &firstBufAddr);
        std::vector<MPI_Aint> displacements;
        displacements.resize(bufs.size());
        for (size_t bufIndex = 0; bufIndex < bufs.size(); ++bufIndex) {
          MPI_Aint addr;
          MPI_Get_address(bufs[bufIndex], &addr);
          displacements[bufIndex] = MPI_Aint_diff(addr, firstBufAddr);
        }
        assert(displacements[0] == 0);
        MPI_Type_create_hindexed_block(static_cast<int>(indices.size()), bsize,
                                       displacements.data(), dfg.getMPIDatatype(), &myHBlock);
        MPI_Type_commit(&myHBlock);
        // start recv operation, use first global index as tag
        {
          int src = static_cast<int>(r);
          int tag = static_cast<int>(*(indices.begin()));
          size_t recvIndex;
#pragma omp atomic capture
          recvIndex = numRecv++;
          MPI_Irecv(static_cast<void*>(bufs[0]), 1, myHBlock, src, tag, dfg.getCommunicator(),
                    &recvRequests[recvIndex]);
        }
        MPI_Type_free(&myHBlock);
      }
    }
  }

  // free after parallel region
  MPI_Type_free(&mysubarray);

  assert(sendRequests.size() == numSend);
  assert(recvRequests.size() == numRecv);
  MPI_Waitall(static_cast<int>(sendRequests.size()), sendRequests.data(), MPI_STATUSES_IGNORE);
  MPI_Waitall(static_cast<int>(recvRequests.size()), recvRequests.data(), MPI_STATUSES_IGNORE);
}

/**
 * @brief share all data with neighboring processes in one dimension
 *
 * @param dfg : the DistributedFullGrid where the own values are stored
 * @param dim : the dimension in which we want to exchange
 * @param remoteData : the data structure into which the received data will be stored
 */
template <typename FG_ELEMENT>
static void exchangeAllData1d(const DistributedFullGrid<FG_ELEMENT>& dfg, DimType dim,
                              RemoteDataCollector<FG_ELEMENT>& remoteData) {
  // send every index to all neighboring ranks in dimension dim
  const auto globalIdxMax = dfg.globalNumPointsInDimension(dim);
  const IndexType idxMin = dfg.getFirstGlobal1dIndex(dim);
  const IndexType idxMax = dfg.getLastGlobal1dIndex(dim);

  const auto poleNeighbors = dfg.getCartesianUtils().getAllMyPoleNeighborRanks(dim);

  std::set<IndexType> allMyIndices;
  // #pragma omp parallel for shared(idxMin, idxMax, allMyIndices) default(none)
  for (IndexType i = idxMin; i <= idxMax; ++i) {
#pragma omp critical
    allMyIndices.insert(i);
  }
  std::map<RankType, std::set<IndexType>> send1dIndices;
  std::map<RankType, std::set<IndexType>> recv1dIndices;
  /*
  // #pragma omp parallel for schedule(static) default(none) \
  //   shared(poleNeighbors, send1dIndices, allMyIndices)
  */
  for (const auto& r : poleNeighbors) {
#pragma omp critical
    send1dIndices[r] = allMyIndices;
  }
  /*
  // #pragma omp parallel for schedule(static) default(none) \
  //     shared(poleNeighbors, recv1dIndices)
  */
  for (const auto& r : poleNeighbors) {
#pragma omp critical
    recv1dIndices[r] = {};
  }

#pragma omp parallel shared(dfg, recv1dIndices) firstprivate(dim, globalIdxMax, idxMin, idxMax)
  {
    // all other points that are not ours can be received from their owners
#pragma omp for schedule(static) nowait
    for (IndexType i = 0; i < idxMin; ++i) {
      // get rank which has i and add to recv list
      // TODO would be easier to iterate the whole range of each neighbor
      const int r = dfg.getNeighbor1dFromAxisIndex(dim, i);
      if (r >= 0) {
#pragma omp critical
        recv1dIndices.at(r).insert(i);
      }
    }
#pragma omp for schedule(dynamic) nowait
    for (IndexType i = idxMax + 1; i < globalIdxMax; ++i) {
      // get rank which has i and add to recv list
      const int r = dfg.getNeighbor1dFromAxisIndex(dim, i);
      if (r >= 0) {
#pragma omp critical
        recv1dIndices.at(r).insert(i);
      }
    }
  }

  return sendAndReceiveIndicesBlock(send1dIndices, recv1dIndices, dfg, dim, remoteData);
}

/**
 * @brief share data with neighboring processes in one dimension, but such that
 *        every process gets only the data for direct hierarchical predecssors of its own points
 *
 * @param dfg : the DistributedFullGrid where the own values are stored
 * @param dim : the dimension in which we want to exchange
 * @param remoteData : the data structure into which the received data will be stored
 */
template <typename FG_ELEMENT>
static void exchangeData1d(const DistributedFullGrid<FG_ELEMENT>& dfg, DimType dim,
                           RemoteDataCollector<FG_ELEMENT>& remoteData, LevelType lmin = 0) {
  // create buffers for every rank
  std::map<RankType, std::set<IndexType>> recv1dIndices;
  std::map<RankType, std::set<IndexType>> send1dIndices;

  // main loop
  const IndexType idxMin = dfg.getFirstGlobal1dIndex(dim);
  const IndexType idxMax = dfg.getLastGlobal1dIndex(dim);
  const bool oneSidedBoundary = dfg.returnBoundaryFlags()[dim] == 1;
  const auto globalIdxMax = dfg.globalNumPointsInDimension(dim);
  const LevelType lmax = dfg.getLevels()[dim];

  IndexType idx = idxMin;

  // for hierarchization, we only need to exchange the direct predecessors
  while (idx <= idxMax) {
    LevelType lidx = dfg.getLevel(dim, idx);
    // check if direct successors of idx outside local domain
    {
      for (LevelType l = static_cast<LevelType>(lidx + 1); l <= lmax; ++l) {
        if (l > lmin) {
          auto ldiff = static_cast<LevelType>(lmax - l);
          auto idiff = static_cast<IndexType>(powerOfTwo[ldiff]);

          // left successor, right successor
          for (const auto& indexShift : {-1, 1}) {
            IndexType sIdx = idx + indexShift * idiff;
            // if sIdx is outside of my domain, but still in the global domain
            if ((indexShift < 0 && sIdx >= 0 && sIdx < idxMin) ||
                (indexShift > 0 && sIdx > idxMax && sIdx < globalIdxMax)) {
              // get rank which has sIdx and add to send list
              int r = dfg.getNeighbor1dFromAxisIndex(dim, sIdx);
              if (r >= 0) send1dIndices[r].insert(idx);
            }
            // also send to the left if on lower boundary and periodic
            if (oneSidedBoundary && idx == 0 && sIdx < 0) {
              assert(dfg.getCartesianUtils().isOnLowerBoundaryInDimension(dim));
              // wrap around
              auto mirroredShiftedIndex = globalIdxMax + sIdx;
              assert(mirroredShiftedIndex > 0);
              if (mirroredShiftedIndex > idxMax) {
                int r = dfg.getNeighbor1dFromAxisIndex(dim, mirroredShiftedIndex);
                if (r >= 0) send1dIndices[r].insert(idx);
              }
            }
          }
        }
      }
    }
    auto ldiff = static_cast<LevelType>(lmax - lidx);
    auto idiff = static_cast<IndexType>(powerOfTwo[ldiff]);
    // check if predecessors of idx outside local domain
    IndexType pIdx;
    // only recv if this index idx needs to be hierarchized
    if (lidx > lmin) {
      // left, right predecessor
      for (const auto& indexShift : {-1, 1}) {
        pIdx = idx + indexShift * idiff;
        // if we are not on the boundary level, and
        // pIdx is outside of my domain, but still in the global domain
        if (lidx > 0 && ((indexShift < 0 && pIdx >= 0 && pIdx < idxMin) ||
                         (indexShift > 0 && pIdx > idxMax && pIdx < globalIdxMax))) {
          // get rank which has predecessor and add to list of indices to recv
          int r = dfg.getNeighbor1dFromAxisIndex(dim, pIdx);
          if (r >= 0) recv1dIndices[r].insert(pIdx);
        }
        // in case I need the upper boundary index and we have periodic boundary
        if (oneSidedBoundary && indexShift > 0 && pIdx == globalIdxMax) {
          // request index 0
          int r = dfg.getNeighbor1dFromAxisIndex(dim, 0);
          if (r >= 0 && r != dfg.getRank()) {
            recv1dIndices[r].insert(0);
          }
        }
      }
    } else {
      pIdx = idx + idiff;
    }
    if (lidx == 0 || pIdx > idxMax) {
      idx = getNextIndex1d(dfg, dim, idx);
    } else {
      // index of right predecessor
      idx = pIdx;
    }
  }

  return sendAndReceiveIndicesBlock(send1dIndices, recv1dIndices, dfg, dim, remoteData);
}

/**
 * @brief share data with neighboring processes in one dimension, but only such that
 *        every process has the data for all hierarchical predecssors of its own points
 *
 * @param dfg : the DistributedFullGrid where the own values are stored
 * @param dim : the dimension in which we want to exchange
 * @param remoteData : the data structure into which the received data will be stored
 */
template <typename FG_ELEMENT>
static void exchangeData1dDehierarchization(const DistributedFullGrid<FG_ELEMENT>& dfg, DimType dim,
                                            RemoteDataCollector<FG_ELEMENT>& remoteData,
                                            LevelType lmin = 0) {
  // create buffers for every rank
  std::map<RankType, std::set<IndexType>> recv1dIndices;
  std::map<RankType, std::set<IndexType>> send1dIndices;

  // main loop
  IndexType idxMin = dfg.getFirstGlobal1dIndex(dim);
  IndexType idxMax = dfg.getLastGlobal1dIndex(dim);
  // LevelType lmax = dfg.getLevels()[dim];

  IndexType idx = idxMin;

  // for dehierarchization, we need to exchange the full tree of
  // successors and predecessors
  while (idx <= idxMax) {
    checkLeftSuccesors(idx, idx, dim, dfg, send1dIndices, lmin);

    checkRightSuccesors(idx, idx, dim, dfg, send1dIndices, lmin);

    idx = checkPredecessors(idx, dim, dfg, recv1dIndices, lmin);
  }

  return sendAndReceiveIndicesBlock(send1dIndices, recv1dIndices, dfg, dim, remoteData);
}

template <typename FG_ELEMENT>
static void checkLeftSuccesors(IndexType checkIdx, const IndexType& rootIdx, const DimType& dim,
                               const DistributedFullGrid<FG_ELEMENT>& dfg,
                               std::map<RankType, std::set<IndexType>>& sendOneDIndices,
                               const LevelType& lmin) {
  const auto levelOfCheckIndex = dfg.getLevel(dim, checkIdx);
  const auto localMinimalIndex = dfg.getFirstGlobal1dIndex(dim);
  const auto lmax = dfg.getLevels()[dim];

  // check left successors of checkIdx
  for (auto l = static_cast<LevelType>(levelOfCheckIndex + 1); l <= lmax; ++l) {
    const auto levelDifference = static_cast<LevelType>(lmax - l);
    const auto indexDifference = static_cast<IndexType>(powerOfTwo[levelDifference]);

    const IndexType leftSuccessorIndex = checkIdx - indexDifference;
    const auto leftSuccessorLevel = dfg.getLevel(dim, leftSuccessorIndex);

    if (checkIdx == rootIdx) {
      // if we are at recursion depth 0, we need to decide if we need to send this root index
      // through this path at all
      if (leftSuccessorLevel <= lmin) {
        continue;
      }
    }

    if (leftSuccessorIndex >= 0 && leftSuccessorIndex < localMinimalIndex) {
      assert(leftSuccessorLevel <= lmax);
      // only send if my successor needs to be dehierarchized
      if (leftSuccessorLevel > lmin) {
        // get rank which has leftSuccessorIndex and add to send list
        int r = dfg.getNeighbor1dFromAxisIndex(dim, leftSuccessorIndex);
        if (r >= 0) sendOneDIndices[r].insert(rootIdx);
      }
    }
    if (leftSuccessorIndex >= 0)
      checkLeftSuccesors(leftSuccessorIndex, rootIdx, dim, dfg, sendOneDIndices, lmin);
  }
}

template <typename FG_ELEMENT>
static void checkRightSuccesors(IndexType checkIdx, const IndexType& rootIdx, const DimType& dim,
                                const DistributedFullGrid<FG_ELEMENT>& dfg,
                                std::map<RankType, std::set<IndexType>>& sendOneDIndices,
                                const LevelType& lmin) {
  const auto levelOfCheckIndex = dfg.getLevel(dim, checkIdx);
  const auto localMaximalIndex = dfg.getLastGlobal1dIndex(dim);
  const auto lmax = dfg.getLevels()[dim];
  bool oneSidedBoundary = dfg.returnBoundaryFlags()[dim] == 1;
  const auto globalIdxMax = dfg.getGlobalSizes()[dim];

  // check right successors of checkIdx
  for (auto l = static_cast<LevelType>(levelOfCheckIndex + 1); l <= lmax; ++l) {
    const auto levelDifference = static_cast<LevelType>(lmax - l);
    const auto indexDifference = static_cast<IndexType>(powerOfTwo[levelDifference]);

    const IndexType rightSuccessorIndex = checkIdx + indexDifference;
    const auto rightSuccessorLevel = dfg.getLevel(dim, rightSuccessorIndex);
    if (checkIdx == rootIdx) {
      // if we are at recursion depth 0, we need to decide if we need to send this root index
      // through this path at all
      if (rightSuccessorLevel <= lmin) {
        continue;
      }
    }

    // only send if my successor needs to be dehierarchized w.r.t. rootIndex
    if (rightSuccessorLevel > lmin) {
      if (rightSuccessorIndex < globalIdxMax && rightSuccessorIndex > localMaximalIndex) {
        assert(rightSuccessorLevel <= lmax);
        // get rank which has rightSuccessorIndex and add to send list
        int r = dfg.getNeighbor1dFromAxisIndex(dim, rightSuccessorIndex);
        if (r >= 0) sendOneDIndices[r].insert(rootIdx);
      }
      // if we are at the lower boundary and are periodic, we need to mirror this index as well
      if (rootIdx == 0 && oneSidedBoundary) {
        auto mirroredRightSuccessorIndex = globalIdxMax - rightSuccessorIndex;
        if (mirroredRightSuccessorIndex > localMaximalIndex) {
          int r = dfg.getNeighbor1dFromAxisIndex(dim, mirroredRightSuccessorIndex);
          if (r >= 0) sendOneDIndices[r].insert(rootIdx);
        }
      }
    }
    if (rightSuccessorIndex < globalIdxMax) {
      checkRightSuccesors(rightSuccessorIndex, rootIdx, dim, dfg, sendOneDIndices, lmin);
    }
  }
}

template <typename FG_ELEMENT>
static IndexType checkPredecessors(IndexType idx, const DimType& dim,
                                   const DistributedFullGrid<FG_ELEMENT>& dfg,
                                   std::map<RankType, std::set<IndexType>>& recvOneDIndices,
                                   const LevelType& lmin) {
  const auto levelOfIndex = dfg.getLevel(dim, idx);
  const auto localMinimalIndex = dfg.getFirstGlobal1dIndex(dim);
  const auto localMaximalIndex = dfg.getLastGlobal1dIndex(dim);

  // if the current level is already small enough, end recursion here since no further predecessors
  // have to be received
  if (levelOfIndex <= lmin || idx < 0) {
    idx = getNextIndex1d(dfg, dim, idx);
    return idx;
  }

  // check if left predecessor outside local domain
  // if returns negative value there's no left predecessor
  const auto leftPredecessorIndex = dfg.getLeftPredecessor(dim, idx);
  const auto leftPredecessorLevel = dfg.getLevel(dim, leftPredecessorIndex);

  if (leftPredecessorIndex >= 0 && leftPredecessorIndex < localMinimalIndex) {
    // get rank which has left predecessor and add to list of indices
    int r = dfg.getNeighbor1dFromAxisIndex(dim, leftPredecessorIndex);
    recvOneDIndices[r].insert(leftPredecessorIndex);
  }
  // check if further dependencies arise from left predecessor
  if (leftPredecessorIndex >= 0 && leftPredecessorLevel >= lmin) {
    checkPredecessors(leftPredecessorIndex, dim, dfg, recvOneDIndices, lmin);
  }

  // check if right predecessor outside local domain
  // if returns negative value there's no right predecessor
  const auto rightPredecessorIndex = dfg.getRightPredecessor(dim, idx);
  const auto rightPredecessorLevel = dfg.getLevel(dim, rightPredecessorIndex);

  if (rightPredecessorIndex < 0) {
    idx = getNextIndex1d(dfg, dim, idx);
    return idx;
  }

  if (rightPredecessorIndex > localMaximalIndex) {
    // get rank which has right predecessor and add to list of indices to recv
    bool oneSidedBoundary = dfg.returnBoundaryFlags()[dim] == 1;
    if (oneSidedBoundary && rightPredecessorIndex == dfg.getGlobalSizes()[dim]) {
      if (localMinimalIndex > 0) {
        // if we need the highest index, and are periodic, need to require 0 instead
        int r = dfg.getNeighbor1dFromAxisIndex(dim, 0);
        if (r != dfg.getRank()) {
          recvOneDIndices[r].insert(0);
        }
      }
    } else {
      int r = dfg.getNeighbor1dFromAxisIndex(dim, rightPredecessorIndex);
      recvOneDIndices[r].insert(rightPredecessorIndex);
      if (rightPredecessorLevel >= lmin) {
        // check if further dependencies arise from right predecessor
        checkPredecessors(rightPredecessorIndex, dim, dfg, recvOneDIndices, lmin);
      }
    }
    idx = getNextIndex1d(dfg, dim, idx);
  } else {
    idx = rightPredecessorIndex;
  }
  return idx;
}

// returns the next one-dimensional global index which fulfills
// min( lmax, l(idx1d) + 1 )
template <typename FG_ELEMENT>
IndexType getNextIndex1d(const DistributedFullGrid<FG_ELEMENT>& dfg, DimType dim, IndexType idx1d) {
  IndexType idxMax = dfg.getLastGlobal1dIndex(dim);
  LevelType lmax = dfg.getLevels()[dim];

  LevelType lold = dfg.getLevel(dim, idx1d);

  for (IndexType i = idx1d + 1; i <= idxMax; ++i) {
    LevelType l = dfg.getLevel(dim, i);

    if (l == lmax || l == lold + 1) {
      return i;
    }
  }

  // no valid index found -> return right neighbor
  return idx1d + 1;
}

template <typename FG_ELEMENT>
static IndexType getFirstIndexOfLevel1d(const DistributedFullGrid<FG_ELEMENT>& dfg, DimType dim,
                                        LevelType l) {
  IndexType idxMax = dfg.getLastGlobal1dIndex(dim);
  IndexType idxMin = dfg.getFirstGlobal1dIndex(dim);

  for (IndexType i = idxMin; i <= idxMax; ++i) {
    if (dfg.getLevel(dim, i) == l) return i;
  }

  // no index on level l found, return -1
  return -1;
}

template <typename FG_ELEMENT, bool periodic = false>
inline void hierarchize_hat_boundary_kernel(FG_ELEMENT* data, LevelType lmax, LevelType lmin = 0) {
  const int lmaxi = static_cast<int>(lmax);
  int ll = lmaxi - 1;
  int steps = (1 << (lmaxi - 1));
  int offset = 1;  // 1 and not 0 because boundary
  int stepsize = 2;
  int parentOffset = 1;

  for (; ll >= lmin; ll--) {
    FG_ELEMENT parentL = 0.5 * data[offset - parentOffset];
    // #pragma omp simd // slowdown!!
    for (int ctr = 0; ctr < steps; ++ctr) {
      const int centralIndex = offset;
      const FG_ELEMENT parentR = 0.5 * data[centralIndex + parentOffset];
      const FG_ELEMENT val1 = data[centralIndex];
      const FG_ELEMENT val2 = val1 - parentL;
      const FG_ELEMENT val3 = val2 - parentR;
      data[centralIndex] = val3;
      parentL = parentR;
      offset += stepsize;
    }

    steps = steps >> 1;
    offset = (1 << (lmaxi - ll));  // boundary case
    parentOffset = stepsize;
    stepsize = stepsize << 1;
  }

  return;
}

/**
 * @brief mass-conserving hierarchization with full weighting scaling function; the hierarchical
 * basis function W corresponds to the "hat wavelet", that is the stencil [-1/2 1 -1/2]
 *
 * @tparam FG_ELEMENT data type on grid
 * @param data pointer to data begin
 * @param lmax maximum level
 * @param start (unused)
 * @param stride (unused)
 * @param lmin minimum level (if > 0, hierarchization is not performed all the way down)
 */
template <typename FG_ELEMENT, bool periodic = false>
inline void hierarchize_full_weighting_boundary_kernel(FG_ELEMENT* data, LevelType lmax,
                                                       LevelType lmin = 0) {
  const int lmaxi = static_cast<int>(lmax);
  const int idxmax = powerOfTwo[lmaxi];

  for (LevelType ldiff = 0; ldiff < lmax - lmin; ++ldiff) {
    const int step_width = powerOfTwo[ldiff];
    // update f at even indices
    if (periodic) {
      // values at 0 and idxmax will be the same
      data[0] = 0.25 * (data[0] + data[idxmax] + data[step_width] + data[idxmax - step_width]);
      data[idxmax] = data[0];
    } else {
      data[0] = 0.5 * (data[0] + data[step_width]);
      data[idxmax] = 0.5 * (data[idxmax] + data[idxmax - step_width]);
    }
    // todo interleave with lower loop for better cache blocking
    auto leftParent = 0.25 * data[step_width];
    for (int i = 2 * step_width; i < idxmax; i += 2 * step_width) {
      const auto rightParent = 0.25 * data[i + step_width];
      data[i] = leftParent + rightParent + 0.5 * data[i];
      leftParent = rightParent;
    }
    // update alpha / hierarchical surplus at odd indices
    leftParent = -0.5 * data[0];
    for (int i = step_width; i < idxmax; i += 2 * step_width) {
      const auto rightParent = -0.5 * data[i + step_width];
      data[i] = leftParent + (rightParent + data[i]);
      leftParent = rightParent;
    }
  }

  return;
}

/**
 * @brief mass-conserving hierarchization with biorthogonal wavelet and scaling function; the
 * hierarchical basis function W corresponds to the wavelet [-1/8 -1/4 3/4 -1/4 -1/8], cf.
 * Cohen-Daubechies-Feauveau
 *
 * @tparam FG_ELEMENT data type on grid
 * @param data pointer to data begin
 * @param lmax maximum level
 * @param lmin minimum level (if > 0, hierarchization is not performed all the way down)
 */
template <typename FG_ELEMENT, bool periodic = false>
inline void hierarchize_biorthogonal_boundary_kernel(FG_ELEMENT* data, LevelType lmax,
                                                     LevelType lmin = 0) {
  const int lmaxi = static_cast<int>(lmax);
  const int idxmax = powerOfTwo[lmaxi];
  //     for l_hierarchization in range(l_max, l_min, -1):
  //         step_width=2**(l_max-l_hierarchization)
  //         for i in range(step_width,2**(l_max),2*step_width):
  //             y[i]= -0.5*y[i-step_width] + y[i] -0.5*y[i+step_width]
  //         y[0] = y[0] + 0.5*y[1*step_width]
  //         y[-1] = y[-1] + 0.5*y[-1-1*step_width]
  //         for i in range(2*step_width,2**(l_max),2*step_width):
  //             y[i] = 0.25*y[i-step_width] + y[i] + 0.25*y[i+step_width]

  for (LevelType ldiff = 0; ldiff < lmax - lmin; ++ldiff) {
    const int step_width = powerOfTwo[ldiff];
      const auto increment = 2 * step_width;
      if (step_width < idxmax) {
        const auto increment = 2 * step_width;
        // do first update outside loop
        FG_ELEMENT leftParent = data[step_width];
        data[step_width] += -0.5 * (leftParent + data[increment]);
        leftParent = data[increment];
        for (int i = increment + step_width; i < idxmax; i += increment) {
          // update alpha / hierarchical surplus at odd indices
          const auto& rightParent = data[i + step_width];
          data[i] = -0.5 * (leftParent + rightParent) + data[i];
          // update f at even indices, here at left parent position
          const auto& leftLeftParent = data[i - increment];
          data[i - step_width] = 0.25 * (leftLeftParent + data[i]) + leftParent;
          leftParent = rightParent;
        }
      }
    if (periodic) {
      // values at 0 and idxmax will be the same
      data[0] =
          0.5 * (data[0] + data[idxmax]) + 0.25 * (data[step_width] + data[idxmax - step_width]);
      data[idxmax] = data[0];
    } else {
      // mass will build up at the boundary; corresponds to 0-neumann-condition
      data[0] = data[0] + 0.5 * data[step_width];
      data[idxmax] = data[idxmax] + 0.5 * data[idxmax - step_width];
    }
  }
}

template <typename FG_ELEMENT, bool periodic = false>
inline void dehierarchize_hat_boundary_kernel(FG_ELEMENT* data, LevelType lmax,
                                              LevelType lmin = 0) {
  const int lmaxi = static_cast<int>(lmax);
  const int lmini = static_cast<int>(lmin);
  int steps = 1 << (lmini);
  int offset = (1 << (lmaxi - lmini - 1));
  int stepsize = (1 << (lmaxi - lmini));
  int parentOffset = offset;
  for (LevelType ll = lmin + 1; ll <= lmax; ++ll) {
    FG_ELEMENT parentL = 0.5 * data[offset - parentOffset];

    // #pragma omp simd //slowdown!
    for (int ctr = 0; ctr < steps; ++ctr) {
      const int centralIndex = offset;
      const FG_ELEMENT parentR = 0.5 * data[centralIndex + parentOffset];
      const FG_ELEMENT val1 = data[centralIndex];
      const FG_ELEMENT val2 = val1 + parentL;
      const FG_ELEMENT val3 = val2 + parentR;
      data[centralIndex] = val3;
      parentL = parentR;
      offset += stepsize;
    }
    steps = steps << 1;
    offset = (1 << (lmaxi - (ll + 1)));  // boundary case
    parentOffset = parentOffset >> 1;
    stepsize = stepsize >> 1;
  }
  return;
}

/**
 * @brief inverse operation to hierarchize_full_weighting_boundary_kernel
 */
template <typename FG_ELEMENT, bool periodic = false>
inline void dehierarchize_full_weighting_boundary_kernel(FG_ELEMENT* data, LevelType lmax,
                                                         LevelType lmin) {
  const int lmaxi = static_cast<int>(lmax);
  const int idxmax = powerOfTwo[lmaxi];

  for (auto ldiff = static_cast<LevelType>(lmax - lmin - 1); ldiff >= 0; --ldiff) {
    int step_width = powerOfTwo[ldiff];
    // update alpha / hierarchical surplus at odd indices
    auto leftParent = 0.5 * data[0];
    for (int i = step_width; i < idxmax; i += 2 * step_width) {
      const auto rightParent = 0.5 * data[i + step_width];
      data[i] = leftParent + (rightParent + data[i]);
      leftParent = rightParent;
    }
    // update f at even indices
    leftParent = -0.5 * data[step_width];
    for (int i = 2 * step_width; i < idxmax; i += 2 * step_width) {
      const auto rightParent = -0.5 * data[i + step_width];
      data[i] = leftParent + rightParent + 2. * data[i];
      leftParent = rightParent;
    }
    if (periodic) {
      // values at 0 and idxmax will be the same
      data[0] = (data[0] + data[idxmax]) - 0.5 * (data[step_width] + data[idxmax - step_width]);
      data[idxmax] = data[0];
    } else {
      data[0] = 2. * data[0] - data[step_width];
      data[idxmax] = 2. * data[idxmax] - data[idxmax - step_width];
    }
  }

  return;
}

/**
 * @brief inverse operation to hierarchize_biorthogonal_boundary_kernel
 */
template <typename FG_ELEMENT, bool periodic = false>
inline void dehierarchize_biorthogonal_boundary_kernel(FG_ELEMENT* data, LevelType lmax,
                                                       LevelType lmin) {
  const int lmaxi = static_cast<int>(lmax);
  const int idxmax = powerOfTwo[lmaxi];

  // for l_hierarchization in range(l_min, l_max, 1):
  //     step_width=2**(l_max-l_hierarchization-1)
  //     y[0] = y[0] - 0.5*y[0+step_width]
  //     y[-1] = y[-1] - 0.5*y[-1-1*step_width]
  //     for i in range(2*step_width,2**(l_max),2*step_width):
  //         y[i] = -0.25*y[i-step_width] + 1*y[i] - 0.25*y[i+step_width]
  //     for i in range(step_width,2**(l_max),2*step_width):
  //         y[i]= 0.5*y[i-step_width] + y[i] + 0.5*y[i+step_width]

  for (auto ldiff = static_cast<LevelType>(lmax - lmin - 1); ldiff >= 0; --ldiff) {
    const int step_width = powerOfTwo[ldiff];
    // update f at even indices
    if (periodic) {
      data[0] =
          -0.25 * (data[step_width] + data[idxmax - step_width]) + 0.5 * (data[0] + data[idxmax]);
      data[idxmax] = data[0];
    } else {
      data[0] = data[0] - 0.5 * data[step_width];
      data[idxmax] = data[idxmax] - 0.5 * data[idxmax - step_width];
    }

    if (step_width < idxmax) {
      const auto increment = 2 * step_width;
      auto leftParent = data[step_width];
      for (int i = 2 * step_width; i < idxmax; i += increment) {
        // update f at even indices
        const auto& rightParent = data[i + step_width];
        data[i] = -0.25 * (leftParent + rightParent) + data[i];
        // update alpha / hierarchical surplus at odd indices, here: at left parent
        const auto& leftLeftParent = data[i - increment];
        data[i - step_width] = 0.5 * (leftLeftParent + data[i]) + leftParent;
        leftParent = rightParent;
      }
      // update last missing alpha
      data[idxmax - step_width] += 0.5 * (data[idxmax - increment] + data[idxmax]);
    }
  }
}

/**
 * @brief  hierarchize a DFG with boundary points in dimension dim
 */
template <typename FG_ELEMENT, void (*HIERARCHIZATION_FCTN)(FG_ELEMENT[], LevelType, LevelType) =
                                   hierarchize_hat_boundary_kernel>
void hierarchizeWithBoundary(DistributedFullGrid<FG_ELEMENT>& dfg,
                             const RemoteDataCollector<FG_ELEMENT>& remoteData, DimType dim,
                             LevelType lmin_n = 0) {
  assert(dfg.returnBoundaryFlags()[dim] > 0);
  const auto& lmax = dfg.getLevels()[dim];
  const auto& numberOfPolesLowerDimensions = dfg.getLocalOffsets()[dim];
  const IndexType jump = numberOfPolesLowerDimensions * dfg.getLocalSizes()[dim];
  const IndexType numberOfPolesHigherDimensions = dfg.getNrLocalElements() / jump;

  // if we are using periodicity, add an entry to tmp for the virtual last value
  const bool oneSidedBoundary = dfg.returnBoundaryFlags()[dim] == 1;
  const IndexType poleLength = dfg.getGlobalSizes()[dim] + (oneSidedBoundary ? 1 : 0);
  FG_ELEMENT* ldata = dfg.getData();
  const IndexType& gstart = dfg.getLowerBounds()[dim];
  const IndexType localOffsetForThisDimension = dfg.getLocalOffsets()[dim];

  // loop over poles
  static thread_local std::vector<FG_ELEMENT> tmp;
#pragma omp parallel for collapse(2) schedule(static)                                       \
    firstprivate(poleLength, ldata, gstart, localOffsetForThisDimension, lmax, lmin_n, dim, \
                     oneSidedBoundary, jump, numberOfPolesLowerDimensions,                  \
                     numberOfPolesHigherDimensions) shared(dfg, remoteData) default(none)
  for (IndexType nHigher = 0; nHigher < numberOfPolesHigherDimensions; ++nHigher) {
    for (IndexType nLower = 0; nLower < numberOfPolesLowerDimensions; ++nLower) {
      const IndexType poleStart = nHigher * jump + nLower;  // local linear index
      const IndexType poleNumber = nLower + nHigher * numberOfPolesLowerDimensions;
#ifndef NDEBUG
      IndexVector localIndexVector(dfg.getDimension());
      // compute global vector index of poleStart, make sure 0 in this dimension
      dfg.getLocalVectorIndex(poleStart, localIndexVector);
      assert(localIndexVector[dim] == 0);
#endif  // NDEBUG
      tmp.resize(poleLength, std::numeric_limits<FG_ELEMENT>::quiet_NaN());
      // go through remote containers, copy remote data
      for (size_t i = 0; i < remoteData.size(); ++i) {
        tmp[remoteData[i].getKeyIndex()] = *remoteData[i].getData(poleNumber);
      }
      // copy local data
      size_t localDataIndex = static_cast<size_t>(poleStart);
      for (IndexType i = 0; i < dfg.getLocalSizes()[dim]; ++i) {
        tmp[gstart + i] = ldata[localDataIndex];
        localDataIndex += localOffsetForThisDimension;
      }

      if (oneSidedBoundary) {
        // assume periodicity
        //  assert(HIERARCHIZATION_FCTN::periodic); //TODO
        tmp[dfg.getGlobalSizes()[dim]] = tmp[0];
        if (!remoteData.empty() && remoteData[0].getKeyIndex() == 0) {
          assert(!std::isnan(std::real(tmp[0])));
        }
      }
      // hierarchize tmp array with hupp function
      HIERARCHIZATION_FCTN(&tmp[0], lmax, lmin_n);

      if (oneSidedBoundary && !remoteData.empty() && remoteData[0].getKeyIndex() == 0) {
        assert(!std::isnan(std::real(tmp[0])));
        assert(tmp[dfg.getGlobalSizes()[dim]] == tmp[0]);
      }

      // copy pole back
      for (IndexType i = dfg.getLocalSizes()[dim] - 1; i > -1; --i) {
        localDataIndex -= localOffsetForThisDimension;
        ldata[localDataIndex] = tmp[gstart + i];
        assert(!std::isnan(std::real(tmp[gstart + i])));
      }
    }
  }
}

/**
 * @brief  hierarchize a DFG without boundary points in dimension dim
 */
template <typename FG_ELEMENT>
void hierarchizeNoBoundary(DistributedFullGrid<FG_ELEMENT>& dfg,
                           std::vector<RemoteDataSlice<FG_ELEMENT>>& remoteData, DimType dim) {
  assert(dfg.returnBoundaryFlags()[dim] == 0);

  LevelType lmax = dfg.getLevels()[dim];
  IndexType size = dfg.getNrLocalElements();
  IndexType stride = dfg.getLocalOffsets()[dim];
  IndexType ndim = dfg.getLocalSizes()[dim];
  IndexType jump = stride * ndim;
  IndexType nbrOfPoles = size / ndim;

  // loop over poles
  std::vector<FG_ELEMENT> tmp(dfg.getGlobalSizes()[dim], std::numeric_limits<double>::quiet_NaN());
  auto ldata = dfg.getData();
  lldiv_t divresult;
  IndexType start;
  IndexType gstart = dfg.getLowerBounds()[dim];

  for (IndexType nn = 0; nn < nbrOfPoles;
       ++nn) {  // integer operations form bottleneck here -- nested loops are twice as slow
    divresult = std::lldiv(nn, stride);
    start = static_cast<decltype(start)>(divresult.quot * jump +
                                         divresult.rem);  // localer lin index start of pole

#ifndef NDEBUG
    IndexVector localIndexVector(dfg.getDimension());
    // compute global vector index of start
    dfg.getLocalVectorIndex(start, localIndexVector);
    assert(localIndexVector[dim] == 0);
#endif  // NDEBUG

    // go through remote containers
    for (size_t i = 0; i < remoteData.size(); ++i) {
      tmp[remoteData[i].getKeyIndex()] = *remoteData[i].getData(nn);
    }

    // copy local data
    for (IndexType i = 0; i < ndim; ++i) tmp[gstart + i] = ldata[start + stride * i];

    // hierarchization kernel
    IndexType idxMax = dfg.getLastGlobal1dIndex(dim);

    for (LevelType l = lmax; l > 0; --l) {
      // get first local point of level and corresponding stride
      IndexType firstOfLevel = getFirstIndexOfLevel1d(dfg, dim, l);
      IndexType parentOffset = static_cast<IndexType>(powerOfTwo[lmax - l]);
      IndexType levelStride = parentOffset * 2;

      // loop over points of this level with level specific stride
      // as long as inside domain
      if (firstOfLevel > -1) {
        for (IndexType idx = firstOfLevel; idx <= idxMax; idx += levelStride) {
          // when no boundary in this dimension we have to check if
          // 1d indices outside domain
          FG_ELEMENT left(0.0);
          FG_ELEMENT right(0.0);

          if (idx - parentOffset > 0) {
            left = tmp[idx - parentOffset];
          }

          if (idx + parentOffset < dfg.getGlobalSizes()[dim]) {
            right = tmp[idx + parentOffset];
          }

          // do calculation
          FG_ELEMENT buf = -0.5 * left;
          tmp[idx] -= 0.5 * right;
          tmp[idx] += buf;
        }
      }
    }

    // copy pole back
    for (IndexType i = 0; i < ndim; ++i) ldata[start + stride * i] = tmp[gstart + i];
  }
}

/**
 * @brief inverse operation for hierarchizeNoBoundary
 */
template <typename FG_ELEMENT>
void dehierarchizeNoBoundary(DistributedFullGrid<FG_ELEMENT>& dfg,
                             RemoteDataCollector<FG_ELEMENT>& remoteData, DimType dim) {
  assert(dfg.returnBoundaryFlags()[dim] == 0);

  auto lmax = dfg.getLevels()[dim];
  auto size = dfg.getNrLocalElements();
  auto stride = dfg.getLocalOffsets()[dim];
  auto ndim = dfg.getLocalSizes()[dim];
  IndexType jump = stride * ndim;
  IndexType nbrOfPoles = size / ndim;

  // loop over poles
  static thread_local std::vector<FG_ELEMENT> tmp;
  tmp.resize(dfg.getGlobalSizes()[dim], std::numeric_limits<double>::quiet_NaN());
  auto ldata = dfg.getData();
  lldiv_t divresult;
  IndexType start;
  IndexType gstart = dfg.getLowerBounds()[dim];

  for (IndexType nn = 0; nn < nbrOfPoles;
       ++nn) {  // integer operations form bottleneck here -- nested loops are twice as slow
    divresult = std::lldiv(nn, stride);
    start = static_cast<decltype(start)>(divresult.quot * jump +
                                         divresult.rem);  // localer lin index start of pole

#ifndef NDEBUG
    IndexVector localIndexVector(dfg.getDimension());
    // compute global vector index of start
    dfg.getLocalVectorIndex(start, localIndexVector);
    assert(localIndexVector[dim] == 0);
#endif  // NDEBUG

    // go through remote containers
    for (size_t i = 0; i < remoteData.size(); ++i) {
      tmp[remoteData[i].getKeyIndex()] = *remoteData[i].getData(nn);
    }

    // copy local data
    for (IndexType i = 0; i < ndim; ++i) tmp[gstart + i] = ldata[start + stride * i];

    // dehierarchization kernel
    for (LevelType l = 2; l <= lmax; ++l) {
      // get first local point of level and corresponding stride
      IndexType parentOffset = static_cast<IndexType>(powerOfTwo[lmax - l]);
      IndexType first = parentOffset - 1;
      IndexType levelStride = parentOffset * 2;

      // loop over points of this level with level specific stride
      // as long as inside domain
      for (IndexType idx = first; idx < dfg.getGlobalSizes()[dim]; idx += levelStride) {
        // when no boundary in this dimension we have to check if
        // 1d indices outside domain
        FG_ELEMENT left(0.0);
        FG_ELEMENT right(0.0);

        if (idx - parentOffset > 0) {
          left = tmp[idx - parentOffset];
        }

        if (idx + parentOffset < dfg.getGlobalSizes()[dim]) {
          right = tmp[idx + parentOffset];
        }

        // do calculation
        FG_ELEMENT buf = 0.5 * left;
        tmp[idx] += 0.5 * right;
        tmp[idx] += buf;
      }
    }

    // copy pole back
    for (IndexType i = 0; i < ndim; ++i) ldata[start + stride * i] = tmp[gstart + i];
  }
}

class DistributedHierarchization {
 public:
  // inplace hierarchization
  template <typename FG_ELEMENT>
  static void hierarchize(DistributedFullGrid<FG_ELEMENT>& dfg, const std::vector<bool>& dims,
                          const std::vector<BasisFunctionBasis*>& hierarchicalBases,
                          const LevelVector& lmin) {
    assert(dfg.getDimension() > 0);
    assert(dfg.getDimension() == dims.size());
    assert(!lmin.empty());
    // exchange data
    for (DimType dim = 0; dim < dfg.getDimension(); ++dim) {
      if (!dims[dim]) continue;
      if (lmin[dim] >= dfg.getLevels()[dim]) continue;

      RemoteDataCollector<FG_ELEMENT> remoteData;
      if (dynamic_cast<HierarchicalHatBasisFunction*>(hierarchicalBases[dim]) != nullptr ||
          dynamic_cast<HierarchicalHatPeriodicBasisFunction*>(hierarchicalBases[dim]) != nullptr) {
        exchangeData1d(dfg, dim, remoteData, lmin[dim]);
      } else {
        exchangeAllData1d(dfg, dim, remoteData);
      }

      if (dfg.returnBoundaryFlags()[dim] > 0) {
        // sorry for the code duplication, could not figure out a clean way
        if (dynamic_cast<HierarchicalHatBasisFunction*>(hierarchicalBases[dim]) != nullptr) {
          hierarchizeWithBoundary<FG_ELEMENT, hierarchize_hat_boundary_kernel<FG_ELEMENT>>(
              dfg, remoteData, dim, lmin[dim]);
        } else if (dynamic_cast<HierarchicalHatPeriodicBasisFunction*>(hierarchicalBases[dim]) !=
                   nullptr) {
          hierarchizeWithBoundary<FG_ELEMENT, hierarchize_hat_boundary_kernel<FG_ELEMENT, true>>(
              dfg, remoteData, dim, lmin[dim]);
        } else if (dynamic_cast<FullWeightingBasisFunction*>(hierarchicalBases[dim]) != nullptr) {
          hierarchizeWithBoundary<FG_ELEMENT,
                                  hierarchize_full_weighting_boundary_kernel<FG_ELEMENT, false>>(
              dfg, remoteData, dim, lmin[dim]);
        } else if (dynamic_cast<FullWeightingPeriodicBasisFunction*>(hierarchicalBases[dim]) !=
                   nullptr) {
          hierarchizeWithBoundary<FG_ELEMENT,
                                  hierarchize_full_weighting_boundary_kernel<FG_ELEMENT, true>>(
              dfg, remoteData, dim, lmin[dim]);
        } else if (dynamic_cast<BiorthogonalBasisFunction*>(hierarchicalBases[dim]) != nullptr) {
          hierarchizeWithBoundary<FG_ELEMENT,
                                  hierarchize_full_weighting_boundary_kernel<FG_ELEMENT, false>>(
              dfg, remoteData, dim, lmin[dim]);
        } else if (dynamic_cast<BiorthogonalPeriodicBasisFunction*>(hierarchicalBases[dim]) !=
                   nullptr) {
          hierarchizeWithBoundary<FG_ELEMENT,
                                  hierarchize_biorthogonal_boundary_kernel<FG_ELEMENT, true>>(
              dfg, remoteData, dim, lmin[dim]);
        } else {
          throw std::logic_error("Not implemented");
        }
      } else {
        if (dynamic_cast<HierarchicalHatBasisFunction*>(hierarchicalBases[dim]) == nullptr) {
          throw std::logic_error("currently only hats supported for non-boundary grids");
        }
        assert(lmin[dim] == 0);
        hierarchizeNoBoundary(dfg, remoteData, dim);
      }
      remoteData.clear();
    }
  }

  template <typename FG_ELEMENT, class T = HierarchicalHatBasisFunction>
  static void hierarchizeHierachicalBasis(DistributedFullGrid<FG_ELEMENT>& dfg,
                                          const std::vector<bool>& dims,
                                          LevelVector lmin = LevelVector(0)) {
    T basisFctn;
    std::vector<BasisFunctionBasis*> bases(dfg.getDimension(), &basisFctn);
    if (lmin.size() == 0) {
      lmin = LevelVector(dims.size(), 0);
    }
    return hierarchize<FG_ELEMENT>(dfg, dims, bases, lmin);
  }

  template <typename FG_ELEMENT>
  static void hierarchize(DistributedFullGrid<FG_ELEMENT>& dfg, LevelVector lmin = LevelVector(0)) {
    std::vector<bool> dims(dfg.getDimension(), true);
    if (lmin.size() == 0) {
      lmin = LevelVector(dfg.getDimension(), 0);
    }
    return hierarchizeHierachicalBasis<FG_ELEMENT, HierarchicalHatBasisFunction>(dfg, dims, lmin);
  }

  // inplace dehierarchization
  template <typename FG_ELEMENT>
  static void dehierarchize(DistributedFullGrid<FG_ELEMENT>& dfg, const std::vector<bool>& dims,
                            const std::vector<BasisFunctionBasis*>& hierarchicalBases,
                            const LevelVector& lmin) {
    assert(!lmin.empty());
    assert(dfg.getDimension() > 0);
    assert(dfg.getDimension() == dims.size());
    for (DimType dim = 0; dim < dfg.getDimension(); ++dim) {
      if (!dims[dim]) continue;
      if (lmin[dim] >= dfg.getLevels()[dim]) continue;

      RemoteDataCollector<FG_ELEMENT> remoteData;
      if (dynamic_cast<HierarchicalHatBasisFunction*>(hierarchicalBases[dim]) != nullptr ||
          dynamic_cast<HierarchicalHatPeriodicBasisFunction*>(hierarchicalBases[dim]) != nullptr) {
        exchangeData1dDehierarchization(dfg, dim, remoteData, lmin[dim]);
      } else {
        exchangeAllData1d(dfg, dim, remoteData);
      }

      if (dfg.returnBoundaryFlags()[dim] > 0) {
        // sorry for the code duplication, could not figure out a clean way
        if (dynamic_cast<HierarchicalHatBasisFunction*>(hierarchicalBases[dim]) != nullptr) {
          hierarchizeWithBoundary<FG_ELEMENT, dehierarchize_hat_boundary_kernel<FG_ELEMENT>>(
              dfg, remoteData, dim, lmin[dim]);
        } else if (dynamic_cast<HierarchicalHatPeriodicBasisFunction*>(hierarchicalBases[dim]) !=
                   nullptr) {
          hierarchizeWithBoundary<FG_ELEMENT, dehierarchize_hat_boundary_kernel<FG_ELEMENT, true>>(
              dfg, remoteData, dim, lmin[dim]);
        } else if (dynamic_cast<FullWeightingBasisFunction*>(hierarchicalBases[dim]) != nullptr) {
          hierarchizeWithBoundary<FG_ELEMENT,
                                  dehierarchize_full_weighting_boundary_kernel<FG_ELEMENT, false>>(
              dfg, remoteData, dim, lmin[dim]);
        } else if (dynamic_cast<FullWeightingPeriodicBasisFunction*>(hierarchicalBases[dim]) !=
                   nullptr) {
          hierarchizeWithBoundary<FG_ELEMENT,
                                  dehierarchize_full_weighting_boundary_kernel<FG_ELEMENT, true>>(
              dfg, remoteData, dim, lmin[dim]);
        } else if (dynamic_cast<BiorthogonalBasisFunction*>(hierarchicalBases[dim]) != nullptr) {
          hierarchizeWithBoundary<FG_ELEMENT,
                                  dehierarchize_full_weighting_boundary_kernel<FG_ELEMENT, false>>(
              dfg, remoteData, dim, lmin[dim]);
        } else if (dynamic_cast<BiorthogonalPeriodicBasisFunction*>(hierarchicalBases[dim]) !=
                   nullptr) {
          hierarchizeWithBoundary<FG_ELEMENT,
                                  dehierarchize_biorthogonal_boundary_kernel<FG_ELEMENT, true>>(
              dfg, remoteData, dim, lmin[dim]);
        } else {
          throw std::logic_error("Not implemented");
        }
      } else {
        if (dynamic_cast<HierarchicalHatBasisFunction*>(hierarchicalBases[dim]) == nullptr) {
          throw std::logic_error("currently only hats supported for non-boundary grids");
        }
        assert(lmin[dim] == 0);
        dehierarchizeNoBoundary(dfg, remoteData, dim);
      }
      remoteData.clear();
    }
#pragma omp barrier
  }

  template <typename FG_ELEMENT, class T = HierarchicalHatBasisFunction>
  static void dehierarchizeHierachicalBasis(DistributedFullGrid<FG_ELEMENT>& dfg,
                                            const std::vector<bool>& dims,
                                            LevelVector lmin = LevelVector(0)) {
    T basisFctn;
    std::vector<BasisFunctionBasis*> bases(dfg.getDimension(), &basisFctn);
    if (lmin.size() == 0) {
      lmin = LevelVector(dfg.getDimension(), 0);
    }
    return dehierarchize<FG_ELEMENT>(dfg, dims, bases, lmin);
  }

  template <typename FG_ELEMENT>
  static void dehierarchize(DistributedFullGrid<FG_ELEMENT>& dfg,
                            LevelVector lmin = LevelVector(0)) {
    std::vector<bool> dims(dfg.getDimension(), true);
    if (lmin.size() == 0) {
      lmin = LevelVector(dfg.getDimension(), 0);
    }
    return dehierarchizeHierachicalBasis<FG_ELEMENT, HierarchicalHatBasisFunction>(dfg, dims, lmin);
  }

  template <typename FG_ELEMENT>
  static void dehierarchizeDFG(DistributedFullGrid<FG_ELEMENT>& dfg,
                               const std::vector<bool>& hierarchizationDims,
                               const std::vector<BasisFunctionBasis*>& hierarchicalBases,
                               LevelVector lmin = LevelVector(0)) {
    // dehierarchize dfg
    DistributedHierarchization::dehierarchize<FG_ELEMENT>(dfg, hierarchizationDims,
                                                          hierarchicalBases, lmin);
  }

  template <typename FG_ELEMENT>
  using FunctionPointer = void (*)(DistributedFullGrid<FG_ELEMENT>& dfg,
                                   const std::vector<bool>& dims, LevelVector lmin);
  // make template specifications visible by alias, hat is the default
  template <typename FG_ELEMENT>
  constexpr static FunctionPointer<FG_ELEMENT> hierarchizeHierarchicalHat =
      &hierarchizeHierachicalBasis<FG_ELEMENT, HierarchicalHatBasisFunction>;

  template <typename FG_ELEMENT>
  constexpr static FunctionPointer<FG_ELEMENT> hierarchizeHierarchicalHatPeriodic =
      &hierarchizeHierachicalBasis<FG_ELEMENT, HierarchicalHatPeriodicBasisFunction>;

  template <typename FG_ELEMENT>
  constexpr static FunctionPointer<FG_ELEMENT> hierarchizeFullWeighting =
      &hierarchizeHierachicalBasis<FG_ELEMENT, FullWeightingBasisFunction>;

  template <typename FG_ELEMENT>
  constexpr static FunctionPointer<FG_ELEMENT> hierarchizeFullWeightingPeriodic =
      &hierarchizeHierachicalBasis<FG_ELEMENT, FullWeightingPeriodicBasisFunction>;

  template <typename FG_ELEMENT>
  constexpr static FunctionPointer<FG_ELEMENT> hierarchizeBiorthogonal =
      &hierarchizeHierachicalBasis<FG_ELEMENT, BiorthogonalBasisFunction>;

  template <typename FG_ELEMENT>
  constexpr static FunctionPointer<FG_ELEMENT> hierarchizeBiorthogonalPeriodic =
      &hierarchizeHierachicalBasis<FG_ELEMENT, BiorthogonalPeriodicBasisFunction>;

  template <typename FG_ELEMENT>
  constexpr static FunctionPointer<FG_ELEMENT> dehierarchizeHierarchicalHat =
      &dehierarchizeHierachicalBasis<FG_ELEMENT, HierarchicalHatBasisFunction>;

  template <typename FG_ELEMENT>
  constexpr static FunctionPointer<FG_ELEMENT> dehierarchizeFullWeighting =
      &dehierarchizeHierachicalBasis<FG_ELEMENT, FullWeightingBasisFunction>;

  template <typename FG_ELEMENT>
  constexpr static FunctionPointer<FG_ELEMENT> dehierarchizeFullWeightingPeriodic =
      &dehierarchizeHierachicalBasis<FG_ELEMENT, FullWeightingPeriodicBasisFunction>;

  template <typename FG_ELEMENT>
  constexpr static FunctionPointer<FG_ELEMENT> dehierarchizeBiorthogonal =
      &dehierarchizeHierachicalBasis<FG_ELEMENT, BiorthogonalBasisFunction>;

  template <typename FG_ELEMENT>
  constexpr static FunctionPointer<FG_ELEMENT> dehierarchizeBiorthogonalPeriodic =
      &dehierarchizeHierachicalBasis<FG_ELEMENT, BiorthogonalPeriodicBasisFunction>;
};
// class DistributedHierarchization

}  // namespace combigrid

#endif /* DistributedHierarchization_HPP_ */
