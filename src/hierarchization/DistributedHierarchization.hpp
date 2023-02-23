#ifndef DISTRIBUTEDHIERARCHIZATION_HPP_
#define DISTRIBUTEDHIERARCHIZATION_HPP_

#include "boost/lexical_cast.hpp"
#include "fullgrid/DistributedFullGrid.hpp"
#include "utils/IndexVector.hpp"
#include "utils/PowerOfTwo.hpp"
#include "utils/Stats.hpp"

namespace combigrid {

/* The RemoteDataContainer is meant to store a (d-1)-dimensional block of a
 * d-dimensional DistributedFullGrid. The RemoteDataContainer is d-dimensional,
 * but has exactly one point in (at least) one dimension.
 */
template <typename FG_ELEMENT>
class RemoteDataContainer {
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
  RemoteDataContainer(const IndexVector& sizes, DimType dim1d, IndexType keyIndex) {
    assert(sizes.size() > 0);

    for (DimType i = 0; i < sizes.size(); ++i) {
      assert(sizes[i] > 0);
    }

    assert(dim1d < sizes.size());
    assert(sizes[dim1d] == 1);

    assert(keyIndex > -1);

    auto dim = static_cast<DimType>(sizes.size());
    dim1d_ = dim1d;
    index1d_ = keyIndex;

    // compute num of elements and offsets
    IndexType nrElements = 1;
    offsets_.resize(dim);

    for (DimType j = 0; j < dim; j++) {
      offsets_[j] = nrElements;
      nrElements = nrElements * sizes[j];
    }

    data_.resize(nrElements);

    /*
     std::cout << "created remote data container with "
     << "\n " << " size" << data_.size()
     << "data_ adress " << this->getData() << std::endl;
     */
  }

  inline FG_ELEMENT* getData(const IndexVector& localIndexVector) {
    static_assert(uniformDecomposition,
                  "this assumes uniform decomposition, so local index vectors are almost the same "
                  "along pole");
    auto dim = static_cast<DimType>(offsets_.size());
    assert(localIndexVector.size() == dim);

    // we have to find the corresponding local IndexVector of the
    // subdomain where the remoteData comes from and reduce it by the
    // key dimension

    IndexType idx = 0;

    for (DimType i = 0; i < dim; ++i) {
      // only use non-key dimensions
      if (i != dim1d_) {
        idx = idx + offsets_[i] * localIndexVector[i];
      }
    }

    assert(static_cast<size_t>(idx) < data_.size());

    return &data_[idx];
  }

  inline const FG_ELEMENT* getData(size_t idx) const {
    return &data_[idx];
  }
  /** the getters for the full grid vector */
  inline std::vector<FG_ELEMENT>& getElementVector() { return data_; }

  inline const std::vector<FG_ELEMENT>& getElementVector() const { return data_; }

  // return index of (d-1)-dimensional subgrid in the d-dimensional grid
  inline IndexType getKeyIndex() const { return index1d_; }

 private:
  // reduced dimension
  DimType dim1d_;

  // index of (d-1)-dimensional subgrid in the d-dimensional grid
  IndexType index1d_;

  // offsets in each dimension. d-dimensional
  IndexVector offsets_;

  // data vector
  std::vector<FG_ELEMENT> data_;
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
static IndexType getNextIndex1d(const DistributedFullGrid<FG_ELEMENT>& dfg, DimType d, IndexType idx1d);

template <typename FG_ELEMENT>
static IndexType getFirstIndexOfLevel1d(const DistributedFullGrid<FG_ELEMENT>& dfg, DimType d,
                                        LevelType l);

/**
 * @brief helper function for data exchange
 *
 * @param send1dIndices : a vector which holds (for each other rank) a set of own 1D indices to send to
 * @param recv1dIndices : a vector which holds (for each other rank) a set of 1D indices to receive from
 * @param dfg : the DistributedFullGrid where the own values are stored
 * @param dim : the dimension in which we want to exchange
 * @param remoteData : the data structure into which the received data will be stored
 */
template <typename FG_ELEMENT>
void sendAndReceiveIndices(const std::map<RankType, std::set<IndexType>>& send1dIndices,
                           const std::map<RankType, std::set<IndexType>>& recv1dIndices,
                           const DistributedFullGrid<FG_ELEMENT>& dfg, DimType dim,
                           std::vector<RemoteDataContainer<FG_ELEMENT>>& remoteData) {
  assert(remoteData.empty());

  std::vector<MPI_Request> sendRequests;
  std::vector<MPI_Request> recvRequests;
  size_t sendcount = 0;
  size_t recvcount = 0;
  for (const auto& x : send1dIndices) {
    sendcount += x.second.size();
  }
  for (const auto& x : recv1dIndices) {
    recvcount += x.second.size();
  }

  sendRequests.resize(sendcount);
  recvRequests.resize(recvcount);

  // for each rank r in send that has a nonempty index list
  sendcount = 0;
  for (const auto& x : send1dIndices) {
    const auto& r = x.first;
    const auto& indices = x.second;

    size_t k = 0;
    // for each index in index list
    for (const auto& index : indices) {
      // convert global 1d index i to local 1d index
      static IndexVector lidxvec(dfg.getDimension(), 0);
      lidxvec.resize(dfg.getDimension());
      {
        static IndexVector gidxvec;
        gidxvec = dfg.getLowerBounds();
        gidxvec[dim] = index;
        [[maybe_unused]] bool tmp = dfg.getLocalVectorIndex(gidxvec, lidxvec);
        assert(tmp && "index to be send not in local domain");
      }

      // create subarray view on the block with the local index
      MPI_Datatype mysubarray;
      {
        // sizes of local grid
        std::vector<int> sizes(dfg.getLocalSizes().begin(), dfg.getLocalSizes().end());

        // sizes of subarray ( full size except dimension d )
        std::vector<int> subsizes = sizes;
        subsizes[dim] = 1;

        // start
        std::vector<int> starts(dfg.getDimension(), 0);
        starts[dim] = static_cast<int>(lidxvec[dim]);

        // create subarray view on data
        MPI_Type_create_subarray(static_cast<int>(dfg.getDimension()), &sizes[0], &subsizes[0],
                                 &starts[0], MPI_ORDER_FORTRAN, dfg.getMPIDatatype(), &mysubarray);
        MPI_Type_commit(&mysubarray);
      }

      // send to rank r, use global index as tag
      {
        int dest = static_cast<int>(r);
        int tag = static_cast<int>(index);
        MPI_Isend(dfg.getData(), 1, mysubarray, dest, tag, dfg.getCommunicator(),
                &sendRequests[sendcount + k++]);

      }

      MPI_Type_free(&mysubarray);
    }
    sendcount += indices.size();
  }
  recvcount = 0;
  {
    // for each index in recv index list
    for (auto const& x : recv1dIndices) {
      auto const& r = x.first;
      // for each index i in index list
      const auto& indices = x.second;
      const IndexVector& lowerBoundsNeighbor = dfg.getLowerBounds(static_cast<int>(r));

      size_t k = 0;
      for (const auto& index : indices) {
        // create RemoteDataContainer to store the subarray
        IndexVector sizes = dfg.getLocalSizes();
        sizes[dim] = 1;
        remoteData.emplace_back(sizes, dim, index);

        // start recv operation, use global index as tag
        {
          int src = static_cast<int>(r);
          int tag = static_cast<int>(index);

          auto& buf = remoteData.back().getElementVector();
          auto bsize = static_cast<int>(buf.size());

          MPI_Irecv(buf.data(), bsize, dfg.getMPIDatatype(), src, tag, dfg.getCommunicator(),
                    &recvRequests[recvcount + k++]);

        }
      }
      recvcount += indices.size();
    }
  }
  // wait for finish of communication
  MPI_Waitall(static_cast<int>(sendRequests.size()), &sendRequests.front(), MPI_STATUSES_IGNORE);
  MPI_Waitall(static_cast<int>(recvRequests.size()), &recvRequests.front(), MPI_STATUSES_IGNORE);
}

/**
 * @brief same as sendAndReceiveIndices, but have only one MPI_Isend/Irecv per rank
 *
 * @note for small numbers of workers per grid, this seems to be less performant. may be different
 * for higher numbers.
 */
template <typename FG_ELEMENT>
void sendAndReceiveIndicesBlock(const std::map<RankType, std::set<IndexType>>& send1dIndices,
                                const std::map<RankType, std::set<IndexType>>& recv1dIndices,
                                const DistributedFullGrid<FG_ELEMENT>& dfg, DimType dim,
                                std::vector<RemoteDataContainer<FG_ELEMENT>>& remoteData) {
  assert(remoteData.empty());

  // count non-empty elements of input indices
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

  // for each rank r in send1dIndices that has a nonempty index list
  numSend = 0;
  for (const auto& x : send1dIndices) {
    const auto& r = x.first;
    const auto& indices = x.second;
    if (!indices.empty()) {
      // make datatype hblock for all indices
      MPI_Datatype myHBlock;
      std::vector<MPI_Aint> displacements;
      displacements.reserve(indices.size());
      for (const auto& index : indices) {
        // convert global 1d index i to local 1d index
        IndexType localLinearIndex = 0;
        {
          static IndexVector lidxvec(dfg.getDimension(), 0);
          lidxvec.resize(dfg.getDimension());
          static IndexVector gidxvec;
          gidxvec = dfg.getLowerBounds();
          gidxvec[dim] = index;
          [[maybe_unused]] bool tmp = dfg.getLocalVectorIndex(gidxvec, lidxvec);
          assert(tmp && "index to be send not in local domain");
          localLinearIndex = dfg.getLocalLinearIndex(lidxvec);
        }
        MPI_Aint addr;
        MPI_Get_address(&(dfg.getElementVector()[localLinearIndex]), &addr);
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
        MPI_Isend(dfg.getData(), 1, myHBlock, dest, tag, dfg.getCommunicator(),
                  &sendRequests[numSend++]);

      }
      // MPI_Type_free(&mysubarrayBlock);
      MPI_Type_free(&myHBlock);
    }
  }
  MPI_Type_free(&mysubarray);

  // for each rank r in recv1dIndices that has a nonempty index list
  numRecv = 0;
  for (const auto& x : recv1dIndices) {
    const auto& r = x.first;
    const auto& indices = x.second;
    if (!indices.empty()) {
      const IndexVector& lowerBoundsNeighbor = dfg.getLowerBounds(static_cast<int>(r));

      std::vector<FG_ELEMENT*> bufs;
      bufs.reserve(indices.size());
      IndexVector sizes = dfg.getLocalSizes();
      sizes[dim] = 1;
      int bsize = static_cast<int>(
          std::accumulate(sizes.begin(), sizes.end(), 1, std::multiplies<IndexType>()));

      for (const auto& index : indices) {
        // create RemoteDataContainer to store the subarray
        remoteData.emplace_back(sizes, dim, index);

        auto& buf = remoteData.back().getElementVector();
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

          MPI_Irecv(static_cast<void*>(bufs[0]), 1, myHBlock, src, tag, dfg.getCommunicator(),
                    &recvRequests[numRecv++]);
        }
        MPI_Type_free(&myHBlock);
      }
    }
  }
  assert(sendRequests.size() == numSend);
  assert(recvRequests.size() == numRecv);
  // wait for finish of communication
  MPI_Waitall(static_cast<int>(sendRequests.size()), &sendRequests.front(), MPI_STATUSES_IGNORE);
  MPI_Waitall(static_cast<int>(recvRequests.size()), &recvRequests.front(), MPI_STATUSES_IGNORE);
}

// exchange data in dimension dim

/**
 * @brief share all data with neighboring processes in one dimension
 *
 * @param dfg : the DistributedFullGrid where the own values are stored
 * @param dim : the dimension in which we want to exchange
 * @param remoteData : the data structure into which the received data will be stored
 */
template <typename FG_ELEMENT>
static void exchangeAllData1d(const DistributedFullGrid<FG_ELEMENT>& dfg, DimType dim,
                              std::vector<RemoteDataContainer<FG_ELEMENT>>& remoteData) {
    // send every index to all neighboring ranks in dimension dim
  auto globalIdxMax = dfg.length(dim);
  IndexType idxMin = dfg.getFirstGlobal1dIndex(dim);
  IndexType idxMax = dfg.getLastGlobal1dIndex(dim);

  auto poleNeighbors = dfg.getCartesianUtils().getAllMyPoleNeighborRanks(dim);

  std::set<IndexType> allMyIndices;
  for (IndexType i = idxMin; i <= idxMax; ++i) {
    allMyIndices.insert(i);
  }
  std::map<RankType, std::set<IndexType>> send1dIndices;
  std::map<RankType, std::set<IndexType>> recv1dIndices;

  for (auto& r : poleNeighbors) {
    send1dIndices[r] = allMyIndices;
    recv1dIndices[r] = {};
  }

  // all other points that are not ours can be received from their owners
  for (IndexType i = 0; i < idxMin; ++i) {
    // get rank which has i and add to recv list
    // TODO would be easier to iterate the whole range of each neighbor
    int r = dfg.getNeighbor1dFromAxisIndex(dim, i);
    if (r >= 0) recv1dIndices.at(r).insert(i);
  }
  for (IndexType i = idxMax + 1; i < globalIdxMax; ++i) {
    // get rank which has i and add to recv list
    int r = dfg.getNeighbor1dFromAxisIndex(dim, i);
    if (r >= 0) recv1dIndices.at(r).insert(i);
  }

  // sendAndReceiveIndices(send1dIndices, recv1dIndices, dfg, dim, remoteData);
  sendAndReceiveIndicesBlock(send1dIndices, recv1dIndices, dfg, dim, remoteData);
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
                           std::vector<RemoteDataContainer<FG_ELEMENT>>& remoteData,
                           LevelType lmin = 0) {

  // create buffers for every rank
  std::map<RankType, std::set<IndexType>> recv1dIndices;
  std::map<RankType, std::set<IndexType>> send1dIndices;

  // main loop
  IndexType idxMin = dfg.getFirstGlobal1dIndex(dim);
  IndexType idxMax = dfg.getLastGlobal1dIndex(dim);
  bool oneSidedBoundary = dfg.returnBoundaryFlags()[dim] == 1;
  auto globalIdxMax = dfg.length(dim);
  LevelType lmax = dfg.getLevels()[dim];

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

  // sendAndReceiveIndices(send1dIndices, recv1dIndices, dfg, dim, remoteData);
  sendAndReceiveIndicesBlock(send1dIndices, recv1dIndices, dfg, dim, remoteData);
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
static void exchangeData1dDehierarchization(
    const DistributedFullGrid<FG_ELEMENT>& dfg, DimType dim,
    std::vector<RemoteDataContainer<FG_ELEMENT>>& remoteData, LevelType lmin = 0) {

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

  // sendAndReceiveIndices(send1dIndices, recv1dIndices, dfg, dim, remoteData);
  sendAndReceiveIndicesBlock(send1dIndices, recv1dIndices, dfg, dim, remoteData);
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
inline void hierarchize_hat_boundary_kernel(FG_ELEMENT* data, LevelType lmax, int start,
                                             int stride, LevelType lmin = 0) {
  const int lmaxi = static_cast<int>(lmax);
  int ll = lmaxi - 1;
  int steps = (1 << (lmaxi - 1));
  int offset = 1;  // 1 and not 0 because boundary
  int stepsize = 2;
  int parentOffset = 1;

  for (; ll >= lmin; ll--) {
    int parOffsetStrided = parentOffset * stride;
    FG_ELEMENT parentL = 0.5 * data[start + offset * stride - parOffsetStrided];

    for (int ctr = 0; ctr < steps; ++ctr) {
      int centralIndex = start + offset * stride;
      FG_ELEMENT parentR = 0.5 * data[centralIndex + parOffsetStrided];
      FG_ELEMENT val1 = data[centralIndex];
      FG_ELEMENT val2 = val1 - parentL;
      FG_ELEMENT val3 = val2 - parentR;
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
                                                            int start, int stride,
                                                            LevelType lmin = 0) {
  assert(start == 0); // could be used but currently is not
  assert(stride == 1);
  const int lmaxi = static_cast<int>(lmax);
  int idxmax = powerOfTwo[lmaxi];
  // auto length = idxmax + 1;

  for (LevelType ldiff = 0; ldiff < lmax-lmin; ++ldiff) {
    int step_width = powerOfTwo[ldiff];
    // update f at even indices
    if (periodic) {
      // values at 0 and idxmax will be the same
      data[0] = 0.25 * (data[0] + data[idxmax] + data[step_width] + data[idxmax - step_width]);
      data[idxmax] = data[0];
    } else {
      data[0] = 0.5 * (data[0] + data[step_width]);
      data[idxmax] = 0.5 * (data[idxmax] + data[idxmax - step_width]);
    }
    //todo iterate only "our" part
    for (int i = 2*step_width; i < idxmax; i += 2*step_width) {
      // todo reformulate more cache-efficient
      data[i] = 0.25 * (data[i-step_width] + data[i+step_width]) + 0.5*data[i];
    }
    // update alpha / hierarchical surplus at odd indices
    for (int i = step_width; i < idxmax; i += 2*step_width) {
      // todo reformulate more cache-efficient
      data[i] = -0.5 * (data[i-step_width] + data[i+step_width]) + data[i];
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
 * @param start (unused)
 * @param stride (unused)
 * @param lmin minimum level (if > 0, hierarchization is not performed all the way down)
 */
template <typename FG_ELEMENT, bool periodic = false>
inline void hierarchize_biorthogonal_boundary_kernel(FG_ELEMENT* data, LevelType lmax, int start,
                                             int stride, LevelType lmin = 0) {

  assert(start == 0);
  assert(stride == 1);
  const int lmaxi = static_cast<int>(lmax);
  int idxmax = powerOfTwo[lmaxi];
  // auto length = idxmax + 1;
//     for l_hierarchization in range(l_max, l_min, -1):
//         step_width=2**(l_max-l_hierarchization)
//         for i in range(step_width,2**(l_max),2*step_width):
//             y[i]= -0.5*y[i-step_width] + y[i] -0.5*y[i+step_width]
//         y[0] = y[0] + 0.5*y[1*step_width]
//         y[-1] = y[-1] + 0.5*y[-1-1*step_width]
//         for i in range(2*step_width,2**(l_max),2*step_width):
//             y[i] = 0.25*y[i-step_width] + y[i] + 0.25*y[i+step_width]

  for (LevelType ldiff = 0; ldiff < lmax-lmin; ++ldiff) {
    int step_width = powerOfTwo[ldiff];
    // update alpha / hierarchical surplus at odd indices
    for (int i = step_width; i < idxmax; i += 2*step_width) {
      // todo reformulate more cache-efficient
      data[i] = -0.5 * (data[i-step_width] + data[i+step_width]) + data[i];
    }
    // update f at even indices
    if (periodic) {
      // values at 0 and idxmax will be the same
      data[0] = 0.5 * (data[0] + data[idxmax]) + 0.25 * (data[step_width] + data[idxmax - step_width]);
      data[idxmax] = data[0];
    } else {
      // mass will build up at the boundary; corresponds to 0-neumann-condition
      data[0] = data[0] + 0.5 * data[step_width];
      data[idxmax] = data[idxmax] + 0.5 * data[idxmax - step_width];
    }
    //todo iterate only "our" part
    for (int i = 2*step_width; i < idxmax; i += 2*step_width) {
      // todo reformulate more cache-efficient
      data[i] = 0.25 * (data[i-step_width] + data[i+step_width]) + data[i];
    }
  }
}

template <typename FG_ELEMENT, bool periodic = false>
inline void dehierarchize_hat_boundary_kernel(FG_ELEMENT* data, LevelType lmax, int start,
                                               int stride, LevelType lmin = 0) {
  const int lmaxi = static_cast<int>(lmax);
  const int lmini = static_cast<int>(lmin);
  int steps = 1 << (lmini);
  int offset = (1 << (lmaxi - lmini - 1));
  int stepsize = (1 << (lmaxi - lmini));
  int parentOffset = offset;

  for (LevelType ll = lmin + 1; ll <= lmax; ++ll) {
    int parOffsetStrided = parentOffset * stride;
    FG_ELEMENT parentL = 0.5 * data[start + offset * stride - parOffsetStrided];

    for (int ctr = 0; ctr < steps; ++ctr) {
      int centralIndex = start + offset * stride;
      FG_ELEMENT parentR = 0.5 * data[centralIndex + parOffsetStrided];
      FG_ELEMENT val1 = data[centralIndex];
      FG_ELEMENT val2 = val1 + parentL;
      FG_ELEMENT val3 = val2 + parentR;
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
                                                          int start, int stride, LevelType lmin) {
  const int lmaxi = static_cast<int>(lmax);
  int idxmax = powerOfTwo[lmaxi];
  // auto length = idxmax + 1;

  for (auto ldiff = static_cast<LevelType>(lmax - lmin - 1); ldiff >= 0; --ldiff) {
    int step_width = powerOfTwo[ldiff];
    // update alpha / hierarchical surplus at odd indices
    for (int i = step_width; i < idxmax; i += 2 * step_width) {
      // todo reformulate more cache-efficient
      data[i] = 0.5 * (data[i - step_width] + data[i + step_width]) + data[i];
    }
    // update f at even indices
    for (int i = 2 * step_width; i < idxmax; i += 2 * step_width) {
      // todo reformulate more cache-efficient
      data[i] = -0.5 * (data[i - step_width] + data[i + step_width]) + 2. * data[i];
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
inline void dehierarchize_biorthogonal_boundary_kernel(FG_ELEMENT* data, LevelType lmax, int start,
                                                        int stride, LevelType lmin) {
  assert(start == 0);
  assert(stride == 1);
  const int lmaxi = static_cast<int>(lmax);
  int idxmax = powerOfTwo[lmaxi];
  // auto length = idxmax + 1;

  // for l_hierarchization in range(l_min, l_max, 1):
  //     step_width=2**(l_max-l_hierarchization-1)
  //     y[0] = y[0] - 0.5*y[0+step_width]
  //     y[-1] = y[-1] - 0.5*y[-1-1*step_width]
  //     for i in range(2*step_width,2**(l_max),2*step_width):
  //         y[i] = -0.25*y[i-step_width] + 1*y[i] - 0.25*y[i+step_width]
  //     for i in range(step_width,2**(l_max),2*step_width):
  //         y[i]= 0.5*y[i-step_width] + y[i] + 0.5*y[i+step_width]

  for (auto ldiff = static_cast<LevelType>(lmax - lmin - 1); ldiff >= 0; --ldiff) {
    int step_width = powerOfTwo[ldiff];
    // update f at even indices
    if (periodic) {
      data[0] =
          -0.25 * (data[step_width] + data[idxmax - step_width]) + 0.5 * (data[0] + data[idxmax]);
      data[idxmax] = data[0];
    } else {
      data[0] = data[0] - 0.5 * data[step_width];
      data[idxmax] = data[idxmax] - 0.5 * data[idxmax - step_width];
    }
    for (int i = 2 * step_width; i < idxmax; i += 2 * step_width) {
      // todo reformulate more cache-efficient
      data[i] = -0.25 * (data[i - step_width] + data[i + step_width]) + data[i];
    }
    // update alpha / hierarchical surplus at odd indices
    for (int i = step_width; i < idxmax; i += 2 * step_width) {
      // todo reformulate more cache-efficient
      data[i] = 0.5 * (data[i - step_width] + data[i + step_width]) + data[i];
    }
  }
}

/**
 * @brief  hierarchize a DFG with boundary points in dimension dim
 */
template <typename FG_ELEMENT,
          void (*HIERARCHIZATION_FCTN)(FG_ELEMENT[], LevelType, int, int,
                                       LevelType) = hierarchize_hat_boundary_kernel>
void hierarchizeWithBoundary(DistributedFullGrid<FG_ELEMENT>& dfg,
                             std::vector<RemoteDataContainer<FG_ELEMENT>>& remoteData, DimType dim,
                             LevelType lmin_n = 0) {
  assert(dfg.returnBoundaryFlags()[dim] > 0);

  auto lmax = dfg.getLevels()[dim];
  auto size = dfg.getNrLocalElements();
  auto stride = dfg.getLocalOffsets()[dim];
  auto ndim = dfg.getLocalSizes()[dim];
  IndexType jump = stride * ndim;
  IndexType nbrOfPoles = size / ndim;

  static IndexVector localIndexVector(dfg.getDimension());
  localIndexVector.resize(dfg.getDimension());

  // loop over poles
  static std::vector<FG_ELEMENT> tmp(dfg.getGlobalSizes()[dim],
                                     std::numeric_limits<double>::quiet_NaN());
  // if we are using periodicity, add an entry to tmp for the virtual last value
  bool oneSidedBoundary = dfg.returnBoundaryFlags()[dim] == 1;
  tmp.resize(dfg.getGlobalSizes()[dim] + (oneSidedBoundary ? 1 : 0),
             std::numeric_limits<double>::quiet_NaN());
  std::vector<FG_ELEMENT>& ldata = dfg.getElementVector();
  lldiv_t divresult;
  IndexType start;
  IndexType gstart = dfg.getLowerBounds()[dim];

  for (IndexType nn = 0; nn < nbrOfPoles;
       ++nn) {  // integer operations form bottleneck here -- nested loops are twice as slow
    divresult = std::lldiv(nn, stride);
    start = divresult.quot * jump + divresult.rem;  // localer lin index start of pole

#ifndef NDEBUG
    // compute global vector index of start
    dfg.getLocalVectorIndex(start, localIndexVector);
    assert(localIndexVector[dim] == 0);
    for (size_t i = 0; i < remoteData.size(); ++i) {
      assert(*remoteData[i].getData(localIndexVector) == *remoteData[i].getData(nn));
    }
#endif // NDEBUG

    // go through remote containers
    for (size_t i = 0; i < remoteData.size(); ++i) {
      tmp[remoteData[i].getKeyIndex()] = *remoteData[i].getData(nn);
    }

    // copy local data
    for (IndexType i = 0; i < ndim; ++i) tmp[gstart + i] = ldata[start + stride * i];

    if (oneSidedBoundary) {
      // assume periodicity
      //  assert(HIERARCHIZATION_FCTN::periodic); //TODO
      tmp[dfg.getGlobalSizes()[dim]] = tmp[0];
      if (!remoteData.empty() && remoteData[0].getKeyIndex() == 0) {
        assert(!std::isnan(std::real(tmp[0])));
      }
    }
    // hierarchize tmp array with hupp function
    HIERARCHIZATION_FCTN(&tmp[0], lmax, 0, 1, lmin_n);

    if (oneSidedBoundary && !remoteData.empty() && remoteData[0].getKeyIndex() == 0) {
      assert(!std::isnan(std::real(tmp[0])));
      assert(tmp[dfg.getGlobalSizes()[dim]] == tmp[0]);
    }

    // copy pole back
    for (IndexType i = 0; i < ndim; ++i) {
      ldata[start + stride * i] = tmp[gstart + i];
      assert(!std::isnan(std::real(tmp[gstart + i])));
    }
  }
}

/**
 * @brief  hierarchize a DFG without boundary points in dimension dim
 */
template <typename FG_ELEMENT>
void hierarchizeNoBoundary(DistributedFullGrid<FG_ELEMENT>& dfg,
                           std::vector<RemoteDataContainer<FG_ELEMENT>>& remoteData, DimType dim) {
  assert(dfg.returnBoundaryFlags()[dim] == 0);

  LevelType lmax = dfg.getLevels()[dim];
  IndexType size = dfg.getNrLocalElements();
  IndexType stride = dfg.getLocalOffsets()[dim];
  IndexType ndim = dfg.getLocalSizes()[dim];
  IndexType jump = stride * ndim;
  IndexType nbrOfPoles = size / ndim;

  IndexVector localIndexVector(dfg.getDimension());

  // loop over poles
  std::vector<FG_ELEMENT> tmp(dfg.getGlobalSizes()[dim], std::numeric_limits<double>::quiet_NaN());
  std::vector<FG_ELEMENT>& ldata = dfg.getElementVector();
  lldiv_t divresult;
  IndexType start;
  IndexType gstart = dfg.getLowerBounds()[dim];

  for (IndexType nn = 0; nn < nbrOfPoles;
       ++nn) {  // integer operations form bottleneck here -- nested loops are twice as slow
    divresult = std::lldiv(nn, stride);
    start = divresult.quot * jump + divresult.rem;  // localer lin index start of pole

    // compute global vector index of start
    dfg.getLocalVectorIndex(start, localIndexVector);
    assert(localIndexVector[dim] == 0);

    // go through remote containers
    for (size_t i = 0; i < remoteData.size(); ++i) {
      tmp[remoteData[i].getKeyIndex()] = *remoteData[i].getData(localIndexVector);
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
 * @brief inverse operation for hierarchizeWithBoundary
 */
template <typename FG_ELEMENT,
          void (*DEHIERARCHIZATION_FCTN)(FG_ELEMENT[], LevelType, int, int,
                                         LevelType) = dehierarchize_hat_boundary_kernel>
void dehierarchizeWithBoundary(DistributedFullGrid<FG_ELEMENT>& dfg,
                               std::vector<RemoteDataContainer<FG_ELEMENT>>& remoteData,
                               DimType dim, LevelType lmin_n = 0) {
  assert(dfg.returnBoundaryFlags()[dim] > 0);

  const auto& lmax = dfg.getLevels()[dim];
  const auto& size = dfg.getNrLocalElements();
  const auto& stride = dfg.getLocalOffsets()[dim];
  const auto& ndim = dfg.getLocalSizes()[dim];
  IndexType jump = stride * ndim;
  IndexType nbrOfPoles = size / ndim;

  static IndexVector localIndexVector(dfg.getDimension());
  localIndexVector.resize(dfg.getDimension());

  // loop over poles
  static std::vector<FG_ELEMENT> tmp(dfg.getGlobalSizes()[dim],
                                     std::numeric_limits<double>::quiet_NaN());
  // if we are using periodicity, add an entry to tmp for the virtual last value
  bool oneSidedBoundary = dfg.returnBoundaryFlags()[dim] == 1;
  tmp.resize(dfg.getGlobalSizes()[dim] + (oneSidedBoundary ? 1 : 0),
             std::numeric_limits<double>::quiet_NaN());
  std::vector<FG_ELEMENT>& ldata = dfg.getElementVector();
  lldiv_t divresult;
  IndexType start;
  IndexType gstart = dfg.getLowerBounds()[dim];

  for (IndexType nn = 0; nn < nbrOfPoles;
       ++nn) {  // integer operations form bottleneck here -- nested loops are twice as slow
    divresult = std::lldiv(nn, stride);
    start = divresult.quot * jump + divresult.rem;  // localer lin index start of pole

#ifndef NDEBUG
    // compute global vector index of start
    dfg.getLocalVectorIndex(start, localIndexVector);
    assert(localIndexVector[dim] == 0);
    for (size_t i = 0; i < remoteData.size(); ++i) {
      assert(*remoteData[i].getData(localIndexVector) == *remoteData[i].getData(nn));
    }
#endif // NDEBUG

    // go through remote containers
    for (size_t i = 0; i < remoteData.size(); ++i) {
      tmp[remoteData[i].getKeyIndex()] = *remoteData[i].getData(nn);
    }

    // copy local data
    for (IndexType i = 0; i < ndim; ++i) tmp[gstart + i] = ldata[start + stride * i];

    if (oneSidedBoundary) {
      // assume periodicity
      tmp[dfg.getGlobalSizes()[dim]] = tmp[0];
    }
    // hierarchize tmp array with hupp function
    DEHIERARCHIZATION_FCTN(&tmp[0], lmax, 0, 1, lmin_n);

    if (oneSidedBoundary && !remoteData.empty() && remoteData[0].getKeyIndex() == 0) {
      assert(tmp[dfg.getGlobalSizes()[dim]] == tmp[0]);
    }

    // copy pole back
    for (IndexType i = 0; i < ndim; ++i) {
      ldata[start + stride * i] = tmp[gstart + i];
      assert(!std::isnan(std::real(tmp[gstart + i])));
    }
  }
}

/**
 * @brief inverse operation for hierarchizeNoBoundary
 */
template <typename FG_ELEMENT>
void dehierarchizeNoBoundary(DistributedFullGrid<FG_ELEMENT>& dfg,
                             std::vector<RemoteDataContainer<FG_ELEMENT>>& remoteData,
                             DimType dim) {
  assert(dfg.returnBoundaryFlags()[dim] == 0);

  auto lmax = dfg.getLevels()[dim];
  auto size = dfg.getNrLocalElements();
  auto stride = dfg.getLocalOffsets()[dim];
  auto ndim = dfg.getLocalSizes()[dim];
  IndexType jump = stride * ndim;
  IndexType nbrOfPoles = size / ndim;

  IndexVector localIndexVector(dfg.getDimension());

  // loop over poles
  static std::vector<FG_ELEMENT> tmp;
  tmp.resize(dfg.getGlobalSizes()[dim], std::numeric_limits<double>::quiet_NaN());
  std::vector<FG_ELEMENT>& ldata = dfg.getElementVector();
  lldiv_t divresult;
  IndexType start;
  IndexType gstart = dfg.getLowerBounds()[dim];

  for (IndexType nn = 0; nn < nbrOfPoles;
       ++nn) {  // integer operations form bottleneck here -- nested loops are twice as slow
    divresult = std::lldiv(nn, stride);
    start = divresult.quot * jump + divresult.rem;  // localer lin index start of pole

    // compute global vector index of start
    dfg.getLocalVectorIndex(start, localIndexVector);
    assert(localIndexVector[dim] == 0);

    // go through remote containers
    for (size_t i = 0; i < remoteData.size(); ++i) {
      tmp[remoteData[i].getKeyIndex()] = *remoteData[i].getData(localIndexVector);
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
                          LevelVector lmin = LevelVector(0)) {
    assert(dfg.getDimension() > 0);
    assert(dfg.getDimension() == dims.size());
    if (lmin.size() == 0) {
      lmin = LevelVector(dims.size(), 0);
    }
    // hierarchize all dimensions, with special treatment for 0
    for (DimType dim = 0; dim < dfg.getDimension(); ++dim) {
      if (!dims[dim]) continue;

      // exchange data
      std::vector<RemoteDataContainer<FG_ELEMENT>> remoteData;
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
    }
  }

  template <typename FG_ELEMENT, class T = HierarchicalHatBasisFunction>
  static void hierarchizeHierachicalBasis(DistributedFullGrid<FG_ELEMENT>& dfg,
                                          const std::vector<bool>& dims,
                                          LevelVector lmin = LevelVector(0)) {
    T basisFctn;
    std::vector<BasisFunctionBasis*> bases(dfg.getDimension(), &basisFctn);
    return hierarchize<FG_ELEMENT>(dfg, dims, bases, lmin);
  }

  template <typename FG_ELEMENT>
  static void hierarchize(DistributedFullGrid<FG_ELEMENT>& dfg, LevelVector lmin = LevelVector(0)) {
    std::vector<bool> dims(dfg.getDimension(), true);
    return hierarchizeHierachicalBasis<FG_ELEMENT, HierarchicalHatBasisFunction>(dfg, dims, lmin);
  }

  // inplace dehierarchization
  template <typename FG_ELEMENT>
  static void dehierarchize(DistributedFullGrid<FG_ELEMENT>& dfg, const std::vector<bool>& dims,
                            const std::vector<BasisFunctionBasis*>& hierarchicalBases,
                            LevelVector lmin = LevelVector(0)) {
    if (lmin.size() == 0) {
      lmin = LevelVector(dims.size(), 0);
    }
    assert(dfg.getDimension() > 0);
    assert(dfg.getDimension() == dims.size());
    // dehierarchize all dimensions, with special treatment for 0
    for (DimType dim = 0; dim < dfg.getDimension(); ++dim) {
      if (!dims[dim]) continue;

      // exchange data
      std::vector<RemoteDataContainer<FG_ELEMENT>> remoteData;
      if (dynamic_cast<HierarchicalHatBasisFunction*>(hierarchicalBases[dim]) != nullptr ||
          dynamic_cast<HierarchicalHatPeriodicBasisFunction*>(hierarchicalBases[dim]) != nullptr) {
        exchangeData1dDehierarchization(dfg, dim, remoteData, lmin[dim]);
      } else {
        exchangeAllData1d(dfg, dim, remoteData);
      }

      if (dfg.returnBoundaryFlags()[dim] > 0) {
        // sorry for the code duplication, could not figure out a clean way
        if (dynamic_cast<HierarchicalHatBasisFunction*>(hierarchicalBases[dim]) != nullptr) {
          dehierarchizeWithBoundary<FG_ELEMENT, dehierarchize_hat_boundary_kernel<FG_ELEMENT>>(
              dfg, remoteData, dim, lmin[dim]);
        } else if (dynamic_cast<HierarchicalHatPeriodicBasisFunction*>(hierarchicalBases[dim]) !=
                   nullptr) {
          dehierarchizeWithBoundary<FG_ELEMENT,
                                      dehierarchize_hat_boundary_kernel<FG_ELEMENT, true>>(
              dfg, remoteData, dim, lmin[dim]);
        } else if (dynamic_cast<FullWeightingBasisFunction*>(hierarchicalBases[dim]) != nullptr) {
          dehierarchizeWithBoundary<
              FG_ELEMENT, dehierarchize_full_weighting_boundary_kernel<FG_ELEMENT, false>>(
              dfg, remoteData, dim, lmin[dim]);
        } else if (dynamic_cast<FullWeightingPeriodicBasisFunction*>(hierarchicalBases[dim]) !=
                   nullptr) {
          dehierarchizeWithBoundary<
              FG_ELEMENT, dehierarchize_full_weighting_boundary_kernel<FG_ELEMENT, true>>(
              dfg, remoteData, dim, lmin[dim]);
        } else if (dynamic_cast<BiorthogonalBasisFunction*>(hierarchicalBases[dim]) != nullptr) {
          dehierarchizeWithBoundary<
              FG_ELEMENT, dehierarchize_full_weighting_boundary_kernel<FG_ELEMENT, false>>(
              dfg, remoteData, dim, lmin[dim]);
        } else if (dynamic_cast<BiorthogonalPeriodicBasisFunction*>(hierarchicalBases[dim]) !=
                   nullptr) {
          dehierarchizeWithBoundary<
              FG_ELEMENT, dehierarchize_biorthogonal_boundary_kernel<FG_ELEMENT, true>>(
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
    }
  }

  template <typename FG_ELEMENT, class T = HierarchicalHatBasisFunction>
  static void dehierarchizeHierachicalBasis(DistributedFullGrid<FG_ELEMENT>& dfg,
                                            const std::vector<bool>& dims,
                                            LevelVector lmin = LevelVector(0)) {
    T basisFctn;
    std::vector<BasisFunctionBasis*> bases(dfg.getDimension(), &basisFctn);
    return dehierarchize<FG_ELEMENT>(dfg, dims, bases, lmin);
  }

  template <typename FG_ELEMENT>
  static void dehierarchize(DistributedFullGrid<FG_ELEMENT>& dfg,
                            LevelVector lmin = LevelVector(0)) {
    std::vector<bool> dims(dfg.getDimension(), true);
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
