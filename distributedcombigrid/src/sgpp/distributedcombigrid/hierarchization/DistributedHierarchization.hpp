#ifndef DISTRIBUTEDHIERARCHIZATION_HPP_
#define DISTRIBUTEDHIERARCHIZATION_HPP_

//#define DEBUG_OUTPUT

#include "boost/lexical_cast.hpp"
#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGrid.hpp"
#include "sgpp/distributedcombigrid/legacy/combigrid_utils.hpp"
#include "sgpp/distributedcombigrid/utils/IndexVector.hpp"
#include "sgpp/distributedcombigrid/utils/Stats.hpp"

using namespace combigrid;

/*
 * Instead of having private static functions, I put these functions in an
 * unnamed namespace. So, they are not accessible from outside the file, as well.
 * In the general case, this would have the advantage, that we can change
 * the declaration of these functions without changing the declaration of the
 * class. So we avoid recompilation of all files that use the class.
 */
namespace {

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
  RemoteDataContainer(const IndexVector& sizes, DimType dim1d, IndexType keyIndex,
                      const IndexVector& lowerBounds) {
    assert(sizes.size() > 0);
    assert(lowerBounds.size() == sizes.size());

    for (DimType i = 0; i < sizes.size(); ++i) {
      assert(sizes[i] > 0);
    }

    assert(dim1d < sizes.size());
    assert(sizes[dim1d] == 1);

    assert(keyIndex > -1);

    dim_ = static_cast<DimType>(sizes.size());
    dim1d_ = dim1d;
    index1d_ = keyIndex;
    nrPoints_ = sizes;
    lowerBounds_ = lowerBounds;

    // compute num of elements and offsets
    nrElements_ = 1;
    offsets_.resize(dim_);

    for (DimType j = 0; j < dim_; j++) {
      offsets_[j] = nrElements_;
      nrElements_ = nrElements_ * nrPoints_[j];
    }

    data_.resize(nrElements_);

    /*
     std::cout << "created remote data container with "
     << "\n " << " size" << data_.size()
     << "data_ adress " << this->getData() << std::endl;
     */
  }

  inline FG_ELEMENT* getData(const IndexVector& globalIndexVector) {
    assert(globalIndexVector.size() == dim_);

    // we have to find the corresponding local IndexVector of the
    // subdomain where the remoteData comes from and reduce it by the
    // key dimension

    // compute local index in remote domain
    IndexVector localIndexVector = globalIndexVector - lowerBounds_;

    // reduce by key dimension
    localIndexVector[dim1d_] = 0;

    for (DimType i = 0; i < dim_; ++i) {
      assert(localIndexVector[i] < nrPoints_[i]);
    }

    IndexType idx = 0;

    for (DimType i = 0; i < dim_; ++i) {
      idx = idx + offsets_[i] * localIndexVector[i];
    }

    assert(idx < nrElements_);

    return &data_[idx];
  }

  inline IndexType get1dIndex(const IndexVector& globalIndexVector) const {
    assert(globalIndexVector.size() == dim_);

    // we have to find the corresponding local IndexVector of the
    // subdomain where the remoteData comes from and reduce it by the
    // key dimension

    // compute local index in remote domain
    IndexVector localIndexVector = globalIndexVector - lowerBounds_;

    // reduce by key dimension
    localIndexVector[dim1d_] = 0;

    for (DimType i = 0; i < dim_; ++i) {
      assert(localIndexVector[i] < nrPoints_[i]);
    }

    IndexType idx = 0;

    for (DimType i = 0; i < dim_; ++i) {
      idx = idx + offsets_[i] * localIndexVector[i];
    }

    assert(idx < nrElements_);

    return idx;
  }

  inline FG_ELEMENT* getData() { return &data_[0]; }

  /** the getters for the full grid vector */
  inline std::vector<FG_ELEMENT>& getElementVector() { return data_; }

  inline const std::vector<FG_ELEMENT>& getElementVector() const { return data_; }

  // return index of (d-1)-dimensional subgrid in the d-dimensional grid
  inline IndexType getKeyIndex() const { return index1d_; }

  inline IndexType getSize() { return data_.size(); }

  inline DimType getDimension() { return dim_; }

  inline DimType getKeyDimension() { return dim1d_; }

 private:
  // dimensionality of the container. although only one point in dim1d, we
  // always use the full dimensionality
  DimType dim_;

  // reduced dimension
  DimType dim1d_;

  // index of (d-1)-dimensional subgrid in the d-dimensional grid
  IndexType index1d_;

  // total number of points
  IndexType nrElements_;

  // nr of points in each dimension (d-dimensional)
  IndexVector nrPoints_;

  // offsets in each dimension. d-dimensional, but 0 in dim1d dimension
  IndexVector offsets_;

  // data vector
  std::vector<FG_ELEMENT> data_;

  // lower bounds of remote domain
  IndexVector lowerBounds_;
};

/**
 * Lookup table that hides the complexity of having both, local and remote data
 * via a common interface for the data access
 */
template <typename FG_ELEMENT>
class LookupTable {
 public:
  /** Constructor
   *
   * \param[in] remoteData  list with remote data
   * \param[in] dfg         local view of distributed fullgrid
   * \param[in] keyDim      dimension to which the (d-1)-dimensional subgrids
   *                        stored in remoteData have been reduced
   */
  LookupTable(std::vector<RemoteDataContainer<FG_ELEMENT> >& remoteData,
              DistributedFullGrid<FG_ELEMENT>& dfg, DimType keyDim)
      : remoteData_(remoteData), dfg_(dfg), keyDim_(keyDim) {
    if (remoteData_.size() > 0) {
      for (size_t i = 0; i < remoteData.size(); ++i) {
        assert(remoteData_[i].getDimension() == dfg_.getDimension());
        assert(remoteData_[i].getKeyDimension() == keyDim_);
      }
    }
  }

  inline FG_ELEMENT* getData(IndexVector globalIndexVector) {
    assert(globalIndexVector.size() == dfg_.getDimension());

    // check if in local part of the distributed full grid
    if (globalIndexVector >= dfg_.getLowerBounds() && globalIndexVector < dfg_.getUpperBounds()) {
      // return point to value in dfg
      static IndexVector localIndexVector(dfg_.getDimension());
      dfg_.getLocalVectorIndex(globalIndexVector, localIndexVector);

      IndexType localLinearIndex = dfg_.getLocalLinearIndex(localIndexVector);

      return &dfg_.getData()[localLinearIndex];
    } else {
      // find subarray remote data wich corresponds to key index
      bool found = false;

      for (size_t i = 0; i < remoteData_.size(); ++i) {
        DimType keyDim = remoteData_[i].getKeyDimension();

        if (remoteData_[i].getKeyIndex() == globalIndexVector[keyDim]) {
          // translate globalIndexVector to IndexVector for remote Data
          return remoteData_[i].getData(globalIndexVector);
        }
      }

      assert(found && "subarray not found in remote data");
    }
  }

  inline std::vector<RemoteDataContainer<FG_ELEMENT> >& getRDCVector() const { return remoteData_; }

 private:
  std::vector<RemoteDataContainer<FG_ELEMENT> >& remoteData_;

  DistributedFullGrid<FG_ELEMENT>& dfg_;

  DimType keyDim_;
};

template <typename FG_ELEMENT>
static void checkLeftSuccesors(IndexType checkIdx, IndexType rootIdx, DimType dim,
                               const DistributedFullGrid<FG_ELEMENT>& dfg,
                               std::map<RankType, std::set<IndexType>>& OneDIndices);

template <typename FG_ELEMENT>
static void checkRightSuccesors(IndexType checkIdx, IndexType rootIdx, DimType dim,
                                const DistributedFullGrid<FG_ELEMENT>& dfg,
                                std::map<RankType, std::set<IndexType>>& OneDIndices);

template <typename FG_ELEMENT>
static IndexType checkPredecessors(IndexType idx, DimType dim, const DistributedFullGrid<FG_ELEMENT>& dfg,
                                   std::map<RankType, std::set<IndexType>>& OneDIndices);

template <typename FG_ELEMENT>
static RankType getNeighbor1d(const DistributedFullGrid<FG_ELEMENT>& dfg, DimType d, IndexType idx1d);

template <typename FG_ELEMENT>
static IndexType getNextIndex1d(const DistributedFullGrid<FG_ELEMENT>& dfg, DimType d, IndexType idx1d);

template <typename FG_ELEMENT>
static IndexType getFirstIndexOfLevel1d(const DistributedFullGrid<FG_ELEMENT>& dfg, DimType d,
                                        LevelType l);
template <typename FG_ELEMENT>
static IndexType getLastIndexOfLevel1d(const DistributedFullGrid<FG_ELEMENT>& dfg, DimType d,
                                        LevelType l);

template <typename FG_ELEMENT>
static IndexVector getFirstIndexOfEachLevel1d(const DistributedFullGrid<FG_ELEMENT>& dfg, DimType d) {
  LevelType lmax = dfg.getLevels()[d];
  IndexVector firstIndices(lmax, -1);
  for (LevelType lidx = 0; lidx <= lmax; ++lidx) {
    // leftmost point of this level
    // currently leaving out level 0
    firstIndices[lidx] = getFirstIndexOfLevel1d(dfg, d, lidx + 1);
  }
  return firstIndices;
}

template <typename FG_ELEMENT>
static IndexVector getLastIndicesOfLevel1d(DistributedFullGrid<FG_ELEMENT>& dfg, DimType d);

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
#ifdef DEBUG_OUTPUT
  auto commSize = dfg.getCommunicatorSize();
  CommunicatorType comm = dfg.getCommunicator();
  auto rank = dfg.getRank();
  MPI_Barrier(comm);

  // print recvindices
  {
    for (int r = 0; r < commSize; ++r) {
      if (r == rank) {
        std::cout << "rank " << r << " recv1dIndices ";
        for (const auto& r : recv1dIndices) {
          std::cout << r.second;
        }
        std::cout << std::endl;
      }
      MPI_Barrier(comm);
    }
  }

  if (rank == 0) std::cout << std::endl;

  // print sendindices
  {
    for (int r = 0; r < commSize; ++r) {
      if (r == rank) {
        std::cout << "rank " << r << " send1dIndices ";
        for (const auto& s : send1dIndices) {
          std::cout << s.second;
        }
        std::cout << std::endl;
      }
      MPI_Barrier(comm);
    }
  }

  MPI_Barrier(comm);

#endif

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
      IndexVector lidxvec(dfg.getDimension(), 0);
      {
        IndexVector gidxvec = dfg.getLowerBounds();
        gidxvec[dim] = index;
        bool tmp = dfg.getLocalVectorIndex(gidxvec, lidxvec);
        assert(tmp && "index to be send not in local domain");
      }

      // create subarray view on the block with the local index
      MPI_Datatype mysubarray;
      {
        // sizes of local grid
        IndexVector sizes(dfg.getLocalSizes().begin(), dfg.getLocalSizes().end());

        // sizes of subarray ( full size except dimension d )
        IndexVector subsizes = sizes;
        subsizes[dim] = 1;

        // start
        IndexVector starts(dfg.getDimension(), 0);
        starts[dim] = lidxvec[dim];

        std::vector<int> csizes(sizes.begin(), sizes.end());
        std::vector<int> csubsizes(subsizes.begin(), subsizes.end());
        std::vector<int> cstarts(starts.begin(), starts.end());

        // create subarray view on data
        MPI_Type_create_subarray(static_cast<int>(dfg.getDimension()), &csizes[0], &csubsizes[0],
                                 &cstarts[0], MPI_ORDER_FORTRAN, dfg.getMPIDatatype(), &mysubarray);
        MPI_Type_commit(&mysubarray);
      }

      // send to rank r, use global index as tag
      {
        int dest = static_cast<int>(r);
        int tag = static_cast<int>(index);
        MPI_Isend(dfg.getData(), 1, mysubarray, dest, tag, dfg.getCommunicator(),
                &sendRequests[sendcount + k++]);

#ifdef DEBUG_OUTPUT
        auto rank = dfg.getRank();
        // print info: dest, size, index
        std::cout << "rank " << rank << ": send gindex " << index << " dest " << dest
                  << std::endl;
#endif
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
        remoteData.emplace_back(sizes, dim, index, lowerBoundsNeighbor);

        // start recv operation, use global index as tag
        {
          int src = static_cast<int>(r);
          int tag = static_cast<int>(index);

          FG_ELEMENT* buf = remoteData.back().getData();
          int bsize = static_cast<int>(remoteData.back().getSize());

          MPI_Irecv(buf, bsize, dfg.getMPIDatatype(), src, tag, dfg.getCommunicator(),
                    &recvRequests[recvcount + k++]);

#ifdef DEBUG_OUTPUT
          auto rank = dfg.getRank();
          // print info: dest, size, index
          std::cout << "rank " << rank << ": recv gindex " << index << " src " << src
                    << " size: " << bsize << std::endl;
#endif
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
#ifdef DEBUG_OUTPUT
  auto commSize = dfg.getCommunicatorSize();
  CommunicatorType comm = dfg.getCommunicator();
  auto rank = dfg.getRank();
  MPI_Barrier(comm);

  // print recvindices
  {
    for (int r = 0; r < commSize; ++r) {
      if (r == rank) {
        std::cout << "rank " << r << " recv1dIndices ";
        for (const auto& r : recv1dIndices) {
          std::cout << r.second;
        }
        std::cout << std::endl;
      }
      MPI_Barrier(comm);
    }
  }

  if (rank == 0) std::cout << std::endl;

  // print sendindices
  {
    for (int r = 0; r < commSize; ++r) {
      if (r == rank) {
        std::cout << "rank " << r << " send1dIndices ";
        for (const auto& s : send1dIndices) {
          std::cout << s.second;
        }
        std::cout << std::endl;
      }
      MPI_Barrier(comm);
    }
  }

  MPI_Barrier(comm);

#endif
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
    IndexVector sizes(dfg.getLocalSizes().begin(), dfg.getLocalSizes().end());

    // sizes of subarray ( full size except dimension d )
    IndexVector subsizes = sizes;
    subsizes[dim] = 1;

    // start
    IndexVector starts(dfg.getDimension(), 0);

    std::vector<int> csizes(sizes.begin(), sizes.end());
    std::vector<int> csubsizes(subsizes.begin(), subsizes.end());
    std::vector<int> cstarts(starts.begin(), starts.end());

    // create subarray view on data
    MPI_Type_create_subarray(static_cast<int>(dfg.getDimension()), &csizes[0], &csubsizes[0],
                             &cstarts[0], MPI_ORDER_FORTRAN, dfg.getMPIDatatype(), &mysubarray);
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
          IndexVector lidxvec(dfg.getDimension(), 0);
          IndexVector gidxvec = dfg.getLowerBounds();
          gidxvec[dim] = index;
          bool tmp = dfg.getLocalVectorIndex(gidxvec, lidxvec);
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

#ifdef DEBUG_OUTPUT
        auto rank = dfg.getRank();
        int mpiDSize = 0;
        MPI_Type_size(myHBlock, &mpiDSize);
        // MPI_Type_size(mysubarrayBlock, &mpiDSize);
        // print info: dest, size, index
        std::cout << "rank " << rank << ": send gindex " << *(indices.begin()) << " + " << displacements << " dest " << dest
                  // << " data " << *(dfg.getData()) << " " << *(dfg.getData()+1)
                  // << " full data " << dfg.getElementVector()
                  << " size " << mpiDSize
                  << std::endl;

#endif
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
        remoteData.emplace_back(sizes, dim, index, lowerBoundsNeighbor);

        FG_ELEMENT* buf = remoteData.back().getData();
        bufs.push_back(buf);
        assert(bsize == static_cast<int>(remoteData.back().getSize()));
      }
      {
        // make datatype hblock for all indices
        MPI_Datatype myHBlock;
        MPI_Aint firstBufAddr;
        MPI_Get_address(bufs[0], &firstBufAddr);
        std::vector<MPI_Aint> displacements;
        displacements.reserve(bufs.size());
        for (const auto& b : bufs) {
          MPI_Aint addr;
          MPI_Get_address(b, &addr);
          auto d = MPI_Aint_diff(addr, firstBufAddr);
          displacements.push_back(d);
        }
        assert(displacements[0] == 0);
        MPI_Type_create_hindexed_block(static_cast<int>(indices.size()), bsize, displacements.data(),
                                       dfg.getMPIDatatype(), &myHBlock);
        MPI_Type_commit(&myHBlock);
        // start recv operation, use first global index as tag
        {
          int src = static_cast<int>(r);
          int tag = static_cast<int>(*(indices.begin()));

          MPI_Irecv(static_cast<void*>(bufs[0]), 1, myHBlock, src, tag, dfg.getCommunicator(),
                    &recvRequests[numRecv++]);

#ifdef DEBUG_OUTPUT
          auto rank = dfg.getRank();
          int mpiDSize = 0;
          MPI_Type_size(myHBlock, &mpiDSize);
          // print info: dest, size, index
          std::cout << "rank " << rank << ": recv gindex " << *(indices.begin()) << " src " << src
                    << " size: " << bsize << "*" << indices.size()
                    << " size " << mpiDSize
                    << std::endl;
#endif
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

  sendAndReceiveIndices(send1dIndices, recv1dIndices, dfg, dim, remoteData);
  // sendAndReceiveIndicesBlock(send1dIndices, recv1dIndices, dfg, dim, remoteData);
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
                           std::vector<RemoteDataContainer<FG_ELEMENT> >& remoteData) {

#ifdef DEBUG_OUTPUT
  auto commSize = dfg.getCommunicatorSize();
  CommunicatorType comm = dfg.getCommunicator();
  auto rank = dfg.getRank();
  std::vector<int> coords(dfg.getDimension());
  dfg.getCartesianUtils().getPartitionCoordsOfLocalRank(coords);
  {
    std::cout << "in debug output" << std::endl;

    IndexType fidx = dfg.getFirstGlobal1dIndex(dim);
    LevelType flvl = dfg.getLevel(dim, fidx);
    IndexType fleftpre = dfg.getLeftPredecessor(dim, fidx);
    IndexType frightpre = dfg.getRightPredecessor(dim, fidx);
    RankType leftPreRank = dfg.getNeighbor1dFromAxisIndex(dim, fleftpre);
    RankType rightPreRank = dfg.getNeighbor1dFromAxisIndex(dim, frightpre);

    if (rank == 0) std::cout << "first point for dim: " << dim << std::endl;

    for (int r = 0; r < commSize; ++r) {
      if (r == rank) {
        std::cout << "rank " << rank << " "
                  << "coords " << coords << " "
                  << "idx " << fidx << " "
                  << "lvl " << flvl << " "
                  << "leftpre " << fleftpre << " "
                  << "right pre " << frightpre << " "
                  << "rank of left pre " << leftPreRank << " "
                  << "rank of righ pre " << rightPreRank << " " << std::endl;
      }

      MPI_Barrier(comm);
    }
  }

  MPI_Barrier(comm);

  {
    if (rank == 0) std::cout << "\n last point:" << std::endl;

    IndexType lidx = dfg.getLastGlobal1dIndex(dim);
    LevelType llvl = dfg.getLevel(dim, lidx);
    IndexType lleftpre = dfg.getLeftPredecessor(dim, lidx);
    IndexType lrightpre = dfg.getRightPredecessor(dim, lidx);
    RankType leftPreRank = dfg.getNeighbor1dFromAxisIndex(dim, lleftpre);
    RankType rightPreRank = dfg.getNeighbor1dFromAxisIndex(dim, lrightpre);

    for (int r = 0; r < commSize; ++r) {
      if (r == rank) {
        std::cout << "rank " << rank << " "
                  << "coords " << coords << " "
                  << "idx " << lidx << " "
                  << "lvl " << llvl << " "
                  << "leftpre " << lleftpre << " "
                  << "rightpre " << lrightpre << " "
                  << "rank of left pre " << leftPreRank << " "
                  << "rank of righ pre " << rightPreRank << " " << std::endl;
      }

      MPI_Barrier(comm);
    }
  }
#endif

  // create buffers for every rank
  std::map<RankType, std::set<IndexType>> recv1dIndices;
  std::map<RankType, std::set<IndexType>> send1dIndices;

  // main loop
  IndexType idxMin = dfg.getFirstGlobal1dIndex(dim);
  IndexType idxMax = dfg.getLastGlobal1dIndex(dim);
  auto globalIdxMax = dfg.length(dim);
  LevelType lmax = dfg.getLevels()[dim];

  IndexType idx = idxMin;

  // for hierarchization, we only need to exchange the direct predecessors
  while (idx <= idxMax) {
    LevelType lidx = dfg.getLevel(dim, idx);
    // check if successors of idx outside local domain
    {
      for (LevelType l = static_cast<LevelType>(lidx + 1); l <= lmax; ++l) {
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
        }
      }
    }
    auto ldiff = static_cast<LevelType>(lmax - lidx);
    IndexType idiff = static_cast<IndexType>(powerOfTwo[ldiff]);
    // check if predecessors of idx outside local domain
    IndexType pIdx;
    // left, right predecessor
    for (const auto& indexShift : {-1, 1}) {
      pIdx = idx + indexShift * idiff;
      // if we are not on the boundary level, and
      // pIdx is outside of my domain, but still in the global domain
      if (lidx > 0 && ((indexShift < 0 && pIdx >= 0 && pIdx < idxMin) ||
                        (indexShift > 0 && pIdx > idxMax && pIdx < globalIdxMax))) {
        // get rank which has predecessor and add to list of indices to recv
        int r = dfg.getNeighbor1dFromAxisIndex(dim, pIdx);
        recv1dIndices[r].insert(pIdx);
      }
    }
    if (lidx == 0 || pIdx > idxMax) {
      idx = getNextIndex1d(dfg, dim, idx);
    } else {
      // index of right predecessor
      idx = pIdx;
    }
  }

  sendAndReceiveIndices(send1dIndices, recv1dIndices, dfg, dim, remoteData);
  // sendAndReceiveIndicesBlock(send1dIndices, recv1dIndices, dfg, dim, remoteData);
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
    std::vector<RemoteDataContainer<FG_ELEMENT> >& remoteData) {

#ifdef DEBUG_OUTPUT
  auto commSize = dfg.getCommunicatorSize();
  CommunicatorType comm = dfg.getCommunicator();
  auto rank = dfg.getRank();
  std::vector<int> coords(dfg.getDimension());
  dfg.getCartesianUtils().getPartitionCoordsOfLocalRank(coords);
  {
    IndexType fidx = dfg.getFirstGlobal1dIndex(dim);
    LevelType flvl = dfg.getLevel(dim, fidx);
    IndexType fleftpre = dfg.getLeftPredecessor(dim, fidx);
    IndexType frightpre = dfg.getRightPredecessor(dim, fidx);
    RankType leftPreRank = dfg.getNeighbor1dFromAxisIndex(dim, fleftpre);
    RankType rightPreRank = dfg.getNeighbor1dFromAxisIndex(dim, frightpre);

    if (rank == 0) std::cout << "first point:" << std::endl;

    for (int r = 0; r < commSize; ++r) {
      if (r == rank) {
        std::cout << "rank " << rank << " "
                  << "coords " << coords << " "
                  << "idx " << fidx << " "
                  << "lvl " << flvl << " "
                  << "leftpre " << fleftpre << " "
                  << "right pre " << frightpre << " "
                  << "rank of left pre " << leftPreRank << " "
                  << "rank of righ pre " << rightPreRank << " " << std::endl;
      }

      MPI_Barrier(comm);
    }
  }

  MPI_Barrier(comm);

  {
    if (rank == 0) std::cout << "\n last point:" << std::endl;

    IndexType lidx = dfg.getLastGlobal1dIndex(dim);
    LevelType llvl = dfg.getLevel(dim, lidx);
    IndexType lleftpre = dfg.getLeftPredecessor(dim, lidx);
    IndexType lrightpre = dfg.getRightPredecessor(dim, lidx);
    RankType leftPreRank = dfg.getNeighbor1dFromAxisIndex(dim, lleftpre);
    RankType rightPreRank = dfg.getNeighbor1dFromAxisIndex(dim, lrightpre);

    for (int r = 0; r < commSize; ++r) {
      if (r == rank) {
        std::cout << "rank " << rank << " "
                  << "coords " << coords << " "
                  << "idx " << lidx << " "
                  << "lvl " << llvl << " "
                  << "leftpre " << lleftpre << " "
                  << "rightpre " << lrightpre << " "
                  << "rank of left pre " << leftPreRank << " "
                  << "rank of righ pre " << rightPreRank << " " << std::endl;
      }

      MPI_Barrier(comm);
    }
  }
#endif

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
    checkLeftSuccesors(idx, idx, dim, dfg, send1dIndices);

    checkRightSuccesors(idx, idx, dim, dfg, send1dIndices);

    idx = checkPredecessors(idx, dim, dfg, recv1dIndices);
  }

  sendAndReceiveIndices(send1dIndices, recv1dIndices, dfg, dim, remoteData);
  // sendAndReceiveIndicesBlock(send1dIndices, recv1dIndices, dfg, dim, remoteData);
}

template <typename FG_ELEMENT>
static void checkLeftSuccesors(IndexType checkIdx, IndexType rootIdx, DimType dim,
                               const DistributedFullGrid<FG_ELEMENT>& dfg,
                               std::map<RankType, std::set<IndexType>>& OneDIndices) {
  LevelType lidx = dfg.getLevel(dim, checkIdx);
  IndexType idxMin = dfg.getFirstGlobal1dIndex(dim);
  // IndexType idxMax = dfg.getLastGlobal1dIndex(dim);
  LevelType lmax = dfg.getLevels()[dim];

  // check left successors of checkIdx
  for (auto l = static_cast<LevelType>(lidx + 1); l <= lmax; ++l) {
    auto ldiff = static_cast<LevelType>(lmax - l);
    auto idiff = static_cast<IndexType>(powerOfTwo[ldiff]);

    IndexType lsIdx = checkIdx - idiff;

    if (lsIdx >= 0 && lsIdx < idxMin) {
      // get rank which has lsIdx and add to send list
      int r = dfg.getNeighbor1dFromAxisIndex(dim, lsIdx);
      if (r >= 0) OneDIndices[r].insert(rootIdx);
    }

    if (lsIdx >= 0) checkLeftSuccesors(lsIdx, rootIdx, dim, dfg, OneDIndices);
  }
}

template <typename FG_ELEMENT>
static void checkRightSuccesors(IndexType checkIdx, IndexType rootIdx, DimType dim,
                                const DistributedFullGrid<FG_ELEMENT>& dfg,
                                std::map<RankType, std::set<IndexType>>& OneDIndices) {
  LevelType lidx = dfg.getLevel(dim, checkIdx);

  IndexType idxMax = dfg.getLastGlobal1dIndex(dim);
  LevelType lmax = dfg.getLevels()[dim];

  // check right successors of checkIdx
  for (auto l = static_cast<LevelType>(lidx + 1); l <= lmax; ++l) {
    auto ldiff = static_cast<LevelType>(lmax - l);
    auto idiff = static_cast<IndexType>(powerOfTwo[ldiff]);

    IndexType rsIdx = checkIdx + idiff;

    if (rsIdx < dfg.getGlobalSizes()[dim] && rsIdx > idxMax) {
      // get rank which has rsIdx and add to send list
      int r = dfg.getNeighbor1dFromAxisIndex(dim, rsIdx);
      if (r >= 0) OneDIndices[r].insert(rootIdx);
    }

    if (rsIdx < dfg.length(dim)) {
      checkRightSuccesors(rsIdx, rootIdx, dim, dfg, OneDIndices);
    }
  }
}

template <typename FG_ELEMENT>
static IndexType checkPredecessors(IndexType idx, DimType dim,
                                   const DistributedFullGrid<FG_ELEMENT>& dfg,
                                   std::map<RankType, std::set<IndexType>>& OneDIndices) {
  IndexType idxMin = dfg.getFirstGlobal1dIndex(dim);
  IndexType idxMax = dfg.getLastGlobal1dIndex(dim);

  // check if left predecessor outside local domain
  // if returns negative value there's no left predecessor
  IndexType lpIdx = dfg.getLeftPredecessor(dim, idx);

  if (lpIdx >= 0 && lpIdx < idxMin) {
    // get rank which has left predecessor and add to list of indices
    int r = dfg.getNeighbor1dFromAxisIndex(dim, lpIdx);
    OneDIndices[r].insert(lpIdx);
  }
  if (lpIdx >= 0) checkPredecessors(lpIdx, dim, dfg, OneDIndices);

  // check if right predecessor outside local domain
  // if returns negative value there's no right predecessor
  IndexType rpIdx = dfg.getRightPredecessor(dim, idx);

  if (rpIdx < 0) {
    idx = getNextIndex1d(dfg, dim, idx);
    return idx;
  }

  if (rpIdx > idxMax) {
    // get rank which has right predecessor and add to list of indices to recv
    int r = dfg.getNeighbor1dFromAxisIndex(dim, rpIdx);
    OneDIndices[r].insert(rpIdx);
    idx = getNextIndex1d(dfg, dim, idx);
  } else {
    idx = rpIdx;
  }

  checkPredecessors(rpIdx, dim, dfg, OneDIndices);

  return idx;
}

/* returns the neighboring process (in the sense that the neighbor has the same
 * partion coordinates in all other dimensions than d) in dimension d which
 * contains the point with the one-dimensional index idx1d
 */
template <typename FG_ELEMENT>
RankType getNeighbor1d(const DistributedFullGrid<FG_ELEMENT>& dfg, DimType dim, IndexType idx1d) {
  // if global index is outside of domain return negative value
  {
    if (idx1d < 0) return -1;

    IndexType numElementsD = dfg.getGlobalSizes()[dim];

    if (idx1d > numElementsD - 1) return -1;
  }

  IndexVector globalAxisIndex = dfg.getLowerBounds();
  globalAxisIndex[dim] = idx1d;

  std::vector<int> partitionCoords(dfg.getDimension());
  dfg.getPartitionCoords(globalAxisIndex, partitionCoords);
  RankType r = dfg.getCartesianUtils().getRankFromPartitionCoords(partitionCoords);

  // check if global index vector is actually contained in the domain of rank r
  assert(globalAxisIndex >= dfg.getLowerBounds(r));
  assert(globalAxisIndex < dfg.getUpperBounds(r));
  assert(r < dfg.getCommunicatorSize());
  return r;
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

template <typename FG_ELEMENT>
static IndexType getLastIndexOfLevel1d(const DistributedFullGrid<FG_ELEMENT>& dfg, DimType dim,
                                        LevelType l) {
  IndexType idxMax = dfg.getLastGlobal1dIndex(dim);
  IndexType idxMin = dfg.getFirstGlobal1dIndex(dim);

  for (IndexType i = idxMax; i >= idxMin; --i) {
    if (dfg.getLevel(dim, i) == l) return i;
  }

  // no index on level l found
  return -1;
}


template <typename FG_ELEMENT>
static void hierarchizeX_opt_noboundary(DistributedFullGrid<FG_ELEMENT>& dfg,
                                        LookupTable<FG_ELEMENT>& lookupTable) {
  const DimType dim = 0;
  assert(dfg.returnBoundaryFlags()[dim] == false);

  LevelType lmax = dfg.getLevels()[dim];
  IndexType idxMax = dfg.getLastGlobal1dIndex(dim);
  IndexType ndim = dfg.getLocalSizes()[dim];

  // size of xBlcok
  IndexType xSize = dfg.getLocalSizes()[0];

  // FG_ELEMENT zeroVal(0);

  // create tmp array to store xblock
  std::vector<FG_ELEMENT> tmp(dfg.getGlobalSizes()[dim]);
  std::vector<FG_ELEMENT>& localData = dfg.getElementVector();

  // loop over all xBlocks of local domain -> linearIndex with stride localndim[0]
  IndexType nbrxBlocks = dfg.getNrLocalElements() / ndim;

  for (IndexType xBlock = 0; xBlock < nbrxBlocks; ++xBlock) {
    // get globalIndexVector of block start
    // this is the base IndexVector of this block
    // only dim component is varied
    IndexType linIdxBlockStart = xBlock * ndim;

    IndexVector localIndexVector(dfg.getDimension());
    IndexVector baseGlobalIndexVector(dfg.getDimension());
    dfg.getLocalVectorIndex(linIdxBlockStart, localIndexVector);
    dfg.getGlobalVectorIndex(localIndexVector, baseGlobalIndexVector);
    assert(localIndexVector[dim] == 0);

    // copy local data to tmp
    for (IndexType i = 0; i < xSize; ++i)
      tmp[baseGlobalIndexVector[dim] + i] = localData[linIdxBlockStart + i];

    // copy remote data to tmp
    std::vector<RemoteDataContainer<FG_ELEMENT> >& rdcs = lookupTable.getRDCVector();
    IndexVector tmpGlobalIndexVector = baseGlobalIndexVector;

    // go through remote containers
    for (size_t i = 0; i < rdcs.size(); ++i) {
      IndexType global1didx = rdcs[i].getKeyIndex();
      tmpGlobalIndexVector[dim] = global1didx;
      tmp[global1didx] = *rdcs[i].getData(tmpGlobalIndexVector);
    }

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

    // copy local data back
    for (IndexType i = 0; i < xSize; ++i)
      localData[linIdxBlockStart + i] = tmp[baseGlobalIndexVector[dim] + i];
  }
}

template <typename FG_ELEMENT>
static void dehierarchizeX_opt_noboundary(DistributedFullGrid<FG_ELEMENT>& dfg,
                                          LookupTable<FG_ELEMENT>& lookupTable) {
  const DimType dim = 0;
  assert(dfg.returnBoundaryFlags()[dim] == false);

  LevelType lmax = dfg.getLevels()[dim];
  IndexType ndim = dfg.getLocalSizes()[dim];

  // size of xBlcok
  IndexType xSize = dfg.getLocalSizes()[0];

  // create tmp array to store xblock
  std::vector<FG_ELEMENT> tmp(dfg.getGlobalSizes()[dim]);
  std::vector<FG_ELEMENT>& localData = dfg.getElementVector();

  // loop over all xBlocks of local domain -> linearIndex with stride localndim[0]
  IndexType nbrxBlocks = dfg.getNrLocalElements() / ndim;

  for (IndexType xBlock = 0; xBlock < nbrxBlocks; ++xBlock) {
    // get globalIndexVector of block start
    // this is the base IndexVector of this block
    // only dim component is varied
    IndexType linIdxBlockStart = xBlock * ndim;

    IndexVector localIndexVector(dfg.getDimension());
    IndexVector baseGlobalIndexVector(dfg.getDimension());
    dfg.getLocalVectorIndex(linIdxBlockStart, localIndexVector);
    dfg.getGlobalVectorIndex(localIndexVector, baseGlobalIndexVector);
    assert(localIndexVector[dim] == 0);

    // copy local data to tmp
    for (IndexType i = 0; i < xSize; ++i)
      tmp[baseGlobalIndexVector[dim] + i] = localData[linIdxBlockStart + i];

    // copy remote data to tmp
    std::vector<RemoteDataContainer<FG_ELEMENT> >& rdcs = lookupTable.getRDCVector();
    IndexVector tmpGlobalIndexVector = baseGlobalIndexVector;

    // go through remote containers
    for (size_t i = 0; i < rdcs.size(); ++i) {
      IndexType global1didx = rdcs[i].getKeyIndex();
      tmpGlobalIndexVector[dim] = global1didx;
      tmp[global1didx] = *rdcs[i].getData(tmpGlobalIndexVector);
    }

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

    // copy local data back
    for (IndexType i = 0; i < xSize; ++i)
      localData[linIdxBlockStart + i] = tmp[baseGlobalIndexVector[dim] + i];
  }
}

template <typename FG_ELEMENT, bool periodic = false>
inline void hierarchizeX_opt_boundary_kernel(FG_ELEMENT* data, LevelType lmax, int start,
                                             int stride, LevelType lmin) {
  assert(lmin == 0);
  const int lmaxi = static_cast<int>(lmax);
  int ll = lmaxi;
  int steps = (1 << (lmaxi - 1));
  int offset = 1;  // 1 and not 0 because boundary
  int stepsize = 2;
  int parentOffset = 1;

  if (periodic) {
    int idxmax = powerOfTwo[lmaxi];
    data[0] = 0.5 * (data[idxmax] + data[0]);
    data[idxmax] = data[0];
  }

  for (ll--; ll > -1; ll--) {
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
inline void hierarchizeX_full_weighting_boundary_kernel(FG_ELEMENT* data, LevelType lmax,
                                                            int start, int stride,
                                                            LevelType lmin) {
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
inline void hierarchizeX_biorthogonal_boundary_kernel(FG_ELEMENT* data, LevelType lmax, int start,
                                             int stride, LevelType lmin) {

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
inline void dehierarchizeX_opt_boundary_kernel(FG_ELEMENT* data, LevelType lmax, int start,
                                               int stride, LevelType lmin) {
  assert(lmin == 0);
  const int lmaxi = static_cast<int>(lmax);
  int steps = 1;
  int offset = (1 << (lmaxi - 1));  // offset =1 da boundary.
  int stepsize = (1 << lmaxi);
  int parentOffset = (1 << (lmaxi - 1));

  for (LevelType ll = 1; ll <= lmax; ++ll) {
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
 * @brief inverse operation to hierarchizeX_full_weighting_boundary_kernel
 */
template <typename FG_ELEMENT, bool periodic = false>
inline void dehierarchizeX_full_weighting_boundary_kernel(FG_ELEMENT* data, LevelType lmax,
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
 * @brief inverse operation to hierarchizeX_biorthogonal_boundary_kernel
 */
template <typename FG_ELEMENT, bool periodic = false>
inline void dehierarchizeX_biorthogonal_boundary_kernel(FG_ELEMENT* data, LevelType lmax, int start,
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
 * @brief hierarchize a DFG in dimension X (with contiguous access)
 *
 * @param dfg : the DFG to hierarchize
 * @param lookupTable: the lookup table for local and remote data
 */
template <typename FG_ELEMENT,
          void (*HIERARCHIZATION_FCTN)(FG_ELEMENT[], LevelType, int, int, LevelType)>
static void hierarchizeX_opt_boundary(DistributedFullGrid<FG_ELEMENT>& dfg,
                                      LookupTable<FG_ELEMENT>& lookupTable) {
  const DimType dim = 0;
  assert(dfg.returnBoundaryFlags()[dim] == true);

  LevelType lmax = dfg.getLevels()[dim];
  IndexType ndim = dfg.getLocalSizes()[dim];

  // size of xBlcok
  IndexType xSize = ndim;

  // create tmp array to store xblock
  std::vector<FG_ELEMENT> tmp(dfg.getGlobalSizes()[dim]);
  std::vector<FG_ELEMENT>& localData = dfg.getElementVector();

  IndexVector localIndexVector(dfg.getDimension());
  IndexVector baseGlobalIndexVector(dfg.getDimension());
  IndexVector tmpGlobalIndexVector(dfg.getDimension());

  IndexType gstart = dfg.getLowerBounds()[dim];
  IndexType linIdxBlockStart;

  // loop over all xBlocks of local domain -> linearIndex with stride localndim[0]
  IndexType nbrxBlocks = dfg.getNrLocalElements() / ndim;

  for (IndexType xBlock = 0; xBlock < nbrxBlocks; ++xBlock) {
    // get globalIndexVector of block start
    // this is the base IndexVector of this block
    // only dim component is varied
    linIdxBlockStart = xBlock * ndim;

    dfg.getLocalVectorIndex(linIdxBlockStart, localIndexVector);
    dfg.getGlobalVectorIndex(localIndexVector, tmpGlobalIndexVector);
    assert(localIndexVector[dim] == 0);

    // copy local data to tmp
    for (IndexType i = 0; i < xSize; ++i) tmp[gstart + i] = localData[linIdxBlockStart + i];

    // copy remote data to tmp
    std::vector<RemoteDataContainer<FG_ELEMENT> >& rdcs = lookupTable.getRDCVector();

    // go through remote containers
    for (size_t i = 0; i < rdcs.size(); ++i) {
      IndexType global1didx = rdcs[i].getKeyIndex();
      tmpGlobalIndexVector[dim] = global1didx;
      tmp[global1didx] = *rdcs[i].getData(tmpGlobalIndexVector);
    }

    HIERARCHIZATION_FCTN(&tmp[0], lmax, 0, 1, 0);

    // copy local data back
    for (IndexType i = 0; i < xSize; ++i) localData[linIdxBlockStart + i] = tmp[gstart + i];
  }
}

/**
 * @brief inverse operation for hierarchizeX_opt_boundary
 */
template <typename FG_ELEMENT,
          void (*DEHIERARCHIZATION_FCTN)(FG_ELEMENT[], LevelType, int, int,
                                       LevelType) = dehierarchizeX_opt_boundary_kernel>
static void dehierarchizeX_opt_boundary(DistributedFullGrid<FG_ELEMENT>& dfg,
                                        LookupTable<FG_ELEMENT>& lookupTable) {
  const DimType dim = 0;
  assert(dfg.returnBoundaryFlags()[dim] == true);

  LevelType lmax = dfg.getLevels()[dim];
  IndexType ndim = dfg.getLocalSizes()[dim];

  // size of xBlcok
  IndexType xSize = ndim;

  // create tmp array to store xblock
  std::vector<FG_ELEMENT> tmp(dfg.getGlobalSizes()[dim]);
  std::vector<FG_ELEMENT>& localData = dfg.getElementVector();

  IndexVector localIndexVector(dfg.getDimension());
  IndexVector baseGlobalIndexVector(dfg.getDimension());
  IndexVector tmpGlobalIndexVector(dfg.getDimension());

  IndexType gstart = dfg.getLowerBounds()[dim];
  IndexType linIdxBlockStart;

  // loop over all xBlocks of local domain -> linearIndex with stride localndim[0]
  IndexType nbrxBlocks = dfg.getNrLocalElements() / ndim;

  for (IndexType xBlock = 0; xBlock < nbrxBlocks; ++xBlock) {
    // get globalIndexVector of block start
    // this is the base IndexVector of this block
    // only dim component is varied
    linIdxBlockStart = xBlock * ndim;

    dfg.getLocalVectorIndex(linIdxBlockStart, localIndexVector);
    dfg.getGlobalVectorIndex(localIndexVector, tmpGlobalIndexVector);
    assert(localIndexVector[dim] == 0);

    // copy local data to tmp
    for (IndexType i = 0; i < xSize; ++i) tmp[gstart + i] = localData[linIdxBlockStart + i];

    // copy remote data to tmp
    std::vector<RemoteDataContainer<FG_ELEMENT> >& rdcs = lookupTable.getRDCVector();

    // go through remote containers
    for (size_t i = 0; i < rdcs.size(); ++i) {
      IndexType global1didx = rdcs[i].getKeyIndex();
      tmpGlobalIndexVector[dim] = global1didx;
      tmp[global1didx] = *rdcs[i].getData(tmpGlobalIndexVector);
    }

    DEHIERARCHIZATION_FCTN(&tmp[0], lmax, 0, 1, 0);
    // copy local data back
    for (IndexType i = 0; i < xSize; ++i) localData[linIdxBlockStart + i] = tmp[gstart + i];
  }
}

/**
 * @brief  hierarchize a DFG with boundary points in dimension dim (with non-contiguous access)
 */
template <typename FG_ELEMENT,
          void (*HIERARCHIZATION_FCTN)(FG_ELEMENT[], LevelType, int, int,
                                       LevelType) = hierarchizeX_opt_boundary_kernel>
void hierarchizeN_opt_boundary(DistributedFullGrid<FG_ELEMENT>& dfg,
                               LookupTable<FG_ELEMENT>& lookupTable, DimType dim) {
  assert(dfg.returnBoundaryFlags()[dim] == true);

  LevelType lmax = dfg.getLevels()[dim];
  IndexType size = dfg.getNrLocalElements();
  IndexType stride = dfg.getLocalOffsets()[dim];
  IndexType ndim = dfg.getLocalSizes()[dim];
  IndexType jump = stride * ndim;
  IndexType nbrOfPoles = size / ndim;

  IndexVector localIndexVector(dfg.getDimension());
  IndexVector baseGlobalIndexVector(dfg.getDimension());
  IndexVector tmpGlobalIndexVector(dfg.getDimension());

  // loop over poles
  std::vector<FG_ELEMENT> tmp(dfg.getGlobalSizes()[dim]);
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
    dfg.getGlobalVectorIndex(localIndexVector, tmpGlobalIndexVector);
    assert(localIndexVector[dim] == 0);

    // copy remote data to tmp
    std::vector<RemoteDataContainer<FG_ELEMENT> >& rdcs = lookupTable.getRDCVector();

    // go through remote containers
    for (size_t i = 0; i < rdcs.size(); ++i) {
      IndexType global1didx = rdcs[i].getKeyIndex();
      tmpGlobalIndexVector[dim] = global1didx;
      tmp[global1didx] = *rdcs[i].getData(tmpGlobalIndexVector);
    }

    // copy local data
    for (IndexType i = 0; i < ndim; ++i) tmp[gstart + i] = ldata[start + stride * i];

    // hierarchize tmp array with hupp function
    HIERARCHIZATION_FCTN(&tmp[0], lmax, 0, 1, 0);

    // copy pole back
    for (IndexType i = 0; i < ndim; ++i) ldata[start + stride * i] = tmp[gstart + i];
  }
}

/**
 * @brief  hierarchize a DFG without boundary points in dimension dim (with non-contiguous access)
 */
template <typename FG_ELEMENT>
void hierarchizeN_opt_noboundary(DistributedFullGrid<FG_ELEMENT>& dfg,
                                 LookupTable<FG_ELEMENT>& lookupTable, DimType dim) {
  assert(dfg.returnBoundaryFlags()[dim] == false);

  LevelType lmax = dfg.getLevels()[dim];
  IndexType size = dfg.getNrLocalElements();
  IndexType stride = dfg.getLocalOffsets()[dim];
  IndexType ndim = dfg.getLocalSizes()[dim];
  IndexType jump = stride * ndim;
  IndexType nbrOfPoles = size / ndim;

  IndexVector localIndexVector(dfg.getDimension());
  IndexVector baseGlobalIndexVector(dfg.getDimension());
  IndexVector tmpGlobalIndexVector(dfg.getDimension());

  // loop over poles
  std::vector<FG_ELEMENT> tmp(dfg.getGlobalSizes()[dim]);
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
    dfg.getGlobalVectorIndex(localIndexVector, tmpGlobalIndexVector);
    assert(localIndexVector[dim] == 0);

    // copy remote data to tmp
    std::vector<RemoteDataContainer<FG_ELEMENT> >& rdcs = lookupTable.getRDCVector();

    // go through remote containers
    for (size_t i = 0; i < rdcs.size(); ++i) {
      IndexType global1didx = rdcs[i].getKeyIndex();
      tmpGlobalIndexVector[dim] = global1didx;
      tmp[global1didx] = *rdcs[i].getData(tmpGlobalIndexVector);
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
 * @brief inverse operation for hierarchizeN_opt_boundary
 */
template <typename FG_ELEMENT,
          void (*DEHIERARCHIZATION_FCTN)(FG_ELEMENT[], LevelType, int, int,
                                       LevelType) = dehierarchizeX_opt_boundary_kernel>
void dehierarchizeN_opt_boundary(DistributedFullGrid<FG_ELEMENT>& dfg,
                                 LookupTable<FG_ELEMENT>& lookupTable, DimType dim) {
  assert(dfg.returnBoundaryFlags()[dim] == true);

  const auto& lmax = dfg.getLevels()[dim];
  const auto& size = dfg.getNrLocalElements();
  const auto& stride = dfg.getLocalOffsets()[dim];
  const auto& ndim = dfg.getLocalSizes()[dim];
  IndexType jump = stride * ndim;
  IndexType nbrOfPoles = size / ndim;

  IndexVector localIndexVector(dfg.getDimension());
  IndexVector baseGlobalIndexVector(dfg.getDimension());
  IndexVector tmpGlobalIndexVector(dfg.getDimension());

  // loop over poles
  std::vector<FG_ELEMENT> tmp(dfg.getGlobalSizes()[dim]);
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
    dfg.getGlobalVectorIndex(localIndexVector, tmpGlobalIndexVector);
    assert(localIndexVector[dim] == 0);

    // copy remote data to tmp
    std::vector<RemoteDataContainer<FG_ELEMENT> >& rdcs = lookupTable.getRDCVector();

    // go through remote containers
    for (size_t i = 0; i < rdcs.size(); ++i) {
      IndexType global1didx = rdcs[i].getKeyIndex();
      tmpGlobalIndexVector[dim] = global1didx;
      tmp[global1didx] = *rdcs[i].getData(tmpGlobalIndexVector);
    }

    // copy local data
    for (IndexType i = 0; i < ndim; ++i) tmp[gstart + i] = ldata[start + stride * i];

    // hierarchize tmp array with hupp function
    DEHIERARCHIZATION_FCTN(&tmp[0], lmax, 0, 1, 0);

    // copy pole back
    for (IndexType i = 0; i < ndim; ++i) ldata[start + stride * i] = tmp[gstart + i];
  }
}

/**
 * @brief inverse operation for hierarchizeN_opt_noboundary
 */
template <typename FG_ELEMENT>
void dehierarchizeN_opt_noboundary(DistributedFullGrid<FG_ELEMENT>& dfg,
                                   LookupTable<FG_ELEMENT>& lookupTable, DimType dim) {
  assert(dfg.returnBoundaryFlags()[dim] == false);

  LevelType lmax = dfg.getLevels()[dim];
  IndexType size = dfg.getNrLocalElements();
  IndexType stride = dfg.getLocalOffsets()[dim];
  IndexType ndim = dfg.getLocalSizes()[dim];
  IndexType jump = stride * ndim;
  IndexType nbrOfPoles = size / ndim;

  IndexVector localIndexVector(dfg.getDimension());
  IndexVector baseGlobalIndexVector(dfg.getDimension());
  IndexVector tmpGlobalIndexVector(dfg.getDimension());

  // loop over poles
  std::vector<FG_ELEMENT> tmp(dfg.getGlobalSizes()[dim]);
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
    dfg.getGlobalVectorIndex(localIndexVector, tmpGlobalIndexVector);
    assert(localIndexVector[dim] == 0);

    // copy remote data to tmp
    std::vector<RemoteDataContainer<FG_ELEMENT> >& rdcs = lookupTable.getRDCVector();

    // go through remote containers
    for (size_t i = 0; i < rdcs.size(); ++i) {
      IndexType global1didx = rdcs[i].getKeyIndex();
      tmpGlobalIndexVector[dim] = global1didx;
      tmp[global1didx] = *rdcs[i].getData(tmpGlobalIndexVector);
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

}  // unnamed namespace

namespace combigrid {

class DistributedHierarchization {
 public:
  // inplace hierarchization
  template <typename FG_ELEMENT>
  static void hierarchize(DistributedFullGrid<FG_ELEMENT>& dfg, const std::vector<bool>& dims,
                          const std::vector<BasisFunctionBasis*>& hierarchicalBases) {
    assert(dfg.getDimension() > 0);
    assert(dfg.getDimension() == dims.size());
    // hierarchize all dimensions, with special treatment for 0
    for (DimType dim = 0; dim < dfg.getDimension(); ++dim) {
      if (!dims[dim]) continue;

      // exchange data
      std::vector<RemoteDataContainer<FG_ELEMENT>> remoteData;
      if (dynamic_cast<HierarchicalHatBasisFunction*>(hierarchicalBases[dim]) != nullptr) {
        exchangeData1d(dfg, dim, remoteData);
      } else {
        exchangeAllData1d(dfg, dim, remoteData);
      }
      LookupTable<FG_ELEMENT> lookupTable(remoteData, dfg, dim);

      if (dfg.returnBoundaryFlags()[dim] == true) {
        // sorry for the code duplication, could not figure out a clean way
        if (dynamic_cast<HierarchicalHatBasisFunction*>(hierarchicalBases[dim]) != nullptr) {
          if (dim == 0)
            hierarchizeX_opt_boundary<FG_ELEMENT, hierarchizeX_opt_boundary_kernel<FG_ELEMENT>>(
                dfg, lookupTable);
          else
            hierarchizeN_opt_boundary<FG_ELEMENT, hierarchizeX_opt_boundary_kernel<FG_ELEMENT>>(
                dfg, lookupTable, dim);
        } else if (dynamic_cast<FullWeightingBasisFunction*>(hierarchicalBases[dim]) != nullptr) {
          if (dim == 0)
            hierarchizeX_opt_boundary<
                FG_ELEMENT, hierarchizeX_full_weighting_boundary_kernel<FG_ELEMENT, false>>(
                dfg, lookupTable);
          else
            hierarchizeN_opt_boundary<
                FG_ELEMENT, hierarchizeX_full_weighting_boundary_kernel<FG_ELEMENT, false>>(
                dfg, lookupTable, dim);
        } else if (dynamic_cast<FullWeightingPeriodicBasisFunction*>(hierarchicalBases[dim]) !=
                   nullptr) {
          if (dim == 0)
            hierarchizeX_opt_boundary<
                FG_ELEMENT, hierarchizeX_full_weighting_boundary_kernel<FG_ELEMENT, true>>(
                dfg, lookupTable);
          else
            hierarchizeN_opt_boundary<
                FG_ELEMENT, hierarchizeX_full_weighting_boundary_kernel<FG_ELEMENT, true>>(
                dfg, lookupTable, dim);
        } else if (dynamic_cast<BiorthogonalBasisFunction*>(hierarchicalBases[dim]) != nullptr) {
          if (dim == 0)
            hierarchizeX_opt_boundary<
                FG_ELEMENT, hierarchizeX_full_weighting_boundary_kernel<FG_ELEMENT, false>>(
                dfg, lookupTable);
          else
            hierarchizeN_opt_boundary<
                FG_ELEMENT, hierarchizeX_full_weighting_boundary_kernel<FG_ELEMENT, false>>(
                dfg, lookupTable, dim);
        } else if (dynamic_cast<BiorthogonalPeriodicBasisFunction*>(hierarchicalBases[dim]) !=
                   nullptr) {
          if (dim == 0)
            hierarchizeX_opt_boundary<
                FG_ELEMENT, hierarchizeX_biorthogonal_boundary_kernel<FG_ELEMENT, true>>(
                dfg, lookupTable);
          else
            hierarchizeN_opt_boundary<
                FG_ELEMENT, hierarchizeX_biorthogonal_boundary_kernel<FG_ELEMENT, true>>(
                dfg, lookupTable, dim);
        } else {
          throw std::logic_error("Not implemented");
        }
      } else {
        if (dynamic_cast<HierarchicalHatBasisFunction*>(hierarchicalBases[dim]) == nullptr) {
          throw std::logic_error("currently only hats supported for non-boundary grids");
        }
        if (dim == 0) {
          hierarchizeX_opt_noboundary(dfg, lookupTable);
        } else {
          hierarchizeN_opt_noboundary(dfg, lookupTable, dim);
        }
      }
    }
  }

  template <typename FG_ELEMENT, class T = HierarchicalHatBasisFunction>
  static void hierarchizeHierachicalBasis(DistributedFullGrid<FG_ELEMENT>& dfg, const std::vector<bool>& dims) {
    T basisFctn;
    std::vector<BasisFunctionBasis*> bases(dfg.getDimension(), &basisFctn);
    return hierarchize<FG_ELEMENT>(dfg, dims, bases);
  }

  template <typename FG_ELEMENT>
  static void hierarchize(DistributedFullGrid<FG_ELEMENT>& dfg) {
    std::vector<bool> dims(dfg.getDimension(), true);
    return hierarchizeHierachicalBasis<FG_ELEMENT,HierarchicalHatBasisFunction>(dfg, dims);
  }

  // inplace dehierarchization
  template <typename FG_ELEMENT>
  static void dehierarchize(DistributedFullGrid<FG_ELEMENT>& dfg, const std::vector<bool>& dims,
                              const std::vector<BasisFunctionBasis*>& hierarchicalBases) {
    assert(dfg.getDimension() > 0);
    assert(dfg.getDimension() == dims.size());
    // dehierarchize all dimensions, with special treatment for 0
    for (DimType dim = 0; dim < dfg.getDimension(); ++dim) {
      if (!dims[dim]) continue;

      // exchange data
      std::vector<RemoteDataContainer<FG_ELEMENT>> remoteData;
      if (dynamic_cast<HierarchicalHatBasisFunction*>(hierarchicalBases[dim]) != nullptr) {
        exchangeData1dDehierarchization(dfg, dim, remoteData);
      } else {
        exchangeAllData1d(dfg, dim, remoteData);
      }
      LookupTable<FG_ELEMENT> lookupTable(remoteData, dfg, dim);

      if (dfg.returnBoundaryFlags()[dim] == true) {
        // sorry for the code duplication, could not figure out a clean way
        if (dynamic_cast<HierarchicalHatBasisFunction*>(hierarchicalBases[dim]) != nullptr) {
          if (dim == 0)
            dehierarchizeX_opt_boundary<FG_ELEMENT, dehierarchizeX_opt_boundary_kernel<FG_ELEMENT>>(
                dfg, lookupTable);
          else
            dehierarchizeN_opt_boundary<FG_ELEMENT, dehierarchizeX_opt_boundary_kernel<FG_ELEMENT>>(
                dfg, lookupTable, dim);
        } else if (dynamic_cast<FullWeightingBasisFunction*>(hierarchicalBases[dim]) != nullptr) {
          if (dim == 0)
            dehierarchizeX_opt_boundary<
                FG_ELEMENT, dehierarchizeX_full_weighting_boundary_kernel<FG_ELEMENT, false>>(
                dfg, lookupTable);
          else
            dehierarchizeN_opt_boundary<
                FG_ELEMENT, dehierarchizeX_full_weighting_boundary_kernel<FG_ELEMENT, false>>(
                dfg, lookupTable, dim);
        } else if (dynamic_cast<FullWeightingPeriodicBasisFunction*>(hierarchicalBases[dim]) !=
                   nullptr) {
          if (dim == 0)
            dehierarchizeX_opt_boundary<
                FG_ELEMENT, dehierarchizeX_full_weighting_boundary_kernel<FG_ELEMENT, true>>(
                dfg, lookupTable);
          else
            dehierarchizeN_opt_boundary<
                FG_ELEMENT, dehierarchizeX_full_weighting_boundary_kernel<FG_ELEMENT, true>>(
                dfg, lookupTable, dim);
        } else if (dynamic_cast<BiorthogonalBasisFunction*>(hierarchicalBases[dim]) != nullptr) {
          if (dim == 0)
            dehierarchizeX_opt_boundary<
                FG_ELEMENT, dehierarchizeX_full_weighting_boundary_kernel<FG_ELEMENT, false>>(
                dfg, lookupTable);
          else
            dehierarchizeN_opt_boundary<
                FG_ELEMENT, dehierarchizeX_full_weighting_boundary_kernel<FG_ELEMENT, false>>(
                dfg, lookupTable, dim);
        } else if (dynamic_cast<BiorthogonalPeriodicBasisFunction*>(hierarchicalBases[dim]) !=
                   nullptr) {
          if (dim == 0)
            dehierarchizeX_opt_boundary<
                FG_ELEMENT, dehierarchizeX_biorthogonal_boundary_kernel<FG_ELEMENT, true>>(
                dfg, lookupTable);
          else
            dehierarchizeN_opt_boundary<
                FG_ELEMENT, dehierarchizeX_biorthogonal_boundary_kernel<FG_ELEMENT, true>>(
                dfg, lookupTable, dim);
        } else {
          throw std::logic_error("Not implemented");
        }
      } else {
        if (dynamic_cast<HierarchicalHatBasisFunction*>(hierarchicalBases[dim]) == nullptr) {
          throw std::logic_error("currently only hats supported for non-boundary grids");
        }
        if (dim == 0) {
          dehierarchizeX_opt_noboundary(dfg, lookupTable);
        } else {
          dehierarchizeN_opt_noboundary(dfg, lookupTable, dim);
        }
      }
    }
  }

  template <typename FG_ELEMENT, class T = HierarchicalHatBasisFunction>
  static void dehierarchizeHierachicalBasis(DistributedFullGrid<FG_ELEMENT>& dfg, const std::vector<bool>& dims) {
    T basisFctn;
    std::vector<BasisFunctionBasis*> bases(dfg.getDimension(), &basisFctn);
    return dehierarchize<FG_ELEMENT>(dfg, dims, bases);
  }

  template <typename FG_ELEMENT>
  static void dehierarchize(DistributedFullGrid<FG_ELEMENT>& dfg) {
    std::vector<bool> dims(dfg.getDimension(), true);
    return dehierarchizeHierachicalBasis<FG_ELEMENT,HierarchicalHatBasisFunction>(dfg, dims);
  }

  template <typename FG_ELEMENT>
  static void dehierarchizeDFG(DistributedFullGrid<FG_ELEMENT>& dfg,
                              const std::vector<bool>& hierarchizationDims,
                              const std::vector<BasisFunctionBasis*>& hierarchicalBases) {
    // dehierarchize dfg
    DistributedHierarchization::dehierarchize<FG_ELEMENT>(dfg, hierarchizationDims,
                                                          hierarchicalBases);
  }

  template<typename FG_ELEMENT>
  using FunctionPointer = void(*)(DistributedFullGrid<FG_ELEMENT>& dfg,
                                                    const std::vector<bool>& dims);
  // make template specifications visible by alias, hat is the default
  template <typename FG_ELEMENT>
  constexpr static FunctionPointer<FG_ELEMENT> hierarchizeHierarchicalHat =
      &hierarchizeHierachicalBasis<FG_ELEMENT, HierarchicalHatBasisFunction>;

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
