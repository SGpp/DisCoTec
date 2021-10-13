#ifndef DISTRIBUTEDHIERARCHIZATION_HPP_
#define DISTRIBUTEDHIERARCHIZATION_HPP_

//#define DEBUG_OUTPUT

#include "boost/lexical_cast.hpp"
#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGrid.hpp"
#include "sgpp/distributedcombigrid/legacy/combigrid_utils.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/DistributedSparseGridUniform.hpp"
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

    dim_ = sizes.size();
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
      IndexVector localIndexVector(dfg_.getDimension());
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
static void exchangeData1d(DistributedFullGrid<FG_ELEMENT>& dfg, DimType dim,
                           std::vector<RemoteDataContainer<FG_ELEMENT> >& remoteData);

template <typename FG_ELEMENT>
static void exchangeData1dDehierarchization(
    DistributedFullGrid<FG_ELEMENT>& dfg, DimType dim,
    std::vector<RemoteDataContainer<FG_ELEMENT> >& remoteData);

template <typename FG_ELEMENT>
static void checkLeftSuccesors(IndexType checkIdx, IndexType rootIdx, DimType dim,
                               DistributedFullGrid<FG_ELEMENT>& dfg,
                               std::vector<std::set<IndexType>>& OneDIndices);

template <typename FG_ELEMENT>
static void checkRightSuccesors(IndexType checkIdx, IndexType rootIdx, DimType dim,
                                DistributedFullGrid<FG_ELEMENT>& dfg,
                                std::vector<std::set<IndexType>>& OneDIndices);

template <typename FG_ELEMENT>
static IndexType checkPredecessors(IndexType idx, DimType dim, DistributedFullGrid<FG_ELEMENT>& dfg,
                                   std::vector<std::set<IndexType>>& OneDIndices);

template <typename FG_ELEMENT>
static RankType getNeighbor1d(DistributedFullGrid<FG_ELEMENT>& dfg, DimType d, IndexType idx1d);

template <typename FG_ELEMENT>
static IndexType getNextIndex1d(DistributedFullGrid<FG_ELEMENT>& dfg, DimType d, IndexType idx1d);

template <typename FG_ELEMENT>
static IndexType getFirstIndexOfLevel1d(DistributedFullGrid<FG_ELEMENT>& dfg, DimType d,
                                        LevelType l);
template <typename FG_ELEMENT>
static IndexType getLastIndexOfLevel1d(DistributedFullGrid<FG_ELEMENT>& dfg, DimType d,
                                        LevelType l);

template <typename FG_ELEMENT>
static void hierarchizeX(DistributedFullGrid<FG_ELEMENT>& dfg,
                         LookupTable<FG_ELEMENT>& lookupTable);

template <typename FG_ELEMENT>
static void hierarchizeX_opt_boundary(DistributedFullGrid<FG_ELEMENT>& dfg,
                                      LookupTable<FG_ELEMENT>& lookupTable);

template <typename FG_ELEMENT>
static void dehierarchizeX_opt_boundary(DistributedFullGrid<FG_ELEMENT>& dfg,
                                        LookupTable<FG_ELEMENT>& lookupTable);

template <typename FG_ELEMENT>
static void hierarchizeX_opt_noboundary(DistributedFullGrid<FG_ELEMENT>& dfg,
                                        LookupTable<FG_ELEMENT>& lookupTable);

template <typename FG_ELEMENT>
inline void hierarchizeX_opt_boundary_kernel(FG_ELEMENT* data, LevelType lmax, int start,
                                             int stride);

template <typename FG_ELEMENT>
inline void dehierarchizeX_opt_boundary_kernel(FG_ELEMENT* data, LevelType lmax, int start,
                                               int stride);


template <typename FG_ELEMENT>
void dehierarchizeN_opt_boundary(DistributedFullGrid<FG_ELEMENT>& dfg,
                                 LookupTable<FG_ELEMENT>& lookupTable, DimType dim);

/**
 * @brief helper function for data exchange
 */
template <typename FG_ELEMENT>
void sendAndReceiveIndices(std::vector<std::set<IndexType>>& send1dIndices,
                           std::vector<std::set<IndexType>>& recv1dIndices,
                           DistributedFullGrid<FG_ELEMENT>& dfg, DimType dim,
                           std::vector<RemoteDataContainer<FG_ELEMENT>>& remoteData) {
  std::vector<MPI_Request> sendRequests;
  std::vector<MPI_Request> recvRequests;
  size_t sendcount = 0;
  size_t recvcount = 0;
  for (size_t r = 0; r < send1dIndices.size(); ++r) {
    sendcount += send1dIndices[r].size();
  }
  for (size_t r = 0; r < recv1dIndices.size(); ++r) {
    recvcount += recv1dIndices[r].size();
  }

  sendRequests.resize(sendcount);
  recvRequests.resize(recvcount);

  IndexType totalSendSize(0);
  IndexType totalRecvSize(0);

  // for each rank r in send that has a nonempty index list
  sendcount = 0;
  for (size_t r = 0; r < send1dIndices.size(); ++r) {
    // for each index i in index list
    const auto& indices = send1dIndices[r];

    size_t k = 0;
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

        // compute size of subarray for stats
        {
          IndexType subarraySize = 1;

          for (DimType i = 0; i < subsizes.size(); ++i) subarraySize *= subsizes[i];

          totalSendSize += subarraySize;
        }

        // start
        IndexVector starts(dfg.getDimension(), 0);
        starts[dim] = lidxvec[dim];

        // note that if we want MPI to use c ordering, for less confusion as we
        // actually store our data in c format, we have to reverse all size and index
        // vectors
        std::vector<int> csizes(sizes.rbegin(), sizes.rend());
        std::vector<int> csubsizes(subsizes.rbegin(), subsizes.rend());
        std::vector<int> cstarts(starts.rbegin(), starts.rend());

        // create subarray view on data
        MPI_Type_create_subarray(static_cast<int>(dfg.getDimension()), &csizes[0], &csubsizes[0],
                                 &cstarts[0], MPI_ORDER_C, dfg.getMPIDatatype(), &mysubarray);
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
    sendcount += send1dIndices[r].size();
  }
  recvcount = 0;
  // IndexVector recv1dIndicesUnique;
  {
    // for each index in recv index list
    for (size_t r = 0; r < recv1dIndices.size(); ++r) {
      const auto& indices = recv1dIndices[r];
      const IndexVector& lowerBoundsNeighbor = dfg.getLowerBounds(static_cast<int>(r));

      size_t k = 0;
      for (const auto& index : indices) {
        // create RemoteDataContainer to store the subarray
        IndexVector sizes = dfg.getLocalSizes();
        sizes[dim] = 1;
        remoteData.emplace_back(sizes, dim, index, lowerBoundsNeighbor);
        // recv1dIndicesUnique.push_back(index);

        // start recv operation, use global index as tag
        {
          int src = static_cast<int>(r);
          int tag = static_cast<int>(index);

          FG_ELEMENT* buf = remoteData.back().getData();
          int bsize = static_cast<int>(remoteData.back().getSize());

          MPI_Irecv(buf, bsize, dfg.getMPIDatatype(), src, tag, dfg.getCommunicator(),
                    &recvRequests[recvcount + k++]);

          { totalRecvSize += bsize; }

#ifdef DEBUG_OUTPUT
          auto rank = dfg.getRank();
          // print info: dest, size, index
          std::cout << "rank " << rank << ": recv gindex " << index << " src " << src
                    << " size: " << bsize << std::endl;
#endif
        }
      }
      recvcount += recv1dIndices[r].size();
    }
  }
  // wait for finish of communication
  MPI_Waitall(static_cast<int>(sendRequests.size()), &sendRequests.front(), MPI_STATUSES_IGNORE);
  MPI_Waitall(static_cast<int>(recvRequests.size()), &recvRequests.front(), MPI_STATUSES_IGNORE);

  /*
  #ifdef DEBUG_OUTPUT
    MPI_Barrier( comm );

    // print array after data exchange
    {
      int size = dfg.getCommunicatorSize();

      for ( int r = 0; r < size; ++r ) {
        if ( r == rank ) {
          std::cout << "rank " << r << ":" << std::endl;
          if(dfg.getDimension() <= 3)
            std::cout << dfg;
        }

        MPI_Barrier(comm);
      }
    }

    // print buffers
    {
      int size = dfg.getCommunicatorSize();

      for ( int r = 0; r < size; ++r ) {
        if ( r == rank ) {
          std::cout << "rank " << r << ":" << std::endl;

          for ( size_t i = 0; i < remoteData.size(); ++i ) {
            assert( remoteData[i].getKeyIndex() == recv1dIndicesUnique[i] );
            std::cout << "\t" << recv1dIndicesUnique[i] << ": ";

            for ( IndexType k = 0; k < remoteData[i].getElementVector().size(); ++k)
              std::cout << "\t " << remoteData[i].getElementVector()[k];

            std::cout << std::endl;
          }
        }

        MPI_Barrier(comm);
      }
    }
  #endif */
}

// exchange data in dimension dim
template <typename FG_ELEMENT>
static void exchangeData1d(DistributedFullGrid<FG_ELEMENT>& dfg, DimType dim,
                           std::vector<RemoteDataContainer<FG_ELEMENT> >& remoteData) {

  auto commSize = dfg.getCommunicatorSize();
#ifdef DEBUG_OUTPUT
  CommunicatorType comm = dfg.getCommunicator();
  auto rank = dfg.getRank();
  std::vector<int> coords(dfg.getDimension());
  dfg.getPartitionCoords(coords);
  {
    std::cout << "in debug output" << std::endl;

    IndexType fidx = dfg.getFirstGlobal1dIndex(dim);
    LevelType flvl = dfg.getLevel(dim, fidx);
    IndexType fleftpre = dfg.getLeftPredecessor(dim, fidx);
    IndexType frightpre = dfg.getRightPredecessor(dim, fidx);
    RankType leftPreRank = getNeighbor1d(dfg, dim, fleftpre);
    RankType rightPreRank = getNeighbor1d(dfg, dim, frightpre);

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
    RankType leftPreRank = getNeighbor1d(dfg, dim, lleftpre);
    RankType rightPreRank = getNeighbor1d(dfg, dim, lrightpre);

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
  std::vector<std::set<IndexType>> recv1dIndices(commSize);
  std::vector<std::set<IndexType>> send1dIndices(commSize);

  // main loop
  IndexType idxMin = dfg.getFirstGlobal1dIndex(dim);
  IndexType idxMax = dfg.getLastGlobal1dIndex(dim);
  auto globalIdxMax = dfg.length(dim);
  LevelType lmax = dfg.getLevels()[dim];

  IndexType idx = idxMin;

  // mass-conserving stencil operations may also need neighbors
  bool exchangeParentsNeighbors = false; //TODO make work //TODO put this somewhere else
  // for hierarchization, we only need to exchange the direct predecessors
  while (idx <= idxMax) {
    LevelType lidx = dfg.getLevel(dim, idx);
    // check if successors of idx outside local domain
    {
      for (LevelType l = lidx + 1; l <= lmax; ++l) {
        LevelType ldiff = lmax - l;
        IndexType idiff = static_cast<IndexType>(std::pow(2, ldiff));

        // left successor, right successor
        for (const auto& indexShift : {-1, 1}) {
          IndexType sIdx = idx + indexShift * idiff;
          // if sIdx is outside of my domain, but still in the global domain
          if ((indexShift < 0 && sIdx >= 0 && sIdx < idxMin) ||
              (indexShift > 0 && sIdx > idxMax && sIdx < globalIdxMax)) {
            // get rank which has sIdx and add to send list
            int r = getNeighbor1d(dfg, dim, sIdx);
            if (r >= 0) send1dIndices[r].insert(idx);
          }
        }
      }
    }
    LevelType ldiff = lmax - lidx;
    IndexType idiff = static_cast<IndexType>(std::pow(2, ldiff));
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
        int r = getNeighbor1d(dfg, dim, pIdx);
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

  // if larger stencil: exchange with same-level neighbors
  if (exchangeParentsNeighbors) {
    for (LevelType lidx = 1; lidx <= lmax; ++lidx) {
      LevelType ldiff = lmax - lidx;
      IndexType idiff = static_cast<IndexType>(std::pow(2, ldiff + 1));
      // leftmost point of this level
      IndexType idx = getFirstIndexOfLevel1d(dfg, dim, lidx);
      // left neighbor
      IndexType nIdx = idx - idiff;
      assert(nIdx < idxMin);
      // if nIdx is in the global domain
      if (idx <= idxMax && nIdx >= 0) {
        // get rank which has same-level neighbor and add to list of indices to recv
        int r = getNeighbor1d(dfg, dim, nIdx);
        recv1dIndices[r].insert(nIdx);
        // it is mutual
        send1dIndices[r].insert(idx);
      }
      // rightmost point of this level
      idx = getLastIndexOfLevel1d(dfg, dim, lidx);
      // right neighbor
      nIdx = idx + idiff;
      assert(nIdx > idxMax);
      // if nIdx is in the global domain
      if (idx > 0 && nIdx < globalIdxMax) {
        // get rank which has same-level neighbor and add to list of indices to recv
        int r = getNeighbor1d(dfg, dim, nIdx);
        recv1dIndices[r].insert(nIdx);
        // it is mutual
        send1dIndices[r].insert(idx);
      }
    }
  }

#ifdef DEBUG_OUTPUT
  // print recvindices
  {
    for (int r = 0; r < commSize; ++r) {
      if (r == rank) {
        std::cout << "rank " << r << " recv1dIndices ";
        for (const auto& r : recv1dIndices) {
          std::cout << r;
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
          std::cout << s;
        }
        std::cout << std::endl;
      }
      MPI_Barrier(comm);
    }
  }

  MPI_Barrier(comm);

  /*
   if( dim == 4 ){
   std::string filename = "debug" + boost::lexical_cast<std::string>(rank);
   std::ofstream ofs( filename.c_str() );

   ofs << "rank " << rank << " recv: " << std::endl;
   for( RankType k = 0; k < dfg.getCommunicatorSize(); ++k ){
   if( recv1dIndices[k].size() > 0 )
   ofs << "\t" << k << ": " << recv1dIndices[k] << std::endl;
   }

   // print sendindices
   ofs << "rank " << rank << " send: " << std::endl;
   for( RankType k = 0; k < dfg.getCommunicatorSize(); ++k ){
   if( send1dIndices[k].size() > 0 )
   ofs << "\t" << k << ": " << send1dIndices[k] << std::endl;
   }

   ofs.close();
   }*/
#endif
  sendAndReceiveIndices(send1dIndices, recv1dIndices, dfg, dim, remoteData);
}


// exchange data in dimension dim
template <typename FG_ELEMENT>
static void exchangeData1dDehierarchization(
    DistributedFullGrid<FG_ELEMENT>& dfg, DimType dim,
    std::vector<RemoteDataContainer<FG_ELEMENT> >& remoteData) {

  auto commSize = dfg.getCommunicatorSize();
#ifdef DEBUG_OUTPUT
  CommunicatorType comm = dfg.getCommunicator();
  auto rank = dfg.getRank();
  std::vector<int> coords(dfg.getDimension());
  dfg.getPartitionCoords(coords);
  {
    IndexType fidx = dfg.getFirstGlobal1dIndex(dim);
    LevelType flvl = dfg.getLevel(dim, fidx);
    IndexType fleftpre = dfg.getLeftPredecessor(dim, fidx);
    IndexType frightpre = dfg.getRightPredecessor(dim, fidx);
    RankType leftPreRank = getNeighbor1d(dfg, dim, fleftpre);
    RankType rightPreRank = getNeighbor1d(dfg, dim, frightpre);

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
    RankType leftPreRank = getNeighbor1d(dfg, dim, lleftpre);
    RankType rightPreRank = getNeighbor1d(dfg, dim, lrightpre);

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
  std::vector<std::set<IndexType>> recv1dIndices(commSize);
  std::vector<std::set<IndexType>> send1dIndices(commSize);

  // main loop
  IndexType idxMin = dfg.getFirstGlobal1dIndex(dim);
  IndexType idxMax = dfg.getLastGlobal1dIndex(dim);
  // LevelType lmax = dfg.getLevels()[dim];

  IndexType idx = idxMin;

  // for dehierarchization, we need to exchange the full tree of predecessors
  while (idx <= idxMax) {
    checkLeftSuccesors(idx, idx, dim, dfg, send1dIndices);

    checkRightSuccesors(idx, idx, dim, dfg, send1dIndices);

    idx = checkPredecessors(idx, dim, dfg, recv1dIndices);
  }

#ifdef DEBUG_OUTPUT
  MPI_Barrier(comm);

  // print recvindices
  {
    for (int r = 0; r < commSize; ++r) {
      if (r == rank) {
        std::cout << "rank " << rank << " recv: " << std::endl;

        for (RankType k = 0; k < dfg.getCommunicatorSize(); ++k) {
          if (recv1dIndices[k].size() > 0)
            std::cout << "\t" << k << ": " << recv1dIndices[k] << std::endl;
        }
      }

      MPI_Barrier(comm);
    }
  }

  if (rank == 0) std::cout << std::endl;

  // print sendindices
  {
    for (int r = 0; r < commSize; ++r) {
      if (r == rank) {
        std::cout << "rank " << rank << " send: " << std::endl;

        for (RankType k = 0; k < dfg.getCommunicatorSize(); ++k) {
          if (send1dIndices[k].size() > 0)
            std::cout << "\t" << k << ": " << send1dIndices[k] << std::endl;
        }
      }

      MPI_Barrier(comm);
    }
  }

  MPI_Barrier(comm);

  /*
   if( dim == 4 ){
   std::string filename = "debug" + boost::lexical_cast<std::string>(rank);
   std::ofstream ofs( filename.c_str() );

   ofs << "rank " << rank << " recv: " << std::endl;
   for( RankType k = 0; k < dfg.getCommunicatorSize(); ++k ){
   if( recv1dIndices[k].size() > 0 )
   ofs << "\t" << k << ": " << recv1dIndices[k] << std::endl;
   }

   // print sendindices
   ofs << "rank " << rank << " send: " << std::endl;
   for( RankType k = 0; k < dfg.getCommunicatorSize(); ++k ){
   if( send1dIndices[k].size() > 0 )
   ofs << "\t" << k << ": " << send1dIndices[k] << std::endl;
   }

   ofs.close();
   }*/
#endif

  sendAndReceiveIndices(send1dIndices, recv1dIndices, dfg, dim, remoteData);
}

template <typename FG_ELEMENT>
void checkLeftSuccesors(IndexType checkIdx, IndexType rootIdx, DimType dim,
                        DistributedFullGrid<FG_ELEMENT>& dfg,
                        std::vector<std::set<IndexType>>& OneDIndices) {
  LevelType lidx = dfg.getLevel(dim, checkIdx);
  IndexType idxMin = dfg.getFirstGlobal1dIndex(dim);
  // IndexType idxMax = dfg.getLastGlobal1dIndex(dim);
  LevelType lmax = dfg.getLevels()[dim];

  // check left successors of checkIdx
  for (LevelType l = lidx + 1; l <= lmax; ++l) {
    LevelType ldiff = lmax - l;
    IndexType idiff = static_cast<IndexType>(std::pow(2, ldiff));

    IndexType lsIdx = checkIdx - idiff;

    if (lsIdx >= 0 && lsIdx < idxMin) {
      // get rank which has lsIdx and add to send list
      int r = getNeighbor1d(dfg, dim, lsIdx);
      if (r >= 0) OneDIndices[r].insert(rootIdx);
    }

    if (lsIdx >= 0) checkLeftSuccesors(lsIdx, rootIdx, dim, dfg, OneDIndices);
  }
}

template <typename FG_ELEMENT>
void checkRightSuccesors(IndexType checkIdx, IndexType rootIdx, DimType dim,
                         DistributedFullGrid<FG_ELEMENT>& dfg,
                         std::vector<std::set<IndexType>>& OneDIndices) {
  LevelType lidx = dfg.getLevel(dim, checkIdx);

  IndexType idxMax = dfg.getLastGlobal1dIndex(dim);
  LevelType lmax = dfg.getLevels()[dim];

  // check right successors of checkIdx
  for (LevelType l = lidx + 1; l <= lmax; ++l) {
    LevelType ldiff = lmax - l;
    IndexType idiff = static_cast<IndexType>(std::pow(2, ldiff));

    IndexType rsIdx = checkIdx + idiff;

    if (rsIdx < dfg.getGlobalSizes()[dim] && rsIdx > idxMax) {
      // get rank which has rsIdx and add to send list
      int r = getNeighbor1d(dfg, dim, rsIdx);
      if (r >= 0) OneDIndices[r].insert(rootIdx);
    }

    if (rsIdx < dfg.length(dim))
      checkRightSuccesors(rsIdx, rootIdx, dim, dfg, OneDIndices);
  }
}

template <typename FG_ELEMENT>
IndexType checkPredecessors(IndexType idx, DimType dim, DistributedFullGrid<FG_ELEMENT>& dfg,
                            std::vector<std::set<IndexType>>& OneDIndices) {
  IndexType idxMin = dfg.getFirstGlobal1dIndex(dim);
  IndexType idxMax = dfg.getLastGlobal1dIndex(dim);

  // check if left predecessor outside local domain
  // if returns negative value there's no left predecessor
  IndexType lpIdx = dfg.getLeftPredecessor(dim, idx);

  if (lpIdx >= 0 && lpIdx < idxMin) {
    // get rank which has left predecessor and add to list of indices
    int r = getNeighbor1d(dfg, dim, lpIdx);
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
    // get rank which has left predecessor and add to list of indices to recv
    int r = getNeighbor1d(dfg, dim, rpIdx);
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
RankType getNeighbor1d(DistributedFullGrid<FG_ELEMENT>& dfg, DimType dim, IndexType idx1d) {
  // if global index is outside of domain return negative value
  {
    if (idx1d < 0) return -1;

    IndexType numElementsD = dfg.getGlobalSizes()[dim];

    if (idx1d > numElementsD - 1) return -1;
  }

  IndexVector globalAxisIndex = dfg.getLowerBounds();
  globalAxisIndex[dim] = idx1d;

  IndexVector partitionCoords(dfg.getDimension());
  dfg.getPartitionCoords(globalAxisIndex, partitionCoords);
  RankType r = dfg.getRank(partitionCoords);

  // check if global index vector is actually contained in the domain of rank r
  assert(globalAxisIndex >= dfg.getLowerBounds(r));
  assert(globalAxisIndex < dfg.getUpperBounds(r));
  assert(r < dfg.getCommunicatorSize());
  return r;
}

// returns the next one-dimensional global index which fulfills
// min( lmax, l(idx1d) + 1 )
template <typename FG_ELEMENT>
IndexType getNextIndex1d(DistributedFullGrid<FG_ELEMENT>& dfg, DimType dim, IndexType idx1d) {
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
static IndexType getFirstIndexOfLevel1d(DistributedFullGrid<FG_ELEMENT>& dfg, DimType dim,
                                        LevelType l) {
  IndexType idxMax = dfg.getLastGlobal1dIndex(dim);
  IndexType idxMin = dfg.getFirstGlobal1dIndex(dim);

  for (IndexType i = idxMin; i <= idxMax; ++i) {
    if (dfg.getLevel(dim, i) == l) return i;
  }

  // no index on level l found, return value which is out of local index range
  return idxMax + 1;
}

template <typename FG_ELEMENT>
static IndexType getLastIndexOfLevel1d(DistributedFullGrid<FG_ELEMENT>& dfg, DimType dim,
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
static void hierarchizeX(DistributedFullGrid<FG_ELEMENT>& dfg,
                         LookupTable<FG_ELEMENT>& lookupTable) {
  const DimType dim = 0;

  LevelType lmax = dfg.getLevels()[dim];
  IndexType idxMax = dfg.getLastGlobal1dIndex(dim);
  IndexType ndim = dfg.getLocalSizes()[dim];

  FG_ELEMENT zeroVal(0);

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

    IndexVector globalIndexVectorCenter = baseGlobalIndexVector;
    IndexVector globalIndexVectorLeft = baseGlobalIndexVector;
    IndexVector globalIndexVectorRight = baseGlobalIndexVector;

    for (LevelType l = lmax; l > 0; --l) {
      // get first local point of level and corresponding stride
      IndexType firstOfLevel = getFirstIndexOfLevel1d(dfg, dim, l);
      IndexType parentOffset = static_cast<IndexType>(std::pow(2, lmax - l));
      IndexType levelStride = parentOffset * 2;

      // loop over points of this level with level specific stride
      // as long as inside domain
      for (IndexType idx = firstOfLevel; idx <= idxMax; idx += levelStride) {
        // compute global index vector of center, left and right predecessor
        globalIndexVectorCenter[dim] = idx;
        globalIndexVectorLeft[dim] = idx - parentOffset;
        globalIndexVectorRight[dim] = idx + parentOffset;

        // translate global indices vector to pointers
        if (dfg.returnBoundaryFlags()[dim] == true) {
          FG_ELEMENT* center = lookupTable.getData(globalIndexVectorCenter);
          FG_ELEMENT* left = lookupTable.getData(globalIndexVectorLeft);
          FG_ELEMENT* right = lookupTable.getData(globalIndexVectorRight);

          // do calculation
          *center -= 0.5 * (*left + *right);
        } else {
          FG_ELEMENT* center = lookupTable.getData(globalIndexVectorCenter);

          // when no boundary in this dimension we have to check if
          // 1d indices outside domain
          FG_ELEMENT* left = &zeroVal;
          FG_ELEMENT* right = &zeroVal;

          if (globalIndexVectorLeft[dim] > 0) {
            left = lookupTable.getData(globalIndexVectorLeft);
          }

          if (globalIndexVectorRight[dim] < dfg.getGlobalSizes()[dim]) {
            right = lookupTable.getData(globalIndexVectorRight);
          }

          // do calculation
          *center -= 0.5 * (*left + *right);
        }
      }
    }
  }
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
      IndexType parentOffset = static_cast<IndexType>(std::pow(2, lmax - l));
      IndexType levelStride = parentOffset * 2;

      // loop over points of this level with level specific stride
      // as long as inside domain
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
      IndexType parentOffset = static_cast<IndexType>(std::pow(2, lmax - l));
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

/**
 * @brief hierarchize a DFG in dimension X (with contiguous access)
 * 
 * @param dfg : the DFG to hierarchize
 * @param lookupTable: the lookup table for local and remote data
 */
template <typename FG_ELEMENT>
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

  // first global index for hierarchization kernel
  // first nonboundary point
  IndexType idxstart = gstart;

  if (gstart == 0) idxstart += 1;

  // last global index inside subdomain.
  // IndexType idxend = dfg.getUpperBounds()[dim] - 1;

  // level of gend
  // LevelType level_idxend = dfg.getLevel(dim, idxend);
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

    hierarchizeX_opt_boundary_kernel(&tmp[0], lmax, 0, 1);

    // copy local data back
    for (IndexType i = 0; i < xSize; ++i) localData[linIdxBlockStart + i] = tmp[gstart + i];
  }
}

template <typename FG_ELEMENT>
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

  // first global index for hierarchization kernel
  // first nonboundary point
  IndexType idxstart = gstart;

  if (gstart == 0) idxstart += 1;

  // last global index inside subdomain.
  // IndexType idxend = dfg.getUpperBounds()[dim] - 1;

  // level of gend
  // LevelType level_idxend = dfg.getLevel(dim, idxend);
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

    dehierarchizeX_opt_boundary_kernel(&tmp[0], lmax, 0, 1);

    // copy local data back
    for (IndexType i = 0; i < xSize; ++i) localData[linIdxBlockStart + i] = tmp[gstart + i];
  }
}

template <typename FG_ELEMENT>
inline void hierarchizeX_opt_boundary_kernel(FG_ELEMENT* data, LevelType lmax, int start,
                                             int stride) {
  const int lmaxi = static_cast<int>(lmax);
  int ll = lmaxi;
  int steps = (1 << (lmaxi - 1));
  int offset = 1;  // 1 and not 0 because boundary
  int stepsize = 2;
  int parentOffset = 1;

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

template <typename FG_ELEMENT>
inline void dehierarchizeX_opt_boundary_kernel(FG_ELEMENT* data, LevelType lmax, int start,
                                               int stride) {
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


template <typename FG_ELEMENT>
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

  // first global index for hierarchization kernel. may not be a boundary point
  IndexType idxstart = gstart;

  if (gstart == 0) idxstart += 1;

  // last global index inside subdomain and corresponding level
  // IndexType idxend = dfg.getUpperBounds()[dim] - 1;
  // LevelType level_idxend = dfg.getLevel(dim, idxend);

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
    hierarchizeX_opt_boundary_kernel(&tmp[0], lmax, 0, 1);

    // copy pole back
    for (IndexType i = 0; i < ndim; ++i) ldata[start + stride * i] = tmp[gstart + i];
  }
}

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

  // first global index for hierarchization kernel. may not be a boundary point
  IndexType idxstart = gstart;

  if (gstart == 0) idxstart += 1;

  // last global index inside subdomain and corresponding level
  // IndexType idxend = dfg.getUpperBounds()[dim] - 1;
  // LevelType level_idxend = dfg.getLevel(dim, idxend);

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
      IndexType parentOffset = static_cast<IndexType>(std::pow(2, lmax - l));
      IndexType levelStride = parentOffset * 2;

      // loop over points of this level with level specific stride
      // as long as inside domain
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

    // copy pole back
    for (IndexType i = 0; i < ndim; ++i) ldata[start + stride * i] = tmp[gstart + i];
  }
}

template <typename FG_ELEMENT>
void dehierarchizeN_opt_boundary(DistributedFullGrid<FG_ELEMENT>& dfg,
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

  // first global index for hierarchization kernel. may not be a boundary point
  IndexType idxstart = gstart;

  if (gstart == 0) idxstart += 1;

  // last global index inside subdomain and corresponding level
  // IndexType idxend = dfg.getUpperBounds()[dim] - 1;
  // LevelType level_idxend = dfg.getLevel(dim, idxend);

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
    dehierarchizeX_opt_boundary_kernel(&tmp[0], lmax, 0, 1);

    // copy pole back
    for (IndexType i = 0; i < ndim; ++i) ldata[start + stride * i] = tmp[gstart + i];
  }
}

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

  // first global index for hierarchization kernel. may not be a boundary point
  IndexType idxstart = gstart;

  if (gstart == 0) idxstart += 1;

  // last global index inside subdomain and corresponding level
  // IndexType idxend = dfg.getUpperBounds()[dim] - 1;
  // LevelType level_idxend = dfg.getLevel(dim, idxend);

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
      IndexType parentOffset = static_cast<IndexType>(std::pow(2, lmax - l));
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

// template<typename FG_ELEMENT>
class DistributedHierarchization {
 public:
  // inplace hierarchization
  template <typename FG_ELEMENT>
  static void hierarchize(DistributedFullGrid<FG_ELEMENT>& dfg, const std::vector<bool>& dims) {
    assert(dfg.getDimension() > 0);
    assert(dfg.getDimension() == dims.size());

    // hierarchize first dimension
    if (dims[0]) {
      DimType dim = 0;

      // exchange data first dimension
      std::vector<RemoteDataContainer<FG_ELEMENT> > remoteData;
      exchangeData1d(dfg, dim, remoteData);

      LookupTable<FG_ELEMENT> lookupTable(remoteData, dfg, dim);

      if (dfg.returnBoundaryFlags()[dim] == true) {
        hierarchizeX_opt_boundary(dfg, lookupTable);
      } else {
        hierarchizeX_opt_noboundary(dfg, lookupTable);
      }
    }

    // hierarchize other dimensions
    for (DimType dim = 1; dim < dfg.getDimension(); ++dim) {
      if (!dims[dim]) continue;

      // exchange data
      std::vector<RemoteDataContainer<FG_ELEMENT> > remoteData;
      exchangeData1d(dfg, dim, remoteData);
      LookupTable<FG_ELEMENT> lookupTable(remoteData, dfg, dim);

      if (dfg.returnBoundaryFlags()[dim] == true) {
        hierarchizeN_opt_boundary(dfg, lookupTable, dim);
      } else {
        hierarchizeN_opt_noboundary(dfg, lookupTable, dim);
      }
    }
  }

  template <typename FG_ELEMENT>
  static void hierarchize(DistributedFullGrid<FG_ELEMENT>& dfg) {
    std::vector<bool> dims(dfg.getDimension(), true);
    hierarchize<FG_ELEMENT>(dfg, dims);
  }

  // inplace dehierarchization
  template <typename FG_ELEMENT>
  static void dehierarchize(DistributedFullGrid<FG_ELEMENT>& dfg, const std::vector<bool>& dims) {
    assert(dfg.getDimension() > 0);
    assert(dfg.getDimension() == dims.size());

    // dehierarchize first dimension
    if (dims[0]) {
      DimType dim = 0;

      // exchange data first dimension
      std::vector<RemoteDataContainer<FG_ELEMENT> > remoteData;
      exchangeData1dDehierarchization(dfg, dim, remoteData);
      LookupTable<FG_ELEMENT> lookupTable(remoteData, dfg, dim);

      if (dfg.returnBoundaryFlags()[dim] == true) {
        dehierarchizeX_opt_boundary(dfg, lookupTable);
      } else {
        dehierarchizeX_opt_noboundary(dfg, lookupTable);
      }
    }

    // dehierarchize other dimensions
    for (DimType dim = 1; dim < dfg.getDimension(); ++dim) {
      if (!dims[dim]) continue;

      // exchange data
      std::vector<RemoteDataContainer<FG_ELEMENT> > remoteData;
      exchangeData1dDehierarchization(dfg, dim, remoteData);
      LookupTable<FG_ELEMENT> lookupTable(remoteData, dfg, dim);

      if (dfg.returnBoundaryFlags()[dim] == true) {
        dehierarchizeN_opt_boundary(dfg, lookupTable, dim);
      } else {
        dehierarchizeN_opt_noboundary(dfg, lookupTable, dim);
      }
    }
  }

  template <typename FG_ELEMENT>
  static void dehierarchize(DistributedFullGrid<FG_ELEMENT>& dfg) {
    std::vector<bool> dims(dfg.getDimension(), true);
    dehierarchize<FG_ELEMENT>(dfg, dims);
  }

  template <typename FG_ELEMENT>
  static void fillDFGFromDSGU(DistributedFullGrid<FG_ELEMENT>& dfg,
      DistributedSparseGridUniform<FG_ELEMENT>& dsg, const std::vector<bool>& hierarchizationDims){
    // fill dfg with hierarchical coefficients from distributed sparse grid
    dfg.extractFromUniformSG(dsg);

    // dehierarchize dfg
    DistributedHierarchization::dehierarchize<FG_ELEMENT>(
        dfg, hierarchizationDims);
  }
};
// class DistributedHierarchization

}  // namespace combigrid

#endif /* DistributedHierarchization_HPP_ */
