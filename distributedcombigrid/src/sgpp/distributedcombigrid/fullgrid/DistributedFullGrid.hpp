#ifndef DISTRIBUTEDCOMBIFULLGRID_HPP_
#define DISTRIBUTEDCOMBIFULLGRID_HPP_

#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>
#include <string>

#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/DistributedSparseGridUniform.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/SGrid.hpp"
#include "sgpp/distributedcombigrid/utils/IndexVector.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/utils/Stats.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"

//#define DEBUG_OUTPUT
#define UNIFORM_SG

namespace combigrid {

// // a struct to hold information on the hierarchical subspaces contained in this grid
// template <typename FG_ELEMENT>
// struct SubspaceDFG {
//   // the hierarchical level
//   LevelVector level_;

//   // size of the subspace, in each dimension
//   IndexVector sizes_;

//   // // the data in this subspace
//   // std::vector<FG_ELEMENT> data_;

//   // // the total number of grid points in this subspace
//   // size_t targetSize_;

//   // // the number of grid points resident in this partition of the grid
//   // size_t localSize_;
// };

/** The full grid class which is the main building block of the combi grid <br>
 *  The index of a gridpoint in the full grid is given by the formula : <br>
 *  ind = i0 + i1*N0 + i2*N0*N1 + ... + id*N0*N1*N2*...*Nd, where i0,i1,i2,... are the indexes in
 * every dimension,
 *  and Nk is the number of gridpoints in direction k. <br>
 *  For more infos you can look at the "getVectorIndex" and "getLinearIndex" functions and their
 * implementation. <br>
 *  Nk=2^level[k]+1 for every k for the directions with boundary points and Nk=2^level[k]-1 for the
 * directions without boundary. <br>
 *
 *  The layout of the process grid, i.e. the ordering of procs is the same as of
 *  levels. The same holds for every function involving partition coordinates, e.g.
 *  getRank().
 *  However, in the internal implementation, i.e. for the grid of MPI processes,
 *  an inverse ordering for the process coordinates is used. This must be kept
 *  in mind when mapping an application's data structure to the distributed
 *  full grid.
 *  In future, we plan to offer more flexiblity here. At least to cover the two
 *  most natural options to arrange processes in a grid. At the end of the day
 *  this only affects the rank in the cartesian communicator, which is assigned
 *  by MPI in fortran-typical column-major ordering.
 * */
template <typename FG_ELEMENT>
class DistributedFullGrid {
 public:
  /** dimension adaptive Ctor */
  DistributedFullGrid(DimType dim, const LevelVector& levels, CommunicatorType comm,
                      const std::vector<bool>& hasBdrPoints, const IndexVector& procs,
                      bool forwardDecomposition = true,
                      const std::vector<IndexVector>& decomposition = std::vector<IndexVector>(),
                      const BasisFunctionBasis* basis = NULL)
      : dim_(dim), levels_(levels), procs_(procs) {
    assert(levels.size() == dim);
    assert(hasBdrPoints.size() == dim);
    assert(procs.size() == dim);
    hasBoundaryPoints_ = hasBdrPoints;

    InitMPI(comm);  // will also check grids per dim

    // set global num of elements and offsets
    nrElements_ = 1;
    offsets_.resize(dim_);
    nrPoints_.resize(dim_);

    for (DimType j = 0; j < dim_; j++) {
      nrPoints_[j] = ((hasBoundaryPoints_[j] == true) ? (powerOfTwo[levels_[j]] + 1)
                                                      : (powerOfTwo[levels_[j]] - 1));
      offsets_[j] = nrElements_;
      nrElements_ = nrElements_ * nrPoints_[j];
    }

    if (decomposition.size() == 0) {
      decomposition_ = getDefaultDecomposition(forwardDecomposition);
    } else {
      setDecomposition(decomposition);
    }
    myPartitionsLowerBounds_ = getLowerBounds(rank_);
    myPartitionsUpperBounds_ = getUpperBounds(rank_);

    // set local elements and local offsets
    nrLocalPoints_ = getUpperBounds() - getLowerBounds();

    nrLocalElements_ = 1;
    localOffsets_.resize(dim);

    for (DimType j = 0; j < dim_; ++j) {
      localOffsets_[j] = nrLocalElements_;
      nrLocalElements_ *= nrLocalPoints_[j];
    }

    // in contrast to serial implementation we directly create the grid
    fullgridVector_.resize(nrLocalElements_);

    dsg_ = nullptr;

#ifdef DEBUG_OUTPUT
    if (rank_ == 0) {
      for (RankType r = 0; r < size_; ++r) {
        std::cout << "rank " << r << ": \n"
                  << "\t lower bounds " << getLowerBounds(r)
                  << "\t upper bounds " << getUpperBounds(r)
                  << std::endl;
      }
    }
#endif
  }

  // copy construction would need to duplicate communicator_
  DistributedFullGrid(const DistributedFullGrid& other) = delete;
  DistributedFullGrid& operator=( const DistributedFullGrid & ) = delete;
  DistributedFullGrid(DistributedFullGrid&& other) = delete;
  DistributedFullGrid& operator=(DistributedFullGrid&& other) = delete;

  virtual ~DistributedFullGrid() {
    // todo: remove communicators? Yes -> Done
    MPI_Comm_free(&communicator_);
    for (size_t i = 0; i < upwardSubarrays_.size(); ++i)  {
      MPI_Type_free(&upwardSubarrays_[i]);
    }
    for (size_t i = 0; i < downwardSubarrays_.size(); ++i)  {
      MPI_Type_free(&downwardSubarrays_[i]);
    }
  }

  FG_ELEMENT evalLocalIndexOn(const IndexVector& localIndex, const std::vector<real>& coords) const {
    auto firstIndex = IndexVector(dim_, 0);
    auto lastIndex = this->getLastGlobalIndex() - this->getFirstGlobalIndex();
    // if this local index is out of bounds, return 0. (will be contributed by other partial dfg)
    if (! (localIndex >= firstIndex && localIndex <= lastIndex) ) {
      // std::cout << "out of bounds" << localIndex << firstIndex << lastIndex << std::endl;
      return 0.;
    }

    // get coords corresponding to localIndex
    auto localLinearIndex = getLocalLinearIndex(localIndex);
    std::vector<real> pointCoords (this->getDimension());
    getCoordsLocal(localLinearIndex, pointCoords);

    // get product of 1D hat functions on coords
    auto h = getGridSpacing();
    real phi_c = 1.; // value of product of basis function on coords
    for (DimType d = 0 ; d < dim_ ; ++d){
      // get distance between coords and point
      pointCoords[d] -= coords[d];
      if (std::abs(pointCoords[d]) > h[d]){
        std::cout << "assert bounds " << pointCoords << coords <<
         h << d << localIndex << lastIndex << std::endl;
        assert(false &&
          "should only be called for coordinates within the support of this point's basis function");
      }
      phi_c *= 1. - std::abs(pointCoords[d]/h[d]);
    }
    // std::cout << "coords " <<  localIndex << coords << localLinearIndex << h << std::endl;
    // std::cout << "phi_c " << phi_c << this->getElementVector()[localLinearIndex] << std::endl;
    assert(phi_c >= 0.);
    return phi_c * this->getElementVector()[localLinearIndex];
  }

  /**
   * @brief recursive call to evaluate all neighbor points' contributions to the coordinate (on this part of the grid)
   *
   * @param localIndex the (in-all-dimensions lower) neighbor of coords
   * @param dim the current dimension to split on (start with 0)
   * @param coords the coordinate to interpolate on
   * @return FG_ELEMENT the interpolated value at coords
   */
  FG_ELEMENT evalMultiindexRecursively (const IndexVector& localIndex, DimType dim, const std::vector<real>& coords) const {
    assert(!(dim > this->getDimension()));
    if (dim == this->getDimension()){
      // std::cout << "eval " << localIndex << std::endl;
      return evalLocalIndexOn(localIndex, coords);
    } else {
      FG_ELEMENT sum = 0.;
      IndexVector localIndexDimPlusOne = localIndex;
      localIndexDimPlusOne[dim] += 1;
      // std::cout << localIndex << localIndexDimPlusOne << std::endl;
      sum += evalMultiindexRecursively(localIndex, dim+1, coords);
      sum += evalMultiindexRecursively(localIndexDimPlusOne, dim+1, coords);
      return sum;
    }
  }

  /** evaluates the full grid on the specified coordinates
   * @param coords ND coordinates on the unit square [0,1]^D*/
  FG_ELEMENT eval(const std::vector<real>& coords) const {
    FG_ELEMENT value;
    eval(coords, value);
    return value;
  }

  void eval(const std::vector<real>& coords, FG_ELEMENT& value, MPI_Request* request = nullptr) const {
    assert(coords.size() == this->getDimension());

    // get the lowest-index point of the points
    // whose basis functions contribute to the interpolated value
    auto lowerCoords = getLowerBoundsCoords();
    auto h = getGridSpacing();
    IndexVector localIndexLowerNonzeroNeighborPoint (dim_);
    for (DimType d = 0 ; d < dim_ ; ++d){
      assert(coords[d] >= 0. && coords[d] <= 1.);
      localIndexLowerNonzeroNeighborPoint[d] = std::floor((coords[d] - lowerCoords[d]) / h[d]);
    }
    // std::cout <<localIndexLowerNonzeroNeighborPoint << coords << lowerCoords << h << std::endl;

    // evaluate at those points and sum up according to the basis function
    // needs to be recursive in order to be dimensionally adaptive
    value = evalMultiindexRecursively(localIndexLowerNonzeroNeighborPoint, 0, coords);

    if (request == nullptr) {
      MPI_Allreduce(MPI_IN_PLACE, &value, 1, this->getMPIDatatype(), MPI_SUM, this->getCommunicator());
    } else {
      MPI_Iallreduce(MPI_IN_PLACE, &value, 1, this->getMPIDatatype(), MPI_SUM, this->getCommunicator(), request);
    }
  }

  /** evaluates the full grid on the specified coordinates
   * @param interpolationCoords vector of ND coordinates on the unit square [0,1]^D*/
  std::vector<FG_ELEMENT> getInterpolatedValues(const std::vector<std::vector<real>>& interpolationCoords) const {
    auto numValues = interpolationCoords.size();
    std::vector<FG_ELEMENT> values;
    values.resize(numValues);
    // using the asynchronous communication makes the runtime quadratic in numValues -> synchronous MPI for now
    for (size_t i = 0; i < numValues; ++i) {
      this->eval(interpolationCoords[i], values[i]);//, &requests[i]);
    }
    return values;
  }

  /** return the coordinates on the unit square corresponding to global idx
   * @param globalIndex [IN] global linear index of the element i
   * @param coords [OUT] the vector must be resized already */
  inline void getCoordsGlobal(IndexType globalLinearIndex, std::vector<real>& coords) const {
    // temporary variables
    // int verb = 6;
    IndexType ind = 0;
    IndexType tmp_add = 0;

    coords.resize(dim_);

    for (DimType j = 0; j < dim_; j++) {
      ind = globalLinearIndex % nrPoints_[j];
      globalLinearIndex = globalLinearIndex / nrPoints_[j];
      // set the coordinate based on if we have boundary points
      tmp_add = (hasBoundaryPoints_[j] == true) ? (0) : (1);
      coords[j] = static_cast<double>(ind + tmp_add) * getGridSpacing()[j];
    }
  }

  /** return coordinates on the unit square corresponding to local idx
   * @param globalIndex [IN] local linear index of the element i
   * @param coords [OUT] the vector must be resized already */
  inline void getCoordsLocal(IndexType localLinearIndex, std::vector<real>& coords) const {
    // todo: probably very inefficient implementation, if crucial for
    // performance implement more direct way of computing the coordinates

    assert(localLinearIndex < getNrLocalElements());

    IndexType globalLinearIndex = getGlobalLinearIndex(localLinearIndex);

    getCoordsGlobal(globalLinearIndex, coords);
  }

  /** returns the LI (level,index) notation for a given element in the full grid
   * @param elementIndex [IN] the linear index of the element
   * @param levels [OUT] the levels of the point in the LI notation
   * @param indexes [OUT] the indexes of the point in the LI notation */
  inline void getGlobalLI(IndexType elementIndex, LevelVector& levels, IndexVector& indexes) const {
    IndexType startindex, tmp_val;

    assert(elementIndex < nrElements_);
    levels.resize(dim_);
    indexes.resize(dim_);

    tmp_val = elementIndex;

    // first calculate intermediary indexes
    for (DimType k = 0; k < dim_; k++) {
      startindex = (hasBoundaryPoints_[k]) ? 0 : 1;
      indexes[k] = tmp_val % nrPoints_[k] + startindex;
      tmp_val = tmp_val / nrPoints_[k];
    }

    // The level and index of the element in the hashgridstorage are computed dividing by two the
    // index and level in the fullgrid
    // until we obtain an impair number for the index, thus obtaining the level and index in the
    // hierarchical basis (Aliz Nagy)
    // ...
    for (DimType k = 0; k < dim_; k++) {
      tmp_val = levels_[k];

      if (indexes[k] != 0) {
        // todo: these operations can be optimized
        while (indexes[k] % 2 == 0) {
          indexes[k] = indexes[k] / 2;
          tmp_val--;
        }
      } else {
        tmp_val = 0;
      }

      levels[k] = tmp_val;
    }
  }

  /** the global vector index corresponding to a global linear index
   * @param linIndex [IN] the linear index
   * @param axisIndex [OUT] the returned vector index */
  inline void getGlobalVectorIndex(IndexType globLinIndex, IndexVector& globAxisIndex) const {
    assert(globLinIndex < nrElements_);
    assert(globAxisIndex.size() == dim_);

    IndexType tmp = globLinIndex;

    for (int i = static_cast<int>(dim_) - 1; i >= 0; i--) {
      globAxisIndex[i] = tmp / (this->getOffset(i));
      tmp = tmp % this->getOffset(i);
    }
  }

  /** the global vector index corresponding to a local vector index */
  inline void getGlobalVectorIndex(const IndexVector& locAxisIndex,
                                   IndexVector& globAxisIndex) const {
    assert(locAxisIndex.size() == dim_);

    globAxisIndex = this->getLowerBounds() + locAxisIndex;
  }

  /** the local vector index corresponding to a local linear index */
  inline void getLocalVectorIndex(IndexType locLinIndex, IndexVector& locAxisIndex) const {
    assert(locLinIndex < nrElements_);
    assert(locAxisIndex.size() == dim_);

    IndexType tmp = locLinIndex;

    for (int i = static_cast<int>(dim_) - 1; i >= 0; i--) {
      locAxisIndex[i] = tmp / localOffsets_[i];
      tmp = tmp % localOffsets_[i];
    }
  }

  /** the local vector index corresponding to a global vector index
   if global index vector not contained in local domain false will be
   returned */
  inline bool getLocalVectorIndex(const IndexVector& globAxisIndex,
                                  IndexVector& locAxisIndex) const {
    assert(globAxisIndex.size() == dim_);

    if (globAxisIndex >= getLowerBounds() && globAxisIndex < getUpperBounds()) {
      locAxisIndex = globAxisIndex - this->getLowerBounds();
      return true;
    } else {
      return false;
    }
  }

  /** returns the global linear index corresponding to the global index vector
   * @param axisIndex [IN] the vector index */
  inline IndexType getGlobalLinearIndex(const IndexVector& globAxisIndex) const {
    assert(globAxisIndex.size() == dim_);

    IndexType tmp = 0;

    for (int i = static_cast<int>(dim_) - 1; i >= 0; i--) {
      tmp = tmp + offsets_[i] * globAxisIndex[i];
    }

    assert(tmp < nrElements_);

    return tmp;
  }

  // the global linear index corresponding to the local linear index
  inline IndexType getGlobalLinearIndex(IndexType locLinIndex) const {
    assert(locLinIndex < nrLocalElements_);

    // convert to local vector index
    IndexVector locAxisIndex(dim_);
    getLocalVectorIndex(locLinIndex, locAxisIndex);

    // convert to global vector index
    IndexVector globAxisIndex(dim_);
    getGlobalVectorIndex(locAxisIndex, globAxisIndex);

    // convert to global linear index
    IndexType globLinIndex = getGlobalLinearIndex(globAxisIndex);

    assert(globLinIndex < nrElements_);

    return globLinIndex;
  }

  /** the local linear index corresponding to the local index vector */
  inline IndexType getLocalLinearIndex(const IndexVector& locAxisIndex) const {
    assert(locAxisIndex.size() == dim_);

    IndexType tmp = 0;

    for (DimType i = 0; i < dim_; ++i) {
      tmp = tmp + localOffsets_[i] * locAxisIndex[i];
    }

    assert(tmp < nrLocalElements_);

    return tmp;
  }

  /** the global linear index corresponding to the local linear index
   returns negative value if element not inside local partition */
  inline IndexType getLocalLinearIndex(IndexType globLinIndex) const {
    assert(globLinIndex < nrElements_);

    // convert to global vector index
    IndexVector globAxisIndex(dim_);
    getGlobalVectorIndex(globLinIndex, globAxisIndex);

    // convert to local vector index
    IndexVector locAxisIndex(dim_);

    if (getLocalVectorIndex(globAxisIndex, locAxisIndex)) {
      // convert to local linear index
      return getLocalLinearIndex(locAxisIndex);
    } else {
      return -1;
    }
  }

  /** returns the dimension of the full grid */
  inline DimType getDimension() const { return dim_; }

  /** the getters for the full grid vector */
  inline std::vector<FG_ELEMENT>& getElementVector() { return fullgridVector_; }

  inline const std::vector<FG_ELEMENT>& getElementVector() const { return fullgridVector_; }

  /** return the offset in the full grid vector of the dimension */
  inline IndexType getOffset(DimType i) const { return offsets_[i]; }

  inline const IndexVector& getOffsets() const { return offsets_; }

  inline const IndexVector& getLocalOffsets() const { return localOffsets_; }

  /** return the level vector */
  inline const LevelVector& getLevels() const { return levels_; }

  /** returns the grid spacing (sometimes called h) */
  std::vector<double> getGridSpacing() const {
    std::vector<double> h(dim_);
    for (DimType d = 0 ; d < dim_ ; ++d){
      h[d] = oneOverPowOfTwo[levels_[d]];
    }
    return h;
  }

  /** returns the number of elements in the full grid */
  inline IndexType getNrElements() const { return nrElements_; }

  /** number of elements in the local partition */
  inline IndexType getNrLocalElements() const { return nrLocalElements_; }

  /** number of points per dimension i */
  inline IndexType length(DimType i) const { return nrPoints_[i]; }

  /** vector of flags to show if the dimension has boundary points*/
  inline const std::vector<bool>& returnBoundaryFlags() const { return hasBoundaryPoints_; }

  // /** copies the input vector to the full grid vector
  //  * @param in [IN] input vector*/
  // void setElementVector(const std::vector<FG_ELEMENT>& in) {
  //   assert(in.size() == static_cast<IndexType>(nrElements_));
  //   fullgridVector_ = in;
  // }

  inline void setZero(){
    for (auto& element : this->getElementVector()){
      element = 0.;
    }
  }

  inline FG_ELEMENT* getData() { return &fullgridVector_[0]; }

  inline const FG_ELEMENT* getData() const { return &fullgridVector_[0]; }

  /** MPI Communicator*/
  inline CommunicatorType getCommunicator() const { return communicator_; }

  inline RankType getRank() const { return rank_; }

  inline int getCommunicatorSize() const { return size_; }

  /** lower Bounds of this process */
  inline IndexVector getLowerBounds() const {
    // return getLowerBounds(rank_);
    return myPartitionsLowerBounds_;
  }

  /** lower bounds of rank r */
  inline IndexVector getLowerBounds(RankType r) const {
    assert(r >= 0 && r < size_);
    // get coords of r in cart comm
    std::vector<int> coords(dim_);
    IndexVector lowerBounds(dim_);
    getPartitionCoords(r, coords);

    for (DimType i = 0; i < dim_; ++i) {
      lowerBounds[i] = getDecomposition()[i][coords[i]];
    }
    return lowerBounds;
  }

  /** coordinates of this process' lower bounds */
  inline std::vector<real> getLowerBoundsCoords() const { return getLowerBoundsCoords(rank_); }

  /** coordinates of rank r's lower bounds */
  inline std::vector<real> getLowerBoundsCoords(RankType r) const {
    assert(r >= 0 && r < size_);
    const IndexVector& lbvi = getLowerBounds(r);
    IndexType lbli = getGlobalLinearIndex(lbvi);
    std::vector<real> coords(dim_);
    getCoordsGlobal(lbli, coords);
    return coords;
  }

  /** upper Bounds of this process */
  inline IndexVector getUpperBounds() const {
    // return getUpperBounds(rank_);
    return myPartitionsUpperBounds_;
  }

  /** upper bounds of rank r */
  inline IndexVector getUpperBounds(RankType r) const {
    assert(r >= 0 && r < size_);
    std::vector<int> coords(dim_);
    IndexVector upperBounds(dim_);
    getPartitionCoords(r, coords);

    for (DimType i = 0; i < dim_; ++i) {
      RankType n;
      std::vector<int> nc(coords);

      if (nc[i] < procs_[i] - 1) {
        // get rank of next neighbor in dim i
        nc[i] += 1;
        n = this->getRank(nc);
        upperBounds[i] = getLowerBounds(n)[i];
      } else {
        // no neighbor in dim i -> end of domain
        upperBounds[i] = nrPoints_[i];
      }
    }
    return upperBounds;
  }

  /** decomposition coords
   * contains same information as lowerBoundsCoords but in different
   * representation. here for each dim we have coordinate of the 1d lower
   * bound
   */
  inline std::vector<std::vector<real>> getDecompositionCoords() const {
    std::vector<std::vector<real>> decompositionCoords(dim_);
    for (size_t j = 0; j < dim_; ++j) decompositionCoords[j].resize(procs_[j]);
    assert(false && "this is pretty much untested, please add test before using this");

    for (RankType r = 0; r < size_; ++r) {
      // get coords of r in cart comm
      std::vector<int> coords(dim_);
      getPartitionCoords(r, coords);

      for (DimType i = 0; i < dim_; ++i) {
        decompositionCoords[i][coords[i]] = getLowerBoundsCoords(r)[i];
      }
    }
    return decompositionCoords;
  }

  /** coordinates of this process' upper bounds */
  inline std::vector<real> getUpperBoundsCoords() const { return getUpperBoundsCoords(rank_); }

  /** coordinates of rank r' upper bounds */
  inline std::vector<real> getUpperBoundsCoords(RankType r) const {
    assert(r >= 0 && r < size_);

    /* the upper bound index vector can correspond to coordinates outside of
     * the domain, thus we cannot use the getGlobalLinearIndex method here.
     */
    const IndexVector& ubvi = getUpperBounds(r);
    std::vector<real> coords(dim_);
    IndexType tmp_add = 0;

    for (DimType j = 0; j < dim_; ++j) {
      tmp_add = (hasBoundaryPoints_[j] == true) ? (0) : (1);
      coords[j] = static_cast<double>(ubvi[j] + tmp_add) * getGridSpacing()[j];
    }
    return coords;
  }

  /** Number of Grids in every dimension*/
  inline const IndexVector& getParallelization() const { return procs_; }

  /** MPI Rank */
  inline int getMpiRank() { return rank_; }

  /** MPI Size */
  inline int getMpiSize() { return size_; }

  /** position of a process in the grid of processes */
  inline void getPartitionCoords(RankType r, IndexVector& coords) const {
    assert(r >= 0 && r < size_);
    assert (!partitionCoords_.empty());
    coords.assign(partitionCoords_[r].begin(), partitionCoords_[r].end());
  }
  inline void getPartitionCoords(RankType r, std::vector<int>& coords) const {
    assert(r >= 0 && r < size_);
    assert (!partitionCoords_.empty());
    coords.assign(partitionCoords_[r].begin(), partitionCoords_[r].end());
  }

  /** position of the local process in the grid of processes */
  inline void getPartitionCoords(IndexVector& coords) const { getPartitionCoords(rank_, coords); }
  inline void getPartitionCoords(std::vector<int>& coords) const { getPartitionCoords(rank_, coords); }

  /** returns the 1d global index of the first point in the local domain
   *
   */
  inline IndexType getFirstGlobal1dIndex(DimType d) const {
    return getLowerBounds()[d];
  }

  IndexVector getFirstGlobalIndex() const {
    IndexVector firstGlobalIndex(dim_);
    for (DimType d = 0; d < dim_; ++d) {
      firstGlobalIndex[d] = getFirstGlobal1dIndex(d);
    }
    return firstGlobalIndex;
  }

  /** returns the 1d global index of the last point in the local domain
   *
   */
  inline IndexType getLastGlobal1dIndex(DimType d) const {
    return getUpperBounds()[d] - 1;
  }

  IndexVector getLastGlobalIndex() const {
    IndexVector lastGlobalIndex(dim_);
    for (DimType d = 0; d < dim_; ++d) {
      lastGlobalIndex[d] = getLastGlobal1dIndex(d);
    }
    return lastGlobalIndex;
  }

  // returns level of a global 1d index
  inline LevelType getLevel(DimType d, IndexType idx1d) const {
    // IndexVector givec(dim_, 0);
    // givec[d] = idx1d;
    // IndexType idx = getGlobalLinearIndex(givec);
    // LevelVector levels(dim_);
    // IndexVector tmp(dim_);
    // getGlobalLI(idx, levels, tmp);
    // return levels[d];

    // get level of idx1d by rightmost set bit
    // (e.g., all points on the finest level already have a 1 at the rightmost bit)
    if (idx1d == 0) {
      return 0;
    }
    const LevelType lmax = this->getLevels()[d];
    // auto set_bit = idx1d & ~(idx1d-1);
    // LevelType l = lmax - log2(set_bit);
    // builtin is fast and should work with gcc and clang
    // if it is not available, use the one above (at a slight performance hit)
    // or c++20's std::countr_zero
    LevelType l = lmax - __builtin_ctz(idx1d);
    return l;
  }

  // get 1d index of LeftPredecessor of a point
  // returns negative number if point has no left predecessor
  inline IndexType getLeftPredecessor(DimType d, IndexType idx1d) {
    LevelType l = getLevel(d, idx1d);

    // boundary points
    if (l == 0) return -1;

    LevelType ldiff = levels_[d] - l;
    IndexType lpidx = idx1d - static_cast<IndexType>(std::pow(2, ldiff));

    return lpidx;
  }

  inline IndexType getRightPredecessor(DimType d, IndexType idx1d) {
    LevelType l = getLevel(d, idx1d);

    // boundary points
    if (l == 0) return -1;

    LevelType ldiff = levels_[d] - l;
    IndexType rpidx = idx1d + static_cast<IndexType>(std::pow(2, ldiff));

    // check if outside of domain
    IndexType numElementsD = this->getGlobalSizes()[d];

    if (rpidx > numElementsD - 1) rpidx = -1;

    return rpidx;
  }

  // get coordinates of the partition which contains the point specified
  // by the global index vector
  inline void getPartitionCoords(IndexVector& globalAxisIndex, IndexVector& partitionCoords) {
    partitionCoords.resize(dim_);

    for (DimType d = 0; d < dim_; ++d) {
      partitionCoords[d] = -1;
      for (size_t i = 0; i < getDecomposition()[d].size(); ++i) {
        if (globalAxisIndex[d] >= getDecomposition()[d][i]) partitionCoords[d] = i;
      }

      // check whether the partition coordinates are valid
      assert(partitionCoords[d] > -1 && partitionCoords[d] < procs_[d]);
    }
  }

  inline RankType getRank(std::vector<int>& partitionCoordsInt) const {
    // check wheter the partition coords are valid
    assert(partitionCoordsInt.size() == dim_);

    for (DimType d = 0; d < dim_; ++d) assert(partitionCoordsInt[d] < procs_[d]);

    assert(!partitionCoords_.empty());
    auto dim = dim_;
    auto it = std::find_if(partitionCoords_.begin(), partitionCoords_.end(),
                           [&partitionCoordsInt, &dim](std::vector<int> pcoord) {
                             return (pcoord == partitionCoordsInt);
                           });
    bool found = (it != partitionCoords_.end());
    assert(found);

    RankType rank = static_cast<RankType>(it - partitionCoords_.begin());

    // RankType rank_old;
    // MPI_Cart_rank(communicator_, &partitionCoordsInt[0], &rank_old);

    // assert (rank_old == rank);

    return rank;
  }

  inline RankType getRank(IndexVector& partitionCoords) const {
    std::vector<int> partitionCoordsInt(partitionCoords.begin(), partitionCoords.end());

    return getRank(partitionCoordsInt);
  }

  inline void print(std::ostream& os) const {
    if (dim_ == 1)
      print1D(os);
    else if (dim_ == 2)
      print2D(os);
    else if (dim_ == 3)
      print3D(os);
    else
      assert(false && "Maximum dimension for printing is 2");
  }

  // return extents of local grid
  inline const IndexVector& getLocalSizes() const { return nrLocalPoints_; }

  // return extents of global grid
  inline const IndexVector& getGlobalSizes() const { return nrPoints_; }

  // return MPI DataType
  inline MPI_Datatype getMPIDatatype() const {
    return abstraction::getMPIDatatype(abstraction::getabstractionDataType<FG_ELEMENT>());
  }

  // gather fullgrid on rank r
  void gatherFullGrid(FullGrid<FG_ELEMENT>& fg, RankType root) {
    int size = this->getCommunicatorSize();
    int rank = this->getMpiRank();
    CommunicatorType comm = this->getCommunicator();

    // each rank: send dfg to root
    int dest = root;
    MPI_Request sendRequest;
    MPI_Isend(this->getData(), static_cast<int>(this->getNrLocalElements()), this->getMPIDatatype(),
              dest, 0, comm, &sendRequest);

    std::vector<MPI_Request> requests;
    std::vector<MPI_Datatype> subarrayTypes;

    // rank r: for each rank create subarray view on fg
    if (rank == root) {
      if (!fg.isGridCreated()) fg.createFullGrid();

      for (int r = 0; r < size; ++r) {
        IndexVector sizes(fg.getSizes().begin(), fg.getSizes().end());
        IndexVector subsizes = this->getUpperBounds(r) - this->getLowerBounds(r);
        IndexVector starts = this->getLowerBounds(r);

        // we store our data in c format, i.e. first dimension is the innermost
        // dimension. however, we access our data in fortran notation, with the
        // first index in indexvectors being the first dimension.
        // to comply with an ordering that mpi understands, we have to reverse
        // our index vectors
        // also we have to use int as datatype
        std::vector<int> csizes(sizes.rbegin(), sizes.rend());
        std::vector<int> csubsizes(subsizes.rbegin(), subsizes.rend());
        std::vector<int> cstarts(starts.rbegin(), starts.rend());
        // if (!reverseOrderingDFGPartitions) { // not sure why, but this produces the wrong results
        //   csizes.assign(sizes.begin(), sizes.end());
        //   csubsizes.assign(subsizes.begin(), subsizes.end());
        //   cstarts.assign(starts.begin(), starts.end());
        // }

        // create subarray view on data
        MPI_Datatype mysubarray;
        MPI_Type_create_subarray(static_cast<int>(this->getDimension()), &csizes[0], &csubsizes[0],
                                 &cstarts[0], MPI_ORDER_C, this->getMPIDatatype(), &mysubarray);
        MPI_Type_commit(&mysubarray);
        subarrayTypes.push_back(mysubarray);

        int src = r;
        MPI_Request req;
        MPI_Irecv(fg.getData(), 1, mysubarray, src, 0, this->getCommunicator(), &req);
        requests.push_back(req);
      }
    }

    // all
    MPI_Wait(&sendRequest, MPI_STATUS_IGNORE);

    if (rank == root) {
      MPI_Waitall(static_cast<int>(requests.size()), &requests[0], MPI_STATUSES_IGNORE);
    }

    // free subarrays
    for (size_t i = 0; i < subarrayTypes.size(); ++i) MPI_Type_free(&subarrayTypes[i]);
  }

  // /* extract subspaces from dfg
  //  * note that this will double the memory demand!
  //  */
  // void fillSubspaces() {
  //   // resize subspace if necessary
  //   for (size_t j = 0; j < subspaces_.size(); ++j) {
  //     SubspaceDFG<FG_ELEMENT>& subsp = subspaces_[j];

  //     // compute expected data size
  //     IndexVector lsize = subsp.upperBounds_[rank_] - subsp.lowerBounds_[rank_];
  //     IndexType dsize = 1;

  //     for (auto l : lsize) dsize *= l;

  //     if (IndexType(subsp.data_.size()) != dsize) subsp.data_.resize(dsize);
  //   }

  //   // create iterator for each subspace of sgrid and set to begin
  //   typedef typename std::vector<FG_ELEMENT>::iterator SubspaceIterator;
  //   typename std::vector<SubspaceIterator> it_sub(subspaces_.size());

  //   for (size_t i = 0; i < it_sub.size(); ++i) it_sub[i] = subspaces_[i].data_.begin();

  //   for (size_t i = 0; i < fullgridVector_.size(); ++i) {
  //     // get subspace id
  //     size_t subI(assigmentList_[i]);

  //     assert(it_sub[subI] != subspaces_[subI].data_.end());

  //     // copy subspace value back to dfg
  //     *it_sub[subI] = fullgridVector_[i];

  //     ++it_sub[subI];
  //   }

  //   subspacesFilled_ = true;
  // }

  void calcLocalSizeRecursive(DimType d, const LevelVector& lvec, IndexVector& ivec, size_t& localSize) {
    IndexVector oneDIndices;
    get1dIndicesLocal(d, lvec, oneDIndices);

    // stupid additive recursion, could be multiplicative (since hyper-rectangular)
    for (IndexType idx : oneDIndices) {
      ivec[d] = idx;

      if (d > 0)
        calcLocalSizeRecursive(d - 1, lvec, ivec, localSize);
      else {
        ++localSize;
      }
    }
  }

  /**
   * @brief Get the number of points of the subspace on this partition
   *
   * @param l level of hierarchical subspace
   * @return size_t the number of points on this partition
   */
  size_t getLocalSizeOfSubspaceOnThisPartition(LevelVector l) {
    IndexVector ivec(dim_);
    size_t localSize = 0;

    calcLocalSizeRecursive(dim_ - 1, l, ivec, localSize);
    return localSize;
  }

  void getFGPointsOfSubspaceRecursive(DimType d, const LevelVector& lvec, IndexVector& ivec, std::vector<IndexType>& subspaceIndices) {
    IndexVector oneDIndices;
    get1dIndicesLocal(d, lvec, oneDIndices);

    for (IndexType idx : oneDIndices) {
      ivec[d] = idx;
      if (d > 0)
        getFGPointsOfSubspaceRecursive(d - 1, lvec, ivec, subspaceIndices);
      else {
        IndexType j = getLocalLinearIndex(ivec);
        subspaceIndices.emplace_back(j);
      }
    }
  }

  /**
   * @brief Get the indices of the points of the subspace on this partition
   *
   * @param l level of hierarchical subspace
   * @return size_t the indices of points on this partition
   */
  std::vector<IndexType> getFGPointsOfSubspace(LevelVector l) {
    std::vector<IndexType> subspaceIndices;
    IndexVector ivec(dim_);

    getFGPointsOfSubspaceRecursive(dim_ - 1, l, ivec, subspaceIndices);
    return subspaceIndices;
  }

  /**
   * @brief "registers" the DistributedSparseGridUniform with this DistributedFullGrid:
   *        sets the dsg_ member,
   *        and sets the dsg's subspaceSizes where they are not yet set:
   *        If the subspaces in the dsg have zero size, all subspaces
   *        of the dsg that the dfg and dsg have in common are resized. The
   *        size of a subspace in the dsg is chosen according to the corresponding
   *        subspace size in the dfg.
   *
   * @param dsg the DSG to register
   */
  void registerUniformSG(DistributedSparseGridUniform<FG_ELEMENT>& dsg) {
    dsg_ = &dsg;
    assert(dsg_->getDim() == dim_);

    subspaceIndexToFGIndices_.clear();
    // subspaceIndexToFGIndices_.reserve(std::accumulate(levels_.begin(), levels_.end(), 1,
    //                                    std::multiplies<LevelType>()));

    // a sparse grid that knows all the hierarchical subspaces
    //    contained in this full grid
    SGrid<bool> sg(dim_, levels_, levels_, hasBoundaryPoints_);
    // resize all common subspaces in dsg, if necessary
    for (size_t subspaceID = 0; subspaceID < sg.getSize(); ++subspaceID) {
      auto level = sg.getLevelVector(subspaceID);

      if (dsg_->isContained(level)) {
        auto index = dsg_->getIndex(level);
        // auto lsize = getLocalSizeOfSubspaceOnThisPartition(level);
        auto FGIndices = getFGPointsOfSubspace(level);
        auto lsize = FGIndices.size();

        size_t subSgDataSize = dsg_->getDataSize(index);
        // resize DSG subspace if it has zero size
        if (subSgDataSize == 0) {
          dsg_->setDataSize(index, lsize);
          // subSgData.resize(lsize);
        } else {
          ASSERT(subSgDataSize == lsize,
                "subSgDataSize: " << subSgDataSize << ", lsize: "
                                      << lsize << std::endl);
        }
        subspaceIndexToFGIndices_.push_back(std::make_pair(index, FGIndices));
      }
    }
    // if localFGIndexToLocalSGPointerList_ was already initialized,
    // it is invalidated now
    localFGIndexToLocalSGPointerList_.clear();
  }


  void setLocalFGPointersRecursive(DimType d, const LevelVector& lvec, IndexVector& ivec, FG_ELEMENT*& sgDataPointer) {
    IndexVector oneDIndices;
    get1dIndicesLocal(d, lvec, oneDIndices);

    for (IndexType idx : oneDIndices) {
      ivec[d] = idx;

      if (d > 0)
        setLocalFGPointersRecursive(d - 1, lvec, ivec, sgDataPointer);
      else {
        IndexType j = getLocalLinearIndex(ivec);
        localFGIndexToLocalSGPointerList_[j] = sgDataPointer;
        ++sgDataPointer;
      }
    }
  }


  /**
   * @brief set localFGIndexToLocalSGPointerList_ based on the current sizes of
   *        dsg_
   *       (dsg_'s sizes may have changed as other full grids were registered)
   *        and return it
   */
  std::vector<FG_ELEMENT*>& getLocalFGIndexToLocalSGPointerList() {
    if (localFGIndexToLocalSGPointerList_.size() == 0) {
      // make sure that using pointers instead of indices is sensible
      assert(sizeof(FG_ELEMENT*) <= sizeof(IndexType));

      localFGIndexToLocalSGPointerList_.resize(nrLocalElements_, nullptr);

      // a sparse grid that knows all the hierarchical subspaces
      //    contained in this full grid
      SGrid<bool> sg(dim_, levels_, levels_, hasBoundaryPoints_);
      // resize all common subspaces in dsg, if necessary
      for (size_t subspaceID = 0; subspaceID < sg.getSize(); ++subspaceID) {
        auto level = sg.getLevelVector(subspaceID);
        if (dsg_->isContained(level)) {
          FG_ELEMENT* dpointer = dsg_->getData(level);
          IndexVector ivec(dim_);
          setLocalFGPointersRecursive(dim_ - 1, level, ivec, dpointer);
          assert((dpointer-dsg_->getData(level)) == dsg_->getDataSize(level));
        }
      }
    }
    return localFGIndexToLocalSGPointerList_;
  }

  /**
   * @brief adds the (hopefully) hierarchical coefficients from fullgridVector_
   *        to the dsg's data structure, multiplied by coeff
   *
   * @param dsg the DSG to add to
   * @param coeff the coefficient that gets multiplied to all entries
   */
  void addToUniformSG(DistributedSparseGridUniform<FG_ELEMENT>& dsg, real coeff) {
    // test if dsg has already been registered
    // assert (&dsg == dsg_); // if using getLocalFGIndexToLocalSGPointerList,
    // this needs to be asserted.
    if (dsg_ != &dsg) {
      this->registerUniformSG(dsg);
    }

    // auto indexAssignment = getLocalFGIndexToLocalSGPointerList();
    // #ifdef DEBUG_OUTPUT
    //   MASTER_EXCLUSIVE_SECTION { std::cout << "is this where " << indexAssignment[0] << "?\n"; }
    // #endif

    bool anythingWasAdded = false;

    // // loop over all grid points
    // for (size_t i = 0; i < nrLocalElements_; ++i) {
    //   // add grid point to subspace, mul with coeff
    //   if (indexAssignment[i] != nullptr) {
    //     *(indexAssignment[i]) += coeff * fullgridVector_[i];
    //     anythingWasAdded = true;
    //   }
    // }

    assert(!subspaceIndexToFGIndices_.empty());

    // loop over all subspaces
    for (const auto& sToF : subspaceIndexToFGIndices_) {
      const auto sIndex = sToF.first;
      auto sPointer = dsg_->getData(sIndex);
      for (const auto& fIndex : sToF.second) {
        *sPointer += coeff * fullgridVector_[fIndex];
        ++sPointer;
        anythingWasAdded = true;
      }
    }
    // make sure that anything was added -- I can only think of weird setups
    // where that would not be the case
    assert(anythingWasAdded);
  }

  /**
   * @brief extracts the (hopefully) hierarchical coefficients from dsg_
   *        to the full grid's data structure
   *
   * @param dsg the DSG to extract from
   */
  void extractFromUniformSG(DistributedSparseGridUniform<FG_ELEMENT>& dsg) {
    // test if dsg has already been registered
    if (dsg_ != &dsg) {
      this->registerUniformSG(dsg);
    }

    // auto indexAssignment = getLocalFGIndexToLocalSGPointerList();

    // // loop over all grid points
    // for (size_t i = 0; i < nrLocalElements_; ++i) {
    //   // check if fg element is contained in sg
    //   // (this might not be the case when a reduced sg is used)
    //   if (indexAssignment[i] != nullptr) {
    //     fullgridVector_[i] = *(indexAssignment[i]);
    //   }
    // }

    // loop over all subspaces (-> somewhat linear access in the sg)
    for (const auto& sToF : subspaceIndexToFGIndices_) {
      const auto sIndex = sToF.first;
      auto sPointer = dsg_->getData(sIndex);
      for (const auto& fIndex : sToF.second) {
        fullgridVector_[fIndex] = *sPointer;
        ++sPointer;
      }
    }
  }

  // void writeBackSubspaces() {
  //   assert(subspacesFilled_ && "subspaces have not been created");

  //   // create iterator for each subspace of sgrid and set to begin
  //   typedef typename std::vector<FG_ELEMENT>::iterator SubspaceIterator;
  //   typename std::vector<SubspaceIterator> it_sub(subspaces_.size());

  //   for (size_t i = 0; i < it_sub.size(); ++i) it_sub[i] = subspaces_[i].data_.begin();

  //   for (size_t i = 0; i < fullgridVector_.size(); ++i) {
  //     // get subspace id corresponding to lvec
  //     size_t subI(assigmentList_[i]);

  //     assert(it_sub[subI] != subspaces_[subI].data_.end());

  //     // copy subspace value back to dfg
  //     fullgridVector_[i] = *it_sub[subI];

  //     ++it_sub[subI];
  //   }
  // }

  // void clearSubspaces() {
  //   for (auto subsp : subspaces_) subsp.data_.resize(0);

  //   subspacesFilled_ = false;
  // }

  // void gatherSubspace(const LevelVector& l, RankType dst, std::vector<FG_ELEMENT>& buf) {
  //   assert(subspacesFilled_ && "subspaces have not been filled");

  //   // get subspace corresponding to l
  //   size_t subI = this->getSubspaceIndex(l);
  //   SubspaceDFG<FG_ELEMENT>& subsp = subspaces_[subI];

  //   std::vector<MPI_Request> requests;
  //   std::vector<MPI_Datatype> subarrayTypes;

  //   // rank r: for each rank create subarray view on fg
  //   if (rank_ == dst) {
  //     // resize buffer
  //     IndexType bsize = 1;

  //     for (auto s : subsp.sizes_) bsize *= s;

  //     buf.resize(bsize);

  //     for (int r = 0; r < size_; ++r) {
  //       MPI_Datatype mysubarray = subsp.subarrayTypes_[r];

  //       if (mysubarray == MPI_DATATYPE_NULL) continue;

  //       int src = r;
  //       MPI_Request req;
  //       MPI_Irecv(buf.data(), 1, mysubarray, src, 0, this->getCommunicator(), &req);
  //       requests.push_back(req);
  //     }
  //   }

  //   // each rank: send subspace to dst
  //   // important: skip this if send size = 0
  //   if (subsp.data_.size() > 0) {
  //     MPI_Send(subsp.data_.data(), int(subsp.data_.size()), this->getMPIDatatype(), dst, 0,
  //              this->getCommunicator());
  //   }

  //   if (rank_ == dst) {
  //     MPI_Waitall(static_cast<int>(requests.size()), &requests[0], MPI_STATUSES_IGNORE);
  //   }

  //   // free subarrays; only dst has a non-zero container here
  //   for (size_t i = 0; i < subarrayTypes.size(); ++i) MPI_Type_free(&subarrayTypes[i]);
  // }

  // void gatherSubspaceBlock(const LevelVector& l, RankType dst,
  //                          std::vector<std::vector<FG_ELEMENT> >& senddata,
  //                          std::vector<std::vector<int> >& sendsizes,
  //                          std::vector<std::vector<int> >& sendsubspaces,
  //                          std::vector<std::vector<int> >& recvsizes,
  //                          std::vector<std::vector<int> >& recvsubspaces) {
  //   assert(subspacesFilled_ && "subspaces have not been filled");

  //   // get subspace corresponding to l
  //   size_t subI = this->getSubspaceIndex(l);
  //   SubspaceDFG<FG_ELEMENT>& subsp = subspaces_[subI];

  //   // rank r: for each rank create subarray view on fg
  //   if (rank_ == dst) {
  //     for (int r = 0; r < size_; ++r) {
  //       IndexVector sizes(subsp.sizes_);
  //       IndexVector subsizes = subsp.upperBounds_[r] - subsp.lowerBounds_[r];
  //       IndexVector starts = subsp.lowerBounds_[r];

  //       // important: skip r if subsize = 0
  //       IndexType ssize = 1;

  //       for (auto s : subsizes) ssize *= s;

  //       if (ssize == 0) continue;

  //       recvsizes[r].push_back(ssize);
  //       recvsubspaces[r].push_back(subI);
  //     }
  //   }

  //   // send subspace copy subspace data into right buffer
  //   // skip this if send size = 0
  //   if (subsp.data_.size() > 0) {
  //     senddata[dst].insert(senddata[dst].end(), subsp.data_.begin(), subsp.data_.end());
  //     sendsizes[dst].push_back(subsp.data_.size());
  //     sendsubspaces[dst].push_back(subI);
  //   }
  // }

  // void gatherSubspaceRed(const LevelVector& l, RankType dst, std::vector<FG_ELEMENT>& buf) {
  //   assert(subspacesFilled_ && "subspaces have not been filled");

  //   // get subspace corresponding to l
  //   size_t subI = this->getSubspaceIndex(l);
  //   SubspaceDFG<FG_ELEMENT>& subsp = subspaces_[subI];

  //   assert(buf.size() >= subsp.targetSize_);

  //   // set interesting part of buf to zero
  //   for (size_t i = 0; i < subsp.targetSize_; ++i) buf[i] = FG_ELEMENT(0);

  //   // todo: copy data into right position in buffer

  //   // perform reduce with target dst
  //   if (rank_ == dst) {
  //     MPI_Reduce(MPI_IN_PLACE, buf.data(), int(subsp.targetSize_), this->getMPIDatatype(), MPI_SUM,
  //                dst, this->getCommunicator());
  //   } else {
  //     MPI_Reduce(buf.data(), buf.data(), int(subsp.targetSize_), this->getMPIDatatype(), MPI_SUM,
  //                dst, this->getCommunicator());
  //   }
  // }

  // void gatherSubspaceNB(const LevelVector& l, RankType dst, std::vector<FG_ELEMENT>& buf,
  //                       std::vector<MPI_Request>& requests) {
  //   assert(subspacesFilled_ && "subspaces have not been filled");

  //   // get subspace corresponding to l
  //   size_t subI = this->getSubspaceIndex(l);
  //   SubspaceDFG<FG_ELEMENT>& subsp = subspaces_[subI];
  //   int tag = int(subI);

  //   // each rank: send subspace to dst
  //   // important: skip this if send size = 0
  //   if (subsp.data_.size() > 0) {
  //     MPI_Request sendRequest;
  //     MPI_Isend(subsp.data_.data(), int(subsp.data_.size()), this->getMPIDatatype(), dst, tag,
  //               this->getCommunicator(), &sendRequest);
  //     requests.push_back(sendRequest);
  //   }

  //   // rank r: for each rank create subarray view on fg
  //   if (rank_ == dst) {
  //     // resize buffer
  //     IndexType bsize = 1;

  //     for (auto s : subsp.sizes_) bsize *= s;

  //     buf.resize(bsize);

  //     for (int r = 0; r < size_; ++r) {
  //       MPI_Datatype mysubarray = subsp.subarrayTypes_[r];

  //       if (mysubarray == MPI_DATATYPE_NULL) continue;

  //       int src = r;
  //       MPI_Request req;
  //       MPI_Irecv(buf.data(), 1, mysubarray, src, tag, this->getCommunicator(), &req);
  //       requests.push_back(req);
  //     }
  //   }

  //   // std::vector< MPI_Datatype > subarrayTypes;

  //   // free subarrays; only dst has a non-zero container here
  //   // todo: free on destruct
  //   // for( size_t i = 0; i < subarrayTypes.size(); ++i )
  //   //  MPI_Type_free( &subarrayTypes[i] );
  // }

  // /* takes the data of subspace l contained in buf and distributes
  //  * it over all processes of dfg
  //  */
  // void scatterSubspace(const LevelVector& l, RankType src, const std::vector<FG_ELEMENT>& buf) {
  //   assert(subspacesFilled_ && "subspaces have not been created");

  //   // get subspace corresponding to l
  //   size_t subI = this->getSubspaceIndex(l);
  //   SubspaceDFG<FG_ELEMENT>& subsp = subspaces_[subI];

  //   // each rank: recv his subspace data from src
  //   // important: skip this if send size = 0
  //   MPI_Request recvRequest;
  //   bool recvd(false);

  //   if (subsp.data_.size() > 0) {
  //     MPI_Irecv(subsp.data_.data(), int(subsp.data_.size()), this->getMPIDatatype(), src, 0,
  //               this->getCommunicator(), &recvRequest);
  //     recvd = true;
  //   }

  //   std::vector<MPI_Request> requests;
  //   std::vector<MPI_Datatype> subarrayTypes;

  //   // rank r: for each rank create subarray view on fg
  //   if (rank_ == src) {
  //     // check buffer size
  //     IndexType bsize = 1;

  //     for (auto s : subsp.sizes_) bsize *= s;

  //     assert(IndexType(buf.size()) == bsize);

  //     for (int r = 0; r < size_; ++r) {
  //       IndexVector sizes(subsp.sizes_);
  //       IndexVector subsizes = subsp.upperBounds_[r] - subsp.lowerBounds_[r];
  //       IndexVector starts = subsp.lowerBounds_[r];

  //       // important: skip r if subsize = 0
  //       IndexType ssize = 1;

  //       for (auto s : subsizes) ssize *= s;

  //       if (ssize == 0) continue;

  //       // we store our data in c format, i.e. first dimension is the innermost
  //       // dimension. however, we access our data in fortran notation, with the
  //       // first index in indexvectors being the first dimension.
  //       // to comply with an ordering that mpi understands, we have to reverse
  //       // our index vectors
  //       // also we have to use int as datatype
  //       std::vector<int> csizes(sizes.rbegin(), sizes.rend());
  //       std::vector<int> csubsizes(subsizes.rbegin(), subsizes.rend());
  //       std::vector<int> cstarts(starts.rbegin(), starts.rend());
  //       if (!reverseOrderingDFGPartitions) {
  //         csizes.assign(sizes.begin(), sizes.end());
  //         csubsizes.assign(subsizes.begin(), subsizes.end());
  //         cstarts.assign(starts.begin(), starts.end());
  //       }

  //       // create subarray view on data
  //       MPI_Datatype mysubarray;
  //       MPI_Type_create_subarray(int(this->getDimension()), &csizes[0], &csubsizes[0], &cstarts[0],
  //                                MPI_ORDER_C, this->getMPIDatatype(), &mysubarray);
  //       MPI_Type_commit(&mysubarray);
  //       subarrayTypes.push_back(mysubarray);

  //       int src = r;
  //       MPI_Request req;
  //       MPI_Isend(buf.data(), 1, mysubarray, src, 0, this->getCommunicator(), &req);
  //       requests.push_back(req);
  //     }
  //   }

  //   // all
  //   if (recvd) MPI_Wait(&recvRequest, MPI_STATUS_IGNORE);

  //   if (rank_ == src) {
  //     MPI_Waitall(static_cast<int>(requests.size()), &requests[0], MPI_STATUSES_IGNORE);
  //   }

  //   // free subarrays; only dst has a non-zero container here
  //   for (size_t i = 0; i < subarrayTypes.size(); ++i) MPI_Type_free(&subarrayTypes[i]);
  // }

//   size_t getSubspaceIndex(const LevelVector& l) const {
//     /* boundary points can have level 0. however, we do not store them
//      * in additional subspaces, but in the level 1 subspaces. thus,
//      * we have to correct the levelvector.
//      */
//     LevelVector myL(l);

//     for (size_t k = 0; k < myL.size(); ++k)
//       if (myL[k] == 0) myL[k] = 1;

//     // get index of l
//     bool found = false;
//     size_t i;

//     for (i = 0; i < subspaces_.size(); ++i) {
//       if (subspaces_[i].level_ == myL) {
//         found = true;
//         break;
//       }
//     }

// #ifdef DEBUG_OUTPUT
//     if (!found) std::cout << "subspace " << myL << " not included" << std::endl;
// #endif

//     assert(found);

//     return i;
//   }

  // inline const LevelVector& getSubspaceLevelVector(size_t i) const {
  //   assert(i < subspaces_.size());
  //   return subspaces_[i].level_;
  // }

  // inline size_t getSubspaceIndex(LevelVector l) {
  //   for (size_t i = 0; i < subspaces_.size(); ++i) {
  //     if (subspaces_[i].level_ == l) return i;
  //   }
  //   return subspaces_.size();
  // }

  // void getSubspacesLevelVectors(std::vector<LevelVector>& lvecs) {
  //   for (auto subsp : subspaces_) lvecs.push_back(subsp.level_);
  // }

  inline size_t getNumSubspaces() const {
    SGrid<bool> sg(dim_, levels_, levels_, hasBoundaryPoints_);
    return sg.getSize();
  //   return subspaces_.size();
  }

  // inline std::vector<FG_ELEMENT>& getSubspaceData(size_t i) {
  //   assert(i < subspaces_.size());
  //   return subspaces_[i].data_;
  // }

  // inline std::vector<FG_ELEMENT>& getSubspaceData(LevelVector l) {
  //   size_t i = this->getSubspaceIndex(l);
  //   assert(i != subspaces_.size());
  //   return subspaces_[i].data_;
  // }

  /** get the pointwise average of the ("small") lp Norm:
   * (1/numPoints * sum_i(abs(x_i)^p))^1/p
   * (almost Lp-Norm, except for the weight of boundary values)
   *
   * TODO: change this to add boundary values with reduced value,
   *        divide by more in case of no boundary -> big Lp norm on [0; 1]
   */
  real getLpNorm(int p) const {
    assert(p >= 0);
    // special case maximum norm
    MPI_Datatype dtype =
      abstraction::getMPIDatatype(abstraction::getabstractionDataType<real>());
    if (p == 0) {
      auto& data = getElementVector();
      real max = 0.0;

      for (size_t i = 0; i < data.size(); ++i) {
        if (std::abs(data[i]) > max) max = std::abs(data[i]);
      }

      real globalMax(-1);
      MPI_Allreduce(&max, &globalMax, 1, dtype, MPI_MAX, getCommunicator());

      return globalMax;
    } else {
      real p_f = static_cast<real>(p);
      auto& data = getElementVector();
      real res = 0.0;

      for (size_t i = 0; i < data.size(); ++i) {
        real abs = std::abs(data[i]);
        res += std::pow(abs, p_f);
      }
      // res /= data.size();
      real globalRes(0.);
      MPI_Allreduce(&res, &globalRes, 1, dtype, MPI_SUM, getCommunicator());
      return std::pow(globalRes, 1.0 / p_f);
    }
  }

  // normalize with specific norm
  inline real normalizelp(int p) {
    real norm = getLpNorm(p);

    std::vector<FG_ELEMENT>& data = getElementVector();
    for (size_t i = 0; i < data.size(); ++i) data[i] = data[i] / norm;

    return norm;
  }

  // multiply with a constant
  inline void mul(real c) {
    std::vector<FG_ELEMENT>& data = getElementVector();
    for (size_t i = 0; i < data.size(); ++i) data[i] *= c;
  }

  // write data to file using MPI-IO
  void writePlotFile(const char* filename) const {
    auto dim = getDimension();

    // create subarray data type
    IndexVector sizes = getGlobalSizes();
    IndexVector subsizes = getUpperBounds() - getLowerBounds();
    IndexVector starts = getLowerBounds();

    // we store our data in c format, i.e. first dimension is the innermost
    // dimension. however, we access our data in fortran notation, with the
    // first index in indexvectors being the first dimension.
    // to comply with an ordering that mpi understands, we have to reverse
    // our index vectors
    // also we have to use int as datatype
    std::vector<int> csizes(sizes.rbegin(), sizes.rend());
    std::vector<int> csubsizes(subsizes.rbegin(), subsizes.rend());
    std::vector<int> cstarts(starts.rbegin(), starts.rend());

    // create subarray view on data
    MPI_Datatype mysubarray;
    MPI_Type_create_subarray(static_cast<int>(getDimension()), &csizes[0], &csubsizes[0],
                             &cstarts[0], MPI_ORDER_C, getMPIDatatype(), &mysubarray);
    MPI_Type_commit(&mysubarray);

    // open file
    MPI_File fh;
    MPI_File_open(getCommunicator(), filename, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL,
                  &fh);

    MPI_Offset offset = 0;

    // .raw files can be read by paraview, in that case write the header separately
    std::string fn = filename;
    if (fn.find(".raw") != std::string::npos){
      // rank 0 write human-readable header
      if (rank_ == 0) {
        auto headername = fn + "_header";
        std::ofstream ofs(headername);

        // first line: dimension
        ofs << "dimensionality " << getDimension() << std::endl;

        // grid points per dimension in the order
        // x_0, x_1, ... , x_d
        ofs << "Extents ";
        for (auto s : sizes){ ofs << s << " ";}
        ofs << std::endl;

        // data type
        ofs << "Data type size " << sizeof(FG_ELEMENT) << std::endl;
        //TODO (pollinta) write endianness and data spacing
      }
    }else{
      // rank 0 write dim and resolution (and data format?)
      if (rank_ == 0) {
        MPI_File_write(fh, &dim, 1, MPI_INT, MPI_STATUS_IGNORE);

        std::vector<int> res(sizes.begin(), sizes.end());
        MPI_File_write(fh, &res[0], 6, MPI_INT, MPI_STATUS_IGNORE);
      }

      // set file view to right offset (in bytes)
      // 1*int + dim*int = (1+dim)*sizeof(int)
      offset = (1 + dim) * sizeof(int);
    }

    MPI_File_set_view(fh, offset, getMPIDatatype(), mysubarray, "native", MPI_INFO_NULL);

    // write subarray
     MPI_File_write_all(fh, getData(), static_cast<int>(getNrLocalElements()), getMPIDatatype(), MPI_STATUS_IGNORE);
    // close file
    MPI_File_close(&fh);
    MPI_Type_free(&mysubarray);
  }

  // write data to legacy-type VTK file using MPI-IO
  void writePlotFileVTK(const char* filename) const {
    auto dim = getDimension();
    assert(dim < 4);  // vtk supports only up to 3D

    // create subarray data type
    IndexVector sizes = getGlobalSizes();
    IndexVector subsizes = getUpperBounds() - getLowerBounds();
    IndexVector starts = getLowerBounds();

    // we store our data in c format, i.e. first dimension is the innermost
    // dimension. however, we access our data in fortran notation, with the
    // first index in indexvectors being the first dimension.
    // to comply with an ordering that mpi understands, we have to reverse
    // our index vectors
    // also we have to use int as datatype
    std::vector<int> csizes(sizes.rbegin(), sizes.rend());
    std::vector<int> csubsizes(subsizes.rbegin(), subsizes.rend());
    std::vector<int> cstarts(starts.rbegin(), starts.rend());

    // create subarray view on data
    MPI_Datatype mysubarray;
    MPI_Type_create_subarray(static_cast<int>(getDimension()), &csizes[0], &csubsizes[0],
                             &cstarts[0], MPI_ORDER_C, getMPIDatatype(), &mysubarray);
    MPI_Type_commit(&mysubarray);

    // open file
    MPI_File fh;
    MPI_File_open(getCommunicator(), filename, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL,
                  &fh);

    std::stringstream vtk_header;
    // cf https://lorensen.github.io/VTKExamples/site/VTKFileFormats/ => structured points
    vtk_header << "# vtk DataFile Version 2.0.\n"
               << "This file contains the combination solution evaluated on a full grid\n"
               << "BINARY\n"
               << "DATASET STRUCTURED_POINTS\n";
    if (dim == 3) { //TODO change for non-boundary grids using getGridSpacing
      vtk_header << "DIMENSIONS " << sizes[0]  << " " << sizes[1]  << " " << sizes[2] << "\n"
                 << "ORIGIN 0 0 0\n"
                 << "SPACING " << 1. / static_cast<double>(sizes[0]-1)  << " " << 1. / static_cast<double>(sizes[1]-1)  << " " << 1. / static_cast<double>(sizes[2]-1) << "\n";
    } else if (dim == 2) {
      vtk_header << "DIMENSIONS " << sizes[0]  << " " << sizes[1]  << " 1\n"
                 << "ORIGIN 0 0 0\n"
                 << "SPACING " << 1. / static_cast<double>(sizes[0]-1)  << " " << 1. / static_cast<double>(sizes[1]-1)  << " 1\n";
    } else if (dim == 1) {
      vtk_header << "DIMENSIONS " << sizes[0]  << " 1 1\n"
                 << "ORIGIN 0 0 0\n"
                 << "SPACING " << 1. / static_cast<double>(sizes[0]-1)  << " 1 1\n";
    } else {
      assert(false);
    }
    vtk_header << "POINT_DATA "
               << std::accumulate(sizes.begin(), sizes.end(), 1, std::multiplies<IndexType>())
               << "\n"
               << "SCALARS quantity double 1\n"
               << "LOOKUP_TABLE default\n";
    //TODO set the right data type from combidatatype, for now double by default
    bool rightDataType = std::is_same<CombiDataType, double>::value;
    assert(rightDataType);
    auto header_string = vtk_header.str();
    auto header_size = header_string.size();

    // rank 0 write header
    if (rank_ == 0) {
      MPI_File_write(fh, header_string.c_str(), static_cast<int>(header_size), MPI_CHAR, MPI_STATUS_IGNORE);
    }

    // set file view to right offset (in bytes)
    MPI_Offset offset = header_size * sizeof(char);
    // external32 not supported in OpenMPI < 5. -> writes "native" endianness
    // might work with MPICH
    MPI_File_set_view(fh, offset, getMPIDatatype(), mysubarray, "external32", MPI_INFO_NULL);

    // write subarray
    MPI_File_write_all(fh, getData(), static_cast<int>(getNrLocalElements()), getMPIDatatype(), MPI_STATUS_IGNORE);

    // close file
    MPI_File_close(&fh);
    MPI_Type_free(&mysubarray);
  }

  const std::vector<IndexVector>& getDecomposition() const { return decomposition_; }

  std::vector<MPI_Datatype> getUpwardSubarrays() {
    // initialize upwardSubarrays_ only once
    if (upwardSubarrays_.size() == 0){
      upwardSubarrays_.resize(this->getDimension());
      for (DimType d = 0; d < this->getDimension(); ++d) {
        // do index calculations
        // set lower bounds of subarray
        IndexVector subarrayLowerBounds = this->getLowerBounds();
        IndexVector subarrayUpperBounds = this->getUpperBounds();
        subarrayLowerBounds[d] += this->getLocalSizes()[d] - 1;

        auto subarrayExtents = subarrayUpperBounds - subarrayLowerBounds;
        assert(subarrayExtents[d] == 1);
        auto subarrayStarts = subarrayLowerBounds - this->getLowerBounds();

        // create MPI datatype
        // also, the data dimensions are reversed
        std::vector<int> sizes(this->getLocalSizes().rbegin(), this->getLocalSizes().rend());
        std::vector<int> subsizes(subarrayExtents.rbegin(), subarrayExtents.rend());
        // the starts are local indices
        std::vector<int> starts(subarrayStarts.rbegin(), subarrayStarts.rend());

        // create subarray view on data //todo do this only once per dimension
        MPI_Datatype mysubarray;
        MPI_Type_create_subarray(static_cast<int>(this->getDimension()), sizes.data(),
                                subsizes.data(), starts.data(), MPI_ORDER_C, this->getMPIDatatype(),
                                &mysubarray);
        MPI_Type_commit(&mysubarray);
        upwardSubarrays_[d] = mysubarray;
      }
    }
    return upwardSubarrays_;
  }

  std::vector<MPI_Datatype> getDownwardSubarrays() {
    // initialize downwardSubarrays_ only once
    if (downwardSubarrays_.size() == 0){
      downwardSubarrays_.resize(this->getDimension());
      for (DimType d = 0; d < this->getDimension(); ++d) {
        // do index calculations
        // set upper bounds of subarray
        IndexVector subarrayLowerBounds = this->getLowerBounds();
        IndexVector subarrayUpperBounds = this->getUpperBounds();
        subarrayUpperBounds[d] -= this->getLocalSizes()[d] - 1;

        auto subarrayExtents = subarrayUpperBounds - subarrayLowerBounds;
        assert(subarrayExtents[d] == 1);
        auto subarrayStarts = subarrayLowerBounds - this->getLowerBounds();

        // create MPI datatype
        // also, the data dimensions are reversed
        std::vector<int> sizes(this->getLocalSizes().rbegin(), this->getLocalSizes().rend());
        std::vector<int> subsizes(subarrayExtents.rbegin(), subarrayExtents.rend());
        // the starts are local indices
        std::vector<int> starts(subarrayStarts.rbegin(), subarrayStarts.rend());

        // create subarray view on data
        MPI_Datatype mysubarray;
        MPI_Type_create_subarray(static_cast<int>(this->getDimension()), sizes.data(),
                                subsizes.data(), starts.data(), MPI_ORDER_C, this->getMPIDatatype(),
                                &mysubarray);
        MPI_Type_commit(&mysubarray);
        downwardSubarrays_[d] = mysubarray;
      }
    }
    return downwardSubarrays_;
  }
  /**
   * @brief get the ranks of the highest and lowest "neighbor" rank in dimension d
   *    only sets highest and lowest if they actually are my neighbors
   */
  void getHighestAndLowestNeighbor(DimType d, int& highest, int& lowest) {
    //TODO this is not going to work with periodic cartesian communicator
    lowest = MPI_PROC_NULL;
    highest = MPI_PROC_NULL;

    auto d_reverse = this->getDimension() - d - 1;
    if (!reverseOrderingDFGPartitions) {
      d_reverse = d;
    }

    MPI_Cart_shift( this->getCommunicator(), d_reverse, getParallelization()[d] -1 , &lowest, &highest );

    // assert only boundaries have those neighbors (remove in case of periodicity)
    // this assumes no periodicity!
    if(! this->getLowerBounds()[d] == 0){
      assert(highest < 0);
    }
    if(! this->getUpperBounds()[d] == this->getGlobalSizes()[d]){
      assert(lowest < 0);
    }
  }

  void writeUpperBoundaryToLowerBoundary(DimType d) {
    assert(hasBoundaryPoints_[d] == true);

    // create MPI datatypes
    auto downSubarrays = getDownwardSubarrays();
    auto upSubarrays = getUpwardSubarrays();

    // if I have the lowest neighbor (i. e. I am the highest rank), I need to send my highest layer in d to them,
    // if I have the highest neighbor (i. e. I am the lowest rank), I can receive it
    int lower, higher;
    getHighestAndLowestNeighbor(d, higher, lower);

    auto success =
        MPI_Sendrecv(this->getData(), 1, upSubarrays[d], lower, TRANSFER_GHOST_LAYER_TAG,
                     this->getData(), 1, downSubarrays[d], higher, 
                     TRANSFER_GHOST_LAYER_TAG, this->getCommunicator(), MPI_STATUS_IGNORE);
    }

  void writeLowerBoundaryToUpperBoundary(DimType d) {
    assert(hasBoundaryPoints_[d] == true);

    // create MPI datatypes
    auto downSubarrays = getDownwardSubarrays();
    auto upSubarrays = getUpwardSubarrays();

    // if I have the highest neighbor (i. e. I am the lowest rank), I need to send my lowest layer in d to them,
    // if I have the lowest neighbor (i. e. I am the highest rank), I can receive it
    int lower, higher;
    getHighestAndLowestNeighbor(d, higher, lower);

    // TODO asynchronous over d??
    auto success =
        MPI_Sendrecv(this->getData(), 1, downSubarrays[d], higher, TRANSFER_GHOST_LAYER_TAG,
                     this->getData(), 1, upSubarrays[d], lower,
                     TRANSFER_GHOST_LAYER_TAG, this->getCommunicator(), MPI_STATUS_IGNORE);
    assert(success == MPI_SUCCESS);
  }

  // non-RVO dependent version of ghost layer exchange
  void exchangeGhostLayerUpward(DimType d, IndexVector& subarrayExtents, std::vector<FG_ELEMENT>& recvbuffer) {
    subarrayExtents = this->getLocalSizes();
    subarrayExtents[d] = 1;

    // create MPI datatype
    auto subarrays = getUpwardSubarrays();

    // if I have a higher neighbor, I need to send my highest layer in d to them,
    // if I have a lower neighbor, I can receive it
    int lower = MPI_PROC_NULL;
    int higher = MPI_PROC_NULL;

    // somehow the cartesian directions in the communicator are reversed
    // cf InitMPI(...)
    auto d_reverse = this->getDimension() - d - 1;
    if (!reverseOrderingDFGPartitions) {
      d_reverse = d;
    }
    MPI_Cart_shift( this->getCommunicator(), d_reverse, 1, &lower, &higher );

    // assert that boundaries have no neighbors (remove in case of periodicity)
    if(this->getLowerBounds()[d] == 0){
      assert(lower < 0);
    }
    if(this->getUpperBounds()[d] == this->getGlobalSizes()[d]){
      assert(higher < 0);
    }

    // create recvbuffer
    auto numElements = std::accumulate(subarrayExtents.begin(), subarrayExtents.end(), 1,
                                       std::multiplies<IndexType>());
    if (lower < 0) {
      numElements = 0;
      subarrayExtents = IndexVector(this->getDimension(), 0);
    }
    recvbuffer = std::vector<FG_ELEMENT>(numElements);

    // TODO asynchronous over d??
    auto success =
        MPI_Sendrecv(this->getData(), 1, subarrays[d], higher, TRANSFER_GHOST_LAYER_TAG,
                     recvbuffer.data(), numElements, this->getMPIDatatype(), lower,
                     TRANSFER_GHOST_LAYER_TAG, this->getCommunicator(), MPI_STATUS_IGNORE);
    assert(success == MPI_SUCCESS);
  }

  std::vector<FG_ELEMENT> exchangeGhostLayerUpward(DimType d, IndexVector& subarrayExtents) {
    std::vector<FG_ELEMENT> recvbuffer{};
    exchangeGhostLayerUpward(d, subarrayExtents, recvbuffer);
    return recvbuffer;
  }

 private:
  /** dimension of the full grid */
  DimType dim_;

  /** the size of the vector, nr of total elements */
  IndexType nrElements_;

  /** the size of the vector, nr of total elements */
  IndexType nrLocalElements_;

  /** levels for each dimension */
  LevelVector levels_;

  /** number of points per axis boundary included*/
  IndexVector nrPoints_;

  /** flag to show if the dimension has boundary points*/
  std::vector<bool> hasBoundaryPoints_;
  /** the offsets in each direction*/
  IndexVector offsets_;

  /** the local offsets in each direction */
  IndexVector localOffsets_;

  /** my partition's lower 1D bounds */
  IndexVector myPartitionsLowerBounds_;

  /** my partition's upper 1D bounds */
  IndexVector myPartitionsUpperBounds_;

  /** the full grid vector, this contains the elements of the full grid */
  std::vector<FG_ELEMENT> fullgridVector_;

  /** Variables for the distributed Full Grid*/
  /** Cartesien MPI Communicator  */
  CommunicatorType communicator_;

 /**
  * the decomposition of the full grid over processors
  * contains (for every dimension) the grid point indices at
  * which a process boundary is assumed
  */
  std::vector<IndexVector> decomposition_;

  /** number of procs in every dimension */
  IndexVector procs_;

  /** the coordinates of each rank on the cartesian communicator*/
  std::vector<std::vector<int>> partitionCoords_;

  /** mpi rank */
  RankType rank_;

  /** mpi size */
  int size_;

  // the MPI Datatypes representing the boundary layers of the MPI processes' subgrid
  std::vector<MPI_Datatype> downwardSubarrays_;
  std::vector<MPI_Datatype> upwardSubarrays_;

  /** number of local (in this grid cell) points per axis*/
  IndexVector nrLocalPoints_;

  // contains for each (local) gridpoint assigment to memory location on the registered dsg
  // use only one of localFGIndexToLocalSGPointerList_ or subspaceIndexToFGIndices_
  std::vector<FG_ELEMENT*> localFGIndexToLocalSGPointerList_;

  // contains for each (local) dsg subspace (.first) the dfg indices that belong to that subspace (.second)
  // if the subspace is not present in either the dsg or the dfg, the entry is not created
  std::vector<std::pair<IndexType, std::vector<IndexType>>> subspaceIndexToFGIndices_;

  DistributedSparseGridUniform<FG_ELEMENT>* dsg_;

  void InitMPI(MPI_Comm comm) {
    MPI_Comm_rank(comm, &rank_);
    MPI_Comm_size(comm, &size_);

    IndexType numSubgrids =
        std::accumulate(procs_.begin(), procs_.end(), 1, std::multiplies<size_t>());

    ASSERT(size_ == static_cast<int>(numSubgrids),
           " size_: " << size_ << " numSubgrids: " << static_cast<int>(numSubgrids));
    assert(size_ == static_cast<int>(numSubgrids));

    // important: note reverse ordering of dims!
    std::vector<int> dims(procs_.rbegin(), procs_.rend());
    if (!reverseOrderingDFGPartitions) {
      dims.assign(procs_.begin(), procs_.end());
    }

    // check if communicator is already cartesian
    int status;
    MPI_Topo_test(comm, &status);

    if (status == MPI_CART) {
      // check if process grid of comm uses the required ordering
      auto maxdims = static_cast<int>(procs_.size());
      std::vector<int> cartdims(maxdims), periods(maxdims), coords(maxdims);
      MPI_Cart_get(comm, static_cast<int>(maxdims), &cartdims[0], &periods[0], &coords[0]);
      ASSERT(cartdims == dims,
            " cartdims: " << cartdims << " dims: " << dims);
      assert(cartdims == dims);

      MPI_Comm_dup(comm, &communicator_);

#ifdef DEBUG_OUTPUT
      std::cout << "DistributedFullGrid: using given cartcomm: " << communicator_ << std::endl;
#endif
    } else {
      // todo mh: think whether periodic bc will be useful
      std::vector<int> periods(dim_, 0);
      int reorder = 0;
      MPI_Cart_create(comm, static_cast<int>(dim_), &dims[0], &periods[0], reorder, &communicator_);
#ifdef DEBUG_OUTPUT
      std::cout << "DistributedFullGrid: create new cartcomm" << std::endl;
#endif
    }
    // fill the partitionCoords_
    {
      partitionCoords_.resize(this->getCommunicatorSize());
      // fill partition coords vector, only once
      for (size_t i = 0; i < this->getCommunicatorSize(); ++i) {
        std::vector<int> tmp(dim_);
        MPI_Cart_coords(communicator_, i, static_cast<int>(dim_), &tmp[0]);

        // important: reverse ordering of partition coords!
        partitionCoords_[i].assign(tmp.rbegin(), tmp.rend());
        if (!reverseOrderingDFGPartitions) {
          partitionCoords_[i].assign(tmp.begin(), tmp.end());
        }
      }
    }
#ifdef DEBUG_OUTPUT
    std::cout << "dfg partition coords "<< partitionCoords_[size_ - 1] << std::endl;
#endif
  }

  /* a regular (equidistant) domain decompositioning for an even number of processes
   * leads to grid points on the (geometrical) process boundaries.
   * with the forwardDecomposition flag it can be decided if the grid points on
   * the process boundaries belong to the process on the right-hand side (true)
   * of the process boundary, or to the one on the left-hand side (false).
   */
  std::vector<IndexVector> getDefaultDecomposition(bool forwardDecomposition) const {
    // create decomposition vectors
    std::vector<IndexVector> decomposition(dim_);

    for (DimType i = 0; i < dim_; ++i) {
      IndexVector& llbnd = decomposition[i];
      llbnd.resize(procs_[i]);

      for (IndexType j = 0; j < procs_[i]; ++j) {
        double tmp = static_cast<double>(nrPoints_[i]) * static_cast<double>(j) /
                     static_cast<double>(procs_[i]);

        if (forwardDecomposition)
          llbnd[j] = static_cast<IndexType>(std::ceil(tmp));
        else
          llbnd[j] = static_cast<IndexType>(std::floor(tmp));
      }
    }
    return decomposition;
  }

  void setDecomposition(const std::vector<IndexVector>& decomposition) {
    assert(decomposition.size() == dim_);
    for (DimType i = 0; i < dim_; ++i)
      assert(decomposition[i].size() == static_cast<size_t>(procs_[i]));

    // check if 1d bounds given in ascending order
    // if not, this might indicate there's something wrong
    for (DimType i = 0; i < dim_; ++i) {
      IndexVector tmp(decomposition[i]);
      std::sort(tmp.begin(), tmp.end());
      assert(tmp == decomposition[i]);
    }

    // check if partitions size larger zero
    // this can be detected by checking for duplicate entries in
    // the 1d decompositions
    for (DimType i = 0; i < dim_; ++i) {
      IndexVector tmp(decomposition[i]);
      IndexVector::iterator last = std::unique(tmp.begin(), tmp.end());

      // todo: this is not sufficient to check for duplicate entries
      // if last points to the same element as tmp.end()
      // this means all entries in tmp are unique
      assert(last == tmp.end() && "some partition in decomposition is of zero size!");
    }

    decomposition_ = decomposition;
  }

  void get1dIndicesLocal(DimType d, const LevelVector& lvec, IndexVector& oneDIndices) {
    LevelType l = lvec[d];

    // get first local idx which has level l
    IndexType start = -1;
    IndexType firstGlobal1dIdx = getFirstGlobal1dIndex(d);

    for (IndexType i = 0; i < nrLocalPoints_[d]; ++i) {
      IndexType global1dIdx = firstGlobal1dIdx + i;

      LevelType myLevel = getLevel(d, global1dIdx);

      // myLevel can be zero if boundary point, but we have our boundary
      // points in level 1 subspaces
      if (myLevel == 0) myLevel = 1;

      if (myLevel == l) {
        start = i;
        break;
      }
    }

    // no point of this level found
    if (start == -1) return;

    IndexType stride;

    // special treatment for level 1 suspaces with boundary
    if (l == 1 && hasBoundaryPoints_[d]) {
      stride = IndexType(std::pow(2, levels_[d] - 1));
    } else {
      stride = IndexType(std::pow(2, levels_[d] - l + 1));
    }

    for (IndexType idx = start; idx < nrLocalPoints_[d]; idx += stride) {
      oneDIndices.push_back(idx);
    }
  }

  // 2d output
  void print2D(std::ostream& os) const {
    assert(dim_ == 2);

    // loop over rows i -> dim2
    for (IndexType i = 0; i < nrLocalPoints_[1]; ++i) {
      IndexType offset = localOffsets_[1] * i;

      for (IndexType j = 0; j < nrLocalPoints_[0]; ++j) {
        os << fullgridVector_[offset + j] << "\t";
      }

      os << std::endl;
    }
  }

  // 1d output
  void print1D(std::ostream& os) const {
    assert(dim_ == 1);

    for (IndexType j = 0; j < nrLocalPoints_[0]; ++j) {
      os << fullgridVector_[j] << "\t";
    }

    os << std::endl;
  }

  // n-dim output. will print 2-dimensional slices of domain
  void print3D(std::ostream& os) const {
    for (IndexType k = 0; k < nrLocalPoints_[2]; ++k) {
      assert(dim_ == 3);

// output global z index
#ifdef DEBUG_OUTPUT
      std::cout << "z = " << this->getLowerBounds()[2] + k << ":" << std::endl;
#endif

      IndexType offsetZ = localOffsets_[2] * k;

      // loop over rows i -> dim2
      for (IndexType i = 0; i < nrLocalPoints_[1]; ++i) {
        IndexType offset = offsetZ + localOffsets_[1] * i;

        for (IndexType j = 0; j < nrLocalPoints_[0]; ++j) {
          os << fullgridVector_[offset + j] << "\t";
        }

        os << std::endl;
      }
    }
  }


};
// end class

// output operator
template <typename FG_ELEMENT>
inline std::ostream& operator<<(std::ostream& os, const DistributedFullGrid<FG_ELEMENT>& dfg) {
  // os << dfg.getRank() << " " << dfg.getDimension() << " " << dfg.getLevels() << std::endl;
  // os << "bounds " << dfg.getLowerBounds() << " to " << dfg.getUpperBounds()  << " to " << dfg.getUpperBoundsCoords() << std::endl;
  // os << "offsets " << dfg.getOffsets() << " " << dfg.getLocalOffsets() << std::endl;
  // os << "sizes " << dfg.getLocalSizes() << "; " << dfg.getElementVector().size() << " " << dfg.getNrLocalElements() << "; " << dfg.getNrElements() << std::endl;
  // std::vector<std::vector<IndexType>> decomposition = dfg.getDecomposition();
  // for (auto& dec : decomposition) {
  //   os << dec;
  // }
  // os << " decomposition; parallelization " << dfg.getParallelization() << std::endl;
  // if(dfg.getNrLocalElements() < 30)
  //   os << "elements " << dfg.getElementVector() << std::endl;

  dfg.print(os);

  return os;
}

}  // namespace combigrid

#endif /* DISTRIBUTEDCOMBIFULLGRID_HPP_ */
