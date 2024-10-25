#ifndef DISTRIBUTEDCOMBIFULLGRID_HPP_
#define DISTRIBUTEDCOMBIFULLGRID_HPP_

#include <algorithm>
#include <bitset>
#include <cassert>
#include <functional>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "fullgrid/FullGrid.hpp"
#include "fullgrid/Tensor.hpp"
#include "mpi/MPICartesianUtils.hpp"
#include "mpi/MPISystem.hpp"
#include "mpi/MPITags.hpp"
#include "sparsegrid/AnyDistributedSparseGrid.hpp"
#include "sparsegrid/DistributedSparseGridUniform.hpp"
#include "utils/DecompositionUtils.hpp"
#include "utils/IndexVector.hpp"
#include "utils/LevelSetUtils.hpp"
#include "utils/LevelVector.hpp"
#include "utils/PowerOfTwo.hpp"
#include "utils/Stats.hpp"
#include "utils/Types.hpp"

namespace combigrid {

/**
 * @brief DistributedFullGrid : a (non-owning) indexed full grid data structure.
 *
 * The DistributedFullGrid accesses a component grid that lives in distributed memory, typically
 * distributed on a cartesian grid within a DisCoTec process group.
 *
 * Because of the distributed-memory character, there are four different ways of indexing the grid:
 * 1. global vector index: d-dimensional index in the global grid, ranges from 0^d to
 * this->getGlobalSizes()
 * 2. local vector index: d-dimensional index in the local part of the global grid, ranges from 0^d
 * to this->getLocalSizes(), which corresponds to this->getLowerBounds to
 * this->getUpperBounds() in the global grid
 * 3. global linear index: linear index in the global grid, ranges from 0 to this->getNrElements()
 * 4. local linear index: linear index in the local grid, ranges from 0 to
 * this->getNrLocalElements() ; these indices do not have to be contiguous in the global grid!
 *
 * The conversion between vector and linear indices is by Fortran order, i.e. elements that differ
 * by 1 in the vector index are the most contiguous in the linear index.
 *
 * For efficient iteration over all the indices in this grid, you should prefer the local linear
 * index.
 *
 * @tparam FG_ELEMENT the data type to be stored on the grid (e.g. double, complex, etc.)
 */
template <typename FG_ELEMENT>
class DistributedFullGrid {
 public:
  /** DistributedFullGrid constructor
   *
   * @param dim the dimensionality of the full grid (may become a template parameter in future
   * versions)
   * @param levels the level vector describing the full grid
   * @param comm the Cartesian communicator to be used for the distributed full grid (will not be
   * duplicated; ownership is not transferred)
   * @param hasBdrPoints a vector of flags to show if the dimension has boundary points (0 for no
   * points, 1 for one-sided boundary, 2 for both sides)
   * @param dataPointer a pointer to the beginning of the data array
   * @param procs a vector of the number of processes in each dimension (must match the
   * decomposition in comm)
   * @param forwardDecomposition a flag to decide if the middle grid points on the process
   * boundaries belong to the process on the right-hand side (true) of the process boundary, or to
   * the one on the left-hand side (false).
   * @param decomposition a vector of the lower bounds of the grid points in each dimension
   */
  DistributedFullGrid(DimType dim, const LevelVector& levels, CommunicatorType const& comm,
                      const std::vector<BoundaryType>& hasBdrPoints, FG_ELEMENT* dataPointer,
                      const std::vector<int>& procs, bool forwardDecomposition = true,
                      const std::vector<IndexVector>& decomposition = std::vector<IndexVector>());

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  DistributedFullGrid(const DistributedFullGrid& other) = delete;
  DistributedFullGrid& operator=(const DistributedFullGrid&) = delete;
  DistributedFullGrid(DistributedFullGrid&& other) = default;
  DistributedFullGrid& operator=(DistributedFullGrid&& other) = default;
#endif  // DOXYGEN_SHOULD_SKIP_THIS

  virtual ~DistributedFullGrid();

  /** get the coordinates on the unit square corresponding to global linear index
   *
   * @param globalIndex [IN] global linear index of the element
   * @param coords [OUT] coordinates on the unit square [0,1]^D
   */
  inline void getCoordsGlobal(IndexType globalLinearIndex, std::vector<real>& coords) const;

  /** get coordinates on the unit square corresponding to local linear index
   *
   * @param localLinearIndex [IN] local linear index of the element
   * @param coords [OUT] coordinates */
  inline void getCoordsLocal(IndexType localLinearIndex, std::vector<real>& coords) const;

  /** get the LI (level,index) notation for a given element in the full grid
   *
   * @param elementIndex [IN] the linear index of the element
   * @param levels [OUT] the levels of the point in the LI notation
   * @param indices [OUT] the indices of the point in the LI notation
   */
  inline void getGlobalLI(IndexType elementIndex, LevelVector& levels, IndexVector& indices) const;

  /** get the global vector index corresponding to a global linear index
   *
   * @param globLinIndex [IN] the global linear index
   * @param globAxisIndex [OUT] the global vector index
   */
  inline void getGlobalVectorIndex(IndexType globLinIndex, IndexVector& globAxisIndex) const;

  /** get the global vector index corresponding to a local vector index
   *
   * @param locAxisIndex [IN] the local vector index
   * @param globAxisIndex [OUT] the global vector index
   */
  inline void getGlobalVectorIndex(const IndexVector& locAxisIndex,
                                   IndexVector& globAxisIndex) const;

  /** get the local vector index corresponding to a local linear index
   *
   * @param locLinIndex [IN] the local linear index
   * @param locAxisIndex [OUT] the local vector index
   */
  inline void getLocalVectorIndex(IndexType locLinIndex, IndexVector& locAxisIndex) const;

  /** get the local vector index corresponding to a global vector index
   *
   * @param globAxisIndex [IN] the global vector index
   * @param locAxisIndex [OUT] the local vector index
   * @return true if global index vector contained in local domain, false otherwise
   */
  inline bool getLocalVectorIndex(const IndexVector& globAxisIndex,
                                  IndexVector& locAxisIndex) const;

  /** get the global linear index corresponding to the global index vector
   *
   * @param axisIndex the vector index
   * @return the global linear index
   */
  inline IndexType getGlobalLinearIndex(const IndexVector& globAxisIndex) const;

  /** get the global linear index corresponding to the local linear index
   *
   * @param locLinIndex the local linear index
   * @return the global linear index
   */
  inline IndexType getGlobalLinearIndex(IndexType locLinIndex) const;

  /** get the local linear index corresponding to the local index vector
   *
   * @param locAxisIndex the local vector index
   * @return the local linear index
   */
  inline IndexType getLocalLinearIndex(const IndexVector& locAxisIndex) const;

  /** get the local linear index corresponding to the global linear index
   *
   * @param globLinIndex the local linear index
   * @return the global linear index, negative value if element not inside local partition
   */
  inline IndexType getLocalLinearIndex(IndexType globLinIndex) const;

  /** is the global index part of this partition?
   *
   * @param globalVectorIndex the global vector index
   * @return true if the global index is part of this' local partition
   */
  inline bool isGlobalIndexHere(IndexVector globalVectorIndex) const;

  /** is the global index part of this partition?
   *
   * @param globLinIndex the global vector index
   * @return true if the global index is part of this' local partition
   */
  inline bool isGlobalIndexHere(IndexType globLinIndex) const;

  /** get the dimension of the full grid */
  inline DimType getDimension() const { return dim_; }

  /** get the strides per dimension for global vector index conversion */
  inline const IndexVector& getOffsets() const { return globalIndexer_.getOffsetsVector(); }

  /** get the strides per dimension for local vector index conversion */
  inline const IndexVector& getLocalOffsets() const { return localTensor_.getOffsetsVector(); }

  /** get the sizes, subsizes, and starts vectors of the local data
   *
   * as required to create a MPI subarray datatype
   */
  inline std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>
  getSizesSubsizesStartsOfSubtensor() const;

  /** get the level vector */
  inline const LevelVector& getLevels() const { return levels_; }

  /** get the grid spacing (sometimes called h)
   *
   * @return the distance between points in each dimension, assuming equidistant grid on unit
   * hypercube
   */
  inline const std::vector<double>& getGridSpacing() const { return gridSpacing_; }

  /** get the inverse grid spacing (1/h) for a given dimension
   *
   * @param inDimension the dimension
   * @return the inverse of getGridSpacing()[inDimension]
   */
  inline double getInverseGridSpacingIn(DimType inDimension) const;

  /** get a vector of inverse grid spacings
   */
  inline std::vector<double> getInverseGridSpacing() const;

  /** get the integral of a "normal" nodal hat basis function on this grid
   *
   * helper function for integration
   *
   * @return the integral of a nodal basis function
   */
  double getInnerNodalBasisFunctionIntegral() const;

  /** get the number of elements in the entire, global full grid */
  inline IndexType getNrElements() const { return static_cast<IndexType>(globalIndexer_.size()); }

  /** get the number of elements in the local partition */
  inline IndexType getNrLocalElements() const {
    return static_cast<IndexType>(localTensor_.size());
  }

  /** get vector extents of the local grid */
  inline const IndexVector& getLocalSizes() const { return localTensor_.getExtentsVector(); }

  /** get vector extents of the global grid */
  inline const IndexVector& getGlobalSizes() const { return globalIndexer_.getExtentsVector(); }

  /** number of points per dimension i */
  inline IndexType globalNumPointsInDimension(DimType i) const { return this->getGlobalSizes()[i]; }

  /** get vector of flags to show how many boundary points this dimension has */
  inline const std::vector<BoundaryType>& returnBoundaryFlags() const { return hasBoundaryPoints_; }

  /** get pointer to the beginning of the local data */
  inline FG_ELEMENT* getData() { return localTensor_.getData(); }

  /** get const pointer to the beginning of the local data */
  inline const FG_ELEMENT* getData() const { return localTensor_.getData(); }

  /** zero out all values in the local storage
   *
   * uses memset for speed
   */
  inline void setZero();

  /** get local Cartesian MPI_Communicator */
  inline CommunicatorType getCommunicator() const { return this->getCartesianUtils().getComm(); }

  /** get the size of the local Cartesian MPI communicator */
  inline int getCommunicatorSize() const { return this->getCartesianUtils().getCommunicatorSize(); }

  /** get the rank in the local Cartesian MPI communicator */
  inline RankType getRank() const { return this->getCartesianUtils().getCommunicatorRank(); }

  /** get lower global vector bounds of the local data in this process */
  inline const IndexVector& getLowerBounds() const;

  /** get lower global vector bounds of the local data in rank r */
  inline IndexVector getLowerBounds(RankType r) const;

  /** get coordinates of this process' lower bounds in some dimension*/
  inline real getLowerBoundsCoord(DimType inDimension) const;

  /** get upper global vector bounds of the local data in this process */
  inline const IndexVector& getUpperBounds() const;

  /** get upper global vector bounds of the local data in rank r */
  inline IndexVector getUpperBounds(RankType r) const;

  /** get coordinates of this process' upper bounds */
  inline std::vector<real> getUpperBoundsCoords() const;

  /** get coordinates of rank r's upper bounds */
  inline std::vector<real> getUpperBoundsCoords(RankType r) const;

  /**
   * @brief get the neighboring process' rank in a given dimension
   *
   * in the sense that the neighbor has the same
   * partion coordinates in all other dimensions than d; and
   * in dimension d it contains the point with the one-dimensional index idx1d
   *
   * @ param dim the dimension in which the neighbor is
   * @ param idx1d the one-dimensional index of the point in the neighbor's partition
   */
  RankType getNeighbor1dFromAxisIndex(DimType dim, IndexType idx1d) const;

  /** get the number of cartesian ranks in every dimension */
  inline const std::vector<int>& getParallelization() const;

  /** get the 1d global index of the first point in the local domain */
  inline IndexType getFirstGlobal1dIndex(DimType d) const { return getLowerBounds()[d]; }

  /** get the 1d global index of the last point in the local domain */
  inline IndexType getLastGlobal1dIndex(DimType d) const { return getUpperBounds()[d] - 1; }

  /** get the (lowest possible) hierarchical level of a global 1d index
   *
   * @param d the dimension
   * @param idx1d the global 1d index
   */
  inline LevelType getLevel(DimType d, IndexType idx1d) const;

  /** get 1d index of the left hierarchical predecessor of a point
   *
   * if this point has level l, return the 1d index of the closest point with lower 1d index and
   * level l-1
   *
   * @param d the dimension
   * @param idx1d the global 1d index
   * @return negative number if point has no left predecessor, else 1d index of the left predecessor
   */
  inline IndexType getLeftPredecessor(DimType d, IndexType idx1d) const;

  /**
   * @brief get 1d index of the right hierarchical predecessor of a point
   *
   * if this point has level l, return the 1d index of the closest point with higher 1d index and
   * level l-1
   *
   * @param d the dimension
   * @param idx1d the global 1d index
   * @return negative number if point has no right predecessor, else 1d index of the right
   * predecessor
   */
  inline IndexType getRightPredecessor(DimType d, IndexType idx1d) const;

  /** get process coordinates of the partition which contains a given point
   *
   * @param globalAxisIndex [IN] the global vector index of the point
   * @param partitionCoords [OUT] the partition coordinates of the process containing the point
   */
  inline void getPartitionCoords(const IndexVector& globalAxisIndex,
                                 std::vector<int>& partitionCoords) const;

  /** print values on the grid
   *
   * (convenience function up to 3d)
   */
  inline void print() const;

  /** get the MPI_Datatype corresponding to FG_ELEMENT */
  inline MPI_Datatype getMPIDatatype() const;

  /** One-D distance between a grid point (by linear index) and a coordinate
   *
   * used as a helper function for interpolation
   *
   * @param oneDimensionalLocalIndex the (local) linear index of the grid point
   * @param coord ND coordinates within the support of the grid point (param 1)
   * @param d the dimension along which the distance is to be calculated
   * @return the distance between the grid point and the coordinate in the specified dimension
   */
  inline double getPointDistanceToCoordinate(IndexType oneDimensionalLocalIndex, double coord,
                                             DimType d) const;

  /** Performs the interpolation of the full grid at the specified coordinates
   *
   * assumes that the localIndex is the highest index whose coordinates on the unit cube are still
   * <= coords
   *
   * @param localIndex the local linear index of the grid point
   * @param coords ND coordinates on the unit square [0,1]^D
   * @return the interpolated value assuming nodal hat basis functions
   */
  inline FG_ELEMENT evalIndexAndAllUpperNeighbors(const IndexVector& localIndex,
                                                  const std::vector<real>& coords) const;

  /** evaluates the full grid at the specified coordinate
   *
   * interpolates by using the nodal basis functions local to this process
   *
   * @param coords ND coordinates on the unit square [0,1]^D
   * @return the interpolated value
   */
  FG_ELEMENT evalLocal(const std::vector<real>& coords) const;

  /** evaluates the full grid at the specified coordinate
   *
   * interpolates by using the nodal basis functions local to this process
   *
   * @param coords [IN] ND coordinates on the unit square [0,1]^D
   * @param value [OUT] the interpolated value
   */
  void evalLocal(const std::vector<real>& coords, FG_ELEMENT& value) const;

  /** evaluates the full grid on a set of specified coordinates
   *
   * @param interpolationCoords vector of ND coordinates on the unit square [0,1]^D
   * @return vector of interpolated values
   */
  std::vector<FG_ELEMENT> getInterpolatedValues(
      const std::vector<std::vector<real>>& interpolationCoords) const;

  /** Gather a full grid from all processes to the root process
   *
   * The root process will have the full grid.
   * Only possible for small grids!
   *
   * @param fg the full grid to gather to
   */
  void gatherFullGrid(FullGrid<FG_ELEMENT>& fg, RankType root);

  /**
   * @brief Get the indices of the points of the subspace on this partition
   *
   * @param l level of hierarchical subspace
   * @return the indices of points on this partition
   */
  inline std::vector<IndexType> getFGPointsOfSubspace(const LevelVector& l) const;

  /**
   * @brief extracts the (hopefully) hierarchical coefficients from dsg
   *        to the full grid's data structure
   *
   * @param dsg the DSG to extract from
   */
  template <bool sparseGridFullyAllocated = true>
  size_t extractFromUniformSG(const DistributedSparseGridUniform<FG_ELEMENT>& dsg);

  /** get the number of points of a given hiararchical 1d-level on this process' grid partition
   *
   * @param l the level
   * @param d the dimension
   * @return the number of 1d points on this partition
   */
  inline IndexType getNumPointsOnThisPartition(LevelType l, DimType d) const;

  /**
   * @brief get the local 1d indices that correspond to a given level
   *
   * @param d the dimension in which we want the indices for the level
   * @param l the level
   * @param oneDIndices [OUT] a list of the local vector indices, length is
   * getNumPointsOnThisPartition(d,l)
   */
  inline void get1dIndicesLocal(DimType d, LevelType l, IndexVector& oneDIndices) const;

  /**
   * @brief Get the Lp Norm of the data on the dfg:
   *      norm = (sum_i(abs(x_i)^p)*integral(basis_i))^1/p
   *
   * @param p : the (polynomial) parameter, 0 is interpreted as maximum norm
   * @return real : the norm
   */
  real getLpNorm(int p) const;

  /** write data to dense binary file using MPI-IO */
  void writePlotFile(const char* filename) const;

  /** write data to legacy-type VTK file using MPI-IO */
  void writePlotFileVTK(const char* filename) const;

  /** get the decomposition vector
   *
   * gives the exact correspondence of 1d indices to cartesian ranks:
   * a vector of the grid point indices at which a process boundary is assumed (for every dimension)
   *
   * @return the decomposition vector of length this->getDimension()
   */
  const std::vector<IndexVector>& getDecomposition() const { return decomposition_; }

  /**
   * @brief get the ranks of the highest and lowest "neighbor" rank in dimension d
   *    only sets highest and lowest if they actually are my neighbors
   */
  void getHighestAndLowestNeighbor(DimType d, int& highest, int& lowest) const;

  /** periodic boundary exchange from lower to upper boundary
   *
   * @param d the dimension which we want to exchange
   */
  void writeLowerBoundaryToUpperBoundary(DimType d);

  /** send the data on the highest-indexed layer in dimension d to the upward-neighboring process
   *  (and receive the data from the downward-neighboring process)
   *
   * @param d the dimension in which we want to exchange
   * @param subarrayExtents [out] the extents of the sent subarray (1 in dimension d and local size
   * in all other dimensions)
   * @return the received data
   */
  std::vector<FG_ELEMENT> exchangeGhostLayerUpward(DimType d, std::vector<int>& subarrayExtents);

  /** non-RVO dependent version of ghost layer exchange */
  void exchangeGhostLayerUpward(DimType d, std::vector<int>& subarrayExtents,
                                std::vector<FG_ELEMENT>& recvbuffer,
                                MPI_Request* recvRequest = nullptr);

  /**
   * @brief check if given globalLinearIndex is on the boundary of this DistributedFullGrid
   */
  std::vector<bool> isGlobalLinearIndexOnBoundary(IndexType globalLinearIndex) const;

  /**
   * @brief check if given localLinearIndex is on the boundary of this DistributedFullGrid
   */
  std::vector<bool> isLocalLinearIndexOnBoundary(IndexType localLinearIndex) const;

  /**
   * @brief recursive helper function for getCornersGlobal*Indices
   *
   */
  std::vector<IndexVector> getCornersGlobalVectorIndicesRecursive(
      std::vector<IndexVector> indicesSoFar, DimType dim) const;

  /**
   * @brief get a vector containing the global vector indices of the 2^d corners of this dfg
   *
   */
  std::vector<IndexVector> getCornersGlobalVectorIndices() const;

  /** get a vector of the values at the grid's global corners
   *
   * @return a vector of the values at the 2^d corners
   */
  std::vector<FG_ELEMENT> getCornersValues() const;

  /** get the MPICartesianUtils associated with this */
  const MPICartesianUtils& getCartesianUtils() const { return cartesianUtils_; }

 protected:
  /** set the data pointer, data needs to be allocated outside */
  inline void setData(FG_ELEMENT* newData) { return localTensor_.setData(newData); }

 private:
  /**
   * @brief sets the MPI-related member cartesianUtils_
   *
   * @param comm the communicator to use (assumed to be cartesian)
   * @param procs the desired partition (for sanity checking)
   */
  void InitMPI(MPI_Comm comm, const std::vector<int>& procs);

  void setDecomposition(const std::vector<IndexVector>& decomposition);

  inline void getFGPointsOfSubspaceRecursive(DimType d, IndexType localLinearIndexSum,
                                             std::vector<IndexVector>& oneDIndices,
                                             std::vector<IndexType>& subspaceIndices) const;

  inline IndexType getStrideForThisLevel(LevelType l, DimType d) const;

  inline IndexType getLocalStartForThisLevel(LevelType l, DimType d,
                                             IndexType strideForThisLevel) const;

  inline IndexType getNumPointsOnThisPartition(DimType d, IndexType localStart,
                                               IndexType strideForThisLevel) const;

  MPI_Datatype getUpwardSubarray(DimType d);

  std::vector<MPI_Datatype> getUpwardSubarrays();

  MPI_Datatype getDownwardSubarray(DimType d);

  std::vector<MPI_Datatype> getDownwardSubarrays();

  /** dimension of the full grid */
  const DimType dim_;

  /** level for each dimension */
  const LevelVector levels_;

  /** the grid spacing h for each dimension */
  std::vector<double> gridSpacing_;

  // TODO: make these normal templates, and SomeDistributedFullGrid a std::variant ?
  /** TensorIndexer , only populated for the used dimensionality**/
  TensorIndexer globalIndexer_{};

  // /** Tensor -- only populated for the used dimensionality**/
  Tensor<FG_ELEMENT> localTensor_{};

  /** flag to show if the dimension has boundary points*/
  const std::vector<BoundaryType> hasBoundaryPoints_;

  /** my partition's lower bounds in global vector coordinates */
  IndexVector myPartitionsLowerBounds_;

  /** my partition's upper bounds in global vector coordinates */
  IndexVector myPartitionsUpperBounds_;

  /** my partition's lower global linear bound */
  IndexType myPartitionsFirstGlobalIndex_;

  /**
   * the decomposition of the full grid over processors
   * contains (for every dimension) the grid point indices at
   * which a process boundary is assumed
   */
  std::vector<IndexVector> decomposition_;

  /** utility to get info about cartesian communicator  */
  static MPICartesianUtils cartesianUtils_;

  // the MPI Datatypes representing the boundary layers of the MPI processes' subgrid
  std::vector<MPI_Datatype> downwardSubarrays_;
  std::vector<MPI_Datatype> upwardSubarrays_;
};

template <typename FG_ELEMENT>
DistributedFullGrid<FG_ELEMENT>::DistributedFullGrid(DimType dim, const LevelVector& levels,
                                                     CommunicatorType const& comm,
                                                     const std::vector<BoundaryType>& hasBdrPoints,
                                                     FG_ELEMENT* dataPointer,
                                                     const std::vector<int>& procs,
                                                     bool forwardDecomposition,
                                                     const std::vector<IndexVector>& decomposition)
    : dim_(dim), levels_(levels), hasBoundaryPoints_(hasBdrPoints) {
  assert(levels_.size() == dim);
  assert(hasBoundaryPoints_.size() == dim);
  assert(procs.size() == dim);

  InitMPI(comm, procs);

  IndexVector nrPoints(dim_);
  gridSpacing_.resize(dim_);

  // set global num of elements and offsets
  IndexType nrElements = 1;
  for (DimType j = 0; j < dim_; j++) {
    nrPoints[j] = combigrid::getNumDofNodal(levels_[j], hasBoundaryPoints_[j]);
    nrElements = nrElements * nrPoints[j];
    if (hasBoundaryPoints_[j] == 1) {
      assert(!decomposition.empty() || !forwardDecomposition);
    }
    assert(levels_[j] < 30);
    gridSpacing_[j] = oneOverPowOfTwo[levels_[j]];
  }

  if (decomposition.size() == 0) {
    setDecomposition(getDefaultDecomposition(
        nrPoints, this->getCartesianUtils().getCartesianDimensions(), forwardDecomposition));
  } else {
    setDecomposition(decomposition);
  }
  globalIndexer_ = TensorIndexer(std::move(nrPoints));
  myPartitionsLowerBounds_ = getLowerBounds(this->getRank());
  myPartitionsFirstGlobalIndex_ = globalIndexer_.sequentialIndex(myPartitionsLowerBounds_);
  myPartitionsUpperBounds_ = getUpperBounds(this->getRank());

  // set local elements and local offsets
  auto nrLocalPoints = getUpperBounds() - getLowerBounds();
  localTensor_ = Tensor<FG_ELEMENT>(dataPointer, std::move(nrLocalPoints));
}

template <typename FG_ELEMENT>
DistributedFullGrid<FG_ELEMENT>::~DistributedFullGrid() {
  for (size_t i = 0; i < upwardSubarrays_.size(); ++i) {
    MPI_Type_free(&upwardSubarrays_[i]);
  }
  for (size_t i = 0; i < downwardSubarrays_.size(); ++i) {
    MPI_Type_free(&downwardSubarrays_[i]);
  }
}


template <typename FG_ELEMENT>
void DistributedFullGrid<FG_ELEMENT>::getCoordsGlobal(IndexType globalLinearIndex,
                                                      std::vector<real>& coords) const {
  IndexType ind = 0;
  IndexType tmp_add = 0;

  coords.resize(dim_);

  for (DimType j = 0; j < dim_; j++) {
    ind = globalLinearIndex % this->getGlobalSizes()[j];
    globalLinearIndex = globalLinearIndex / this->getGlobalSizes()[j];
    // set the coordinate based on if we have boundary points
    tmp_add = (hasBoundaryPoints_[j] > 0) ? (0) : (1);
    coords[j] = static_cast<double>(ind + tmp_add) * getGridSpacing()[j];
  }
}

template <typename FG_ELEMENT>
void DistributedFullGrid<FG_ELEMENT>::getCoordsLocal(IndexType localLinearIndex,
                                                     std::vector<real>& coords) const {
  // todo: probably very inefficient implementation, if crucial for
  // performance implement more direct way of computing the coordinates
  assert(localLinearIndex < getNrLocalElements());

  IndexType globalLinearIndex = getGlobalLinearIndex(localLinearIndex);

  getCoordsGlobal(globalLinearIndex, coords);
}

template <typename FG_ELEMENT>
void DistributedFullGrid<FG_ELEMENT>::getGlobalLI(IndexType elementIndex, LevelVector& levels,
                                                  IndexVector& indices) const {
  IndexType startindex, tmp_val;

  assert(elementIndex < this->getNrElements());
  levels.resize(dim_);
  indices.resize(dim_);

  tmp_val = elementIndex;

  // first calculate intermediary indices
  for (DimType k = 0; k < dim_; k++) {
    startindex = (hasBoundaryPoints_[k] > 0) ? 0 : 1;
    indices[k] = tmp_val % this->getGlobalSizes()[k] + startindex;
    tmp_val = tmp_val / this->getGlobalSizes()[k];
  }

  // The level and index of the element in the hashgridstorage are computed dividing by two the
  // index and level in the fullgrid
  // until we obtain an impair number for the index, thus obtaining the level and index in the
  // hierarchical basis (Aliz Nagy)
  // ...
  for (DimType k = 0; k < dim_; k++) {
    tmp_val = levels_[k];

    if (indices[k] != 0) {
      // todo: these operations can be optimized
      while (indices[k] % 2 == 0) {
        indices[k] = indices[k] / 2;
        tmp_val--;
      }
    } else {
      tmp_val = 0;
    }

    levels[k] = tmp_val;
  }
}

template <typename FG_ELEMENT>
void DistributedFullGrid<FG_ELEMENT>::getGlobalVectorIndex(IndexType globLinIndex,
                                                           IndexVector& globAxisIndex) const {
  assert(globLinIndex < this->getNrElements());
  assert(globAxisIndex.size() == dim_);

  globAxisIndex = this->globalIndexer_.getVectorIndex(globLinIndex);
}

template <typename FG_ELEMENT>
void DistributedFullGrid<FG_ELEMENT>::getGlobalVectorIndex(const IndexVector& locAxisIndex,
                                                           IndexVector& globAxisIndex) const {
  assert(locAxisIndex.size() == dim_);

  globAxisIndex = this->getLowerBounds() + locAxisIndex;
}

template <typename FG_ELEMENT>
void DistributedFullGrid<FG_ELEMENT>::getLocalVectorIndex(IndexType locLinIndex,
                                                          IndexVector& locAxisIndex) const {
  locAxisIndex = this->localTensor_.getVectorIndex(locLinIndex);
}

template <typename FG_ELEMENT>
bool DistributedFullGrid<FG_ELEMENT>::getLocalVectorIndex(const IndexVector& globAxisIndex,
                                                          IndexVector& locAxisIndex) const {
  assert(globAxisIndex.size() == dim_);

  if (this->isGlobalIndexHere(globAxisIndex)) {
    locAxisIndex.assign(globAxisIndex.begin(), globAxisIndex.end());
    std::transform(locAxisIndex.begin(), locAxisIndex.end(), this->getLowerBounds().begin(),
                   locAxisIndex.begin(), std::minus<IndexType>());
    return true;
  } else {
    return false;
  }
}

template <typename FG_ELEMENT>
IndexType DistributedFullGrid<FG_ELEMENT>::getGlobalLinearIndex(
    const IndexVector& globAxisIndex) const {
  return globalIndexer_.sequentialIndex(globAxisIndex);
}

template <typename FG_ELEMENT>
IndexType DistributedFullGrid<FG_ELEMENT>::getGlobalLinearIndex(IndexType locLinIndex) const {
  assert(locLinIndex < this->getNrLocalElements());

  // convert to local vector index
  static thread_local IndexVector locAxisIndex(dim_);
  locAxisIndex.resize(dim_);
  getLocalVectorIndex(locLinIndex, locAxisIndex);

  // convert to global linear index
  IndexType globLinIndex =
      this->myPartitionsFirstGlobalIndex_ + this->getGlobalLinearIndex(locAxisIndex);
  assert(globLinIndex < this->getNrElements());
  return globLinIndex;
}

template <typename FG_ELEMENT>
IndexType DistributedFullGrid<FG_ELEMENT>::getLocalLinearIndex(
    const IndexVector& locAxisIndex) const {
  return localTensor_.sequentialIndex(locAxisIndex);
}

template <typename FG_ELEMENT>
IndexType DistributedFullGrid<FG_ELEMENT>::getLocalLinearIndex(IndexType globLinIndex) const {
  assert(globLinIndex < this->getNrElements());

  // convert to global vector index
  static thread_local IndexVector globAxisIndex(dim_);
  globAxisIndex.resize(dim_);
  getGlobalVectorIndex(globLinIndex, globAxisIndex);

  if (this->isGlobalIndexHere(globAxisIndex)) {
    // convert to local linear index
    return getLocalLinearIndex(globAxisIndex) - this->myPartitionsFirstGlobalIndex_;
  } else {
    return -1;
  }
}

template <typename FG_ELEMENT>
bool DistributedFullGrid<FG_ELEMENT>::isGlobalIndexHere(IndexVector globalVectorIndex) const {
  return (globalVectorIndex >= this->getLowerBounds()) &&
         (globalVectorIndex < this->getUpperBounds());
}

template <typename FG_ELEMENT>
bool DistributedFullGrid<FG_ELEMENT>::isGlobalIndexHere(IndexType globLinIndex) const {
  return isGlobalIndexHere(this->globalIndexer_.getVectorIndex(globLinIndex));
}

template <typename FG_ELEMENT>
std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>
DistributedFullGrid<FG_ELEMENT>::getSizesSubsizesStartsOfSubtensor() const {
  std::vector<int> csizes(this->getGlobalSizes().begin(), this->getGlobalSizes().end());
  std::vector<int> csubsizes(this->getLocalSizes().begin(), this->getLocalSizes().end());
  std::vector<int> cstarts(this->getLowerBounds().begin(), this->getLowerBounds().end());

  return std::make_tuple(csizes, csubsizes, cstarts);
}

template <typename FG_ELEMENT>
double DistributedFullGrid<FG_ELEMENT>::getInverseGridSpacingIn(DimType inDimension) const {
  return powerOfTwo[levels_[inDimension]];
}

template <typename FG_ELEMENT>
std::vector<double> DistributedFullGrid<FG_ELEMENT>::getInverseGridSpacing() const {
  std::vector<double> oneOverH;
  oneOverH.resize(gridSpacing_.size());
  // should be the same as the number of intervals (N)
  for (DimType j = 0; j < dim_; j++) {
    oneOverH[j] = powerOfTwo[levels_[j]];
  }
#ifndef NDEBUG
  std::vector<double> oneOverHByDivision;
  oneOverHByDivision.resize(gridSpacing_.size());
  std::transform(gridSpacing_.begin(), gridSpacing_.end(), oneOverHByDivision.begin(),
                 std::bind(std::divides<double>(), 1, std::placeholders::_1));
  for (DimType j = 0; j < dim_; j++) {
    assert(std::abs(oneOverHByDivision[j] - oneOverH[j]) < 1e-10);
    assert(oneOverHByDivision[j] == oneOverH[j]);
  }
#endif
  return oneOverH;
}

template <typename FG_ELEMENT>
double DistributedFullGrid<FG_ELEMENT>::getInnerNodalBasisFunctionIntegral() const {
  const auto& h = this->getGridSpacing();
  return std::accumulate(h.begin(), h.end(), 1., std::multiplies<double>());
}

template <typename FG_ELEMENT>
void DistributedFullGrid<FG_ELEMENT>::setZero() {
  std::memset(this->getData(), 0, this->getNrLocalElements() * sizeof(FG_ELEMENT));
}

template <typename FG_ELEMENT>
const IndexVector& DistributedFullGrid<FG_ELEMENT>::getLowerBounds() const {
  // return getLowerBounds(this->getRank());
  return myPartitionsLowerBounds_;
}

template <typename FG_ELEMENT>
IndexVector DistributedFullGrid<FG_ELEMENT>::getLowerBounds(RankType r) const {
  assert(r >= 0 && r < this->getCommunicatorSize());
  // get coords of r in cart comm
  IndexVector lowerBounds(dim_);
  const auto& coords = cartesianUtils_.getPartitionCoordsOfRank(r);

  for (DimType i = 0; i < dim_; ++i) {
    lowerBounds[i] = getDecomposition()[i][coords[i]];
  }
  return lowerBounds;
}

template <typename FG_ELEMENT>
real DistributedFullGrid<FG_ELEMENT>::getLowerBoundsCoord(DimType inDimension) const {
  return static_cast<real>(this->getLowerBounds()[inDimension] +
                           (hasBoundaryPoints_[inDimension] > 0 ? 0 : 1)) *
         this->getGridSpacing()[inDimension];
}

template <typename FG_ELEMENT>
const IndexVector& DistributedFullGrid<FG_ELEMENT>::getUpperBounds() const {
  // return getUpperBounds(this->getRank());
  return myPartitionsUpperBounds_;
}

template <typename FG_ELEMENT>
IndexVector DistributedFullGrid<FG_ELEMENT>::getUpperBounds(RankType r) const {
  assert(r >= 0 && r < this->getCommunicatorSize());
  IndexVector upperBounds(dim_);
  const auto& coords = cartesianUtils_.getPartitionCoordsOfRank(r);

  for (DimType i = 0; i < dim_; ++i) {
    RankType n;
    std::vector<int> nc(coords);

    if (nc[i] < this->getCartesianUtils().getCartesianDimensions()[i] - 1) {
      // get rank of next neighbor in dim i
      nc[i] += 1;
      n = this->getCartesianUtils().getRankFromPartitionCoords(nc);
      upperBounds[i] = getLowerBounds(n)[i];
    } else {
      // no neighbor in dim i -> end of domain
      upperBounds[i] = this->getGlobalSizes()[i];
    }
  }
  return upperBounds;
}

template <typename FG_ELEMENT>
std::vector<real> DistributedFullGrid<FG_ELEMENT>::getUpperBoundsCoords() const {
  return getUpperBoundsCoords(this->getRank());
}

template <typename FG_ELEMENT>
std::vector<real> DistributedFullGrid<FG_ELEMENT>::getUpperBoundsCoords(RankType r) const {
  assert(r >= 0 && r < this->getCommunicatorSize());

  /* the upper bound index vector can correspond to coordinates outside of
   * the domain, thus we cannot use the getGlobalLinearIndex method here.
   */
  const IndexVector& ubvi = getUpperBounds(r);
  std::vector<real> coords(dim_);
  IndexType tmp_add = 0;

  for (DimType j = 0; j < dim_; ++j) {
    tmp_add = (hasBoundaryPoints_[j] > 0) ? (0) : (1);
    coords[j] = static_cast<double>(ubvi[j] + tmp_add) * getGridSpacing()[j];
  }
  return coords;
}

template <typename FG_ELEMENT>
RankType DistributedFullGrid<FG_ELEMENT>::getNeighbor1dFromAxisIndex(DimType dim,
                                                                     IndexType idx1d) const {
  // if global index is outside of domain return negative value
  {
    if (idx1d < 0) return -1;
    if (idx1d > this->getGlobalSizes()[dim] - 1) return -1;
  }
  const auto& decomp1d = this->getDecomposition()[dim];
  auto lower = std::lower_bound(decomp1d.begin(), decomp1d.end(), idx1d + 1);
  int partitionIdx1d = static_cast<int>(std::distance(decomp1d.begin(), lower)) - 1;
  RankType r = this->getCartesianUtils().getNeighbor1dFromPartitionIndex(dim, partitionIdx1d);

  // check if global index vector is actually contained in the domain of rank r
  assert(idx1d >= this->getLowerBounds(r)[dim]);
  assert(idx1d < this->getUpperBounds(r)[dim]);
  assert(r < this->getCommunicatorSize());
  return r;
}

template <typename FG_ELEMENT>
const std::vector<int>& DistributedFullGrid<FG_ELEMENT>::getParallelization() const {
  return this->getCartesianUtils().getCartesianDimensions();
}

template <typename FG_ELEMENT>
LevelType DistributedFullGrid<FG_ELEMENT>::getLevel(DimType d, IndexType idx1d) const {
  // IndexVector givec(dim_, 0);
  // givec[d] = idx1d;
  // IndexType idx = getGlobalLinearIndex(givec);
  // LevelVector levels(dim_);
  // IndexVector tmp(dim_);
  // getGlobalLI(idx, levels, tmp);
  // return levels[d];

  if (this->hasBoundaryPoints_[d] == 0) {
    ++idx1d;
  }
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
  LevelType l;
  if constexpr (sizeof(IndexType) == sizeof(long)) {
    l = static_cast<LevelType>(lmax - __builtin_ctzl(static_cast<unsigned long>(idx1d)));
  } else if constexpr (sizeof(IndexType) == sizeof(long long)) {
    l = static_cast<LevelType>(lmax - __builtin_ctzll(static_cast<unsigned long long>(idx1d)));
  } else {
    l = static_cast<LevelType>(lmax - __builtin_ctz(static_cast<unsigned int>(idx1d)));
  }
  return l;
}

template <typename FG_ELEMENT>
IndexType DistributedFullGrid<FG_ELEMENT>::getLeftPredecessor(DimType d, IndexType idx1d) const {
  LevelType l = getLevel(d, idx1d);

  // boundary points
  if (l == 0) return -1;

  LevelType ldiff = static_cast<LevelType>(levels_[d] - l);
  IndexType lpidx = idx1d - combigrid::powerOfTwoByBitshift(ldiff);

  return lpidx;
}

template <typename FG_ELEMENT>
IndexType DistributedFullGrid<FG_ELEMENT>::getRightPredecessor(DimType d, IndexType idx1d) const {
  LevelType l = getLevel(d, idx1d);

  // boundary points
  if (l == 0) return -1;

  LevelType ldiff = static_cast<LevelType>(levels_[d] - l);
  IndexType rpidx = idx1d + combigrid::powerOfTwoByBitshift(ldiff);

  // check if outside of domain
  IndexType numElementsD = this->getGlobalSizes()[d];

  // in case of periodic, "virtual" domain is larger by one
  bool oneSidedBoundary = this->returnBoundaryFlags()[d] == 1;
  if (oneSidedBoundary && rpidx == numElementsD) return numElementsD;

  if (rpidx > numElementsD - 1) rpidx = -1;

  return rpidx;
}

template <typename FG_ELEMENT>
void DistributedFullGrid<FG_ELEMENT>::getPartitionCoords(const IndexVector& globalAxisIndex,
                                                         std::vector<int>& partitionCoords) const {
  partitionCoords.resize(dim_);

  for (DimType d = 0; d < dim_; ++d) {
    partitionCoords[d] = -1;
    for (int i = 0; i < this->getCartesianUtils().getCartesianDimensions()[d]; ++i) {
      if (globalAxisIndex[d] >= getDecomposition()[d][i]) partitionCoords[d] = i;
    }

    // check whether the partition coordinates are valid
    assert(partitionCoords[d] > -1 &&
           partitionCoords[d] < this->getCartesianUtils().getCartesianDimensions()[d]);
  }
}

template <typename FG_ELEMENT>
void DistributedFullGrid<FG_ELEMENT>::print() const {
  std::visit([&](auto&& arg) { print(arg); }, localTensor_);
}

template <typename FG_ELEMENT>
MPI_Datatype DistributedFullGrid<FG_ELEMENT>::getMPIDatatype() const {
  return abstraction::getMPIDatatype(abstraction::getabstractionDataType<FG_ELEMENT>());
}

template <typename FG_ELEMENT>
double DistributedFullGrid<FG_ELEMENT>::getPointDistanceToCoordinate(
    IndexType oneDimensionalLocalIndex, double coord, DimType d) const {
  const auto& h = getGridSpacing();
  auto coordDistance =
      static_cast<double>(this->getLowerBounds()[d] + (hasBoundaryPoints_[d] > 0 ? 0 : 1) +
                          oneDimensionalLocalIndex) *
          h[d] -
      coord;
#ifndef NDEBUG
  if (std::abs(coordDistance) > h[d]) {
    std::cout << "assert bounds " << coordDistance << coord << h << static_cast<int>(d)
              << oneDimensionalLocalIndex << std::endl;
    assert(false &&
           "should only be called for coordinates within the support of this point's basis "
           "function");
  }
#endif  // ndef NDEBUG
  return std::abs(coordDistance);
}

template <typename FG_ELEMENT>
FG_ELEMENT DistributedFullGrid<FG_ELEMENT>::evalIndexAndAllUpperNeighbors(
    const IndexVector& localIndex, const std::vector<real>& coords) const {
  FG_ELEMENT result = 0.;
  auto localLinearIndex = this->getLocalLinearIndex(localIndex);
  for (size_t localIndexIterate = 0;
       localIndexIterate < combigrid::powerOfTwoByBitshift(this->getDimension());
       ++localIndexIterate) {
#ifndef NDEBUG
    auto neighborVectorIndex = localIndex;
#endif
    auto neighborIndex = localLinearIndex;
    real phi_c = 1.;  // value of product of iterate's basis function on coords
    for (DimType d = 0; d < dim_; ++d) {
      const auto lastIndexInDim = this->getLocalSizes()[d] - 1;
      bool isUpperInDim = std::bitset<sizeof(int) * CHAR_BIT>(localIndexIterate).test(d);
      auto iterateIndexInThisDimension = localIndex[d] + isUpperInDim;
      if (iterateIndexInThisDimension < 0) {
        phi_c = 0.;
        break;
      } else if (isUpperInDim && this->hasBoundaryPoints_[d] == 1 &&
                 //  this->getCartesianUtils().isOnLowerBoundaryInDimension(d) &&
                 iterateIndexInThisDimension == this->getGlobalSizes()[d]) {
        // if we have a wrap-around AND are on the lowest cartesian process in this dimension
        // we need to set the index in this dim to 0 and shift the coordinate by 1.
        neighborIndex -= localIndex[d] * this->getLocalOffsets()[d];
        iterateIndexInThisDimension = 0;
        auto coordDistance =
            getPointDistanceToCoordinate(iterateIndexInThisDimension, coords[d] - 1.0, d);
        auto oneOverHinD = powerOfTwo[levels_[d]];
        phi_c *= 1. - coordDistance * oneOverHinD;
      } else if (iterateIndexInThisDimension > lastIndexInDim) {
        phi_c = 0.;
        break;
      } else {
        auto coordDistance =
            getPointDistanceToCoordinate(iterateIndexInThisDimension, coords[d], d);
        auto oneOverHinD = powerOfTwo[levels_[d]];
        phi_c *= 1. - coordDistance * oneOverHinD;
        neighborIndex += isUpperInDim * this->getLocalOffsets()[d];
      }
#ifndef NDEBUG
      neighborVectorIndex[d] = iterateIndexInThisDimension;
#endif
    }
    if (phi_c > 0.) {
#ifndef NDEBUG
      auto unlinearizedNeighborIndex = neighborVectorIndex;
      this->getLocalVectorIndex(neighborIndex, unlinearizedNeighborIndex);
      if (unlinearizedNeighborIndex != neighborVectorIndex) {
        std::cerr << "expected " << unlinearizedNeighborIndex << " got " << neighborVectorIndex
                  << std::endl;
      }
      if (neighborIndex < 0) {
        std::cerr << "expected " << unlinearizedNeighborIndex << " got " << neighborVectorIndex
                  << " or " << neighborIndex << std::endl;
      }
#endif
      assert(neighborIndex > -1);
      assert(neighborIndex < this->getNrLocalElements());
      result += phi_c * this->getData()[neighborIndex];
    }
  }
  return result;
}

template <typename FG_ELEMENT>
FG_ELEMENT DistributedFullGrid<FG_ELEMENT>::evalLocal(const std::vector<real>& coords) const {
  FG_ELEMENT value;
  evalLocal(coords, value);
  return value;
}
template <typename FG_ELEMENT>
void DistributedFullGrid<FG_ELEMENT>::evalLocal(const std::vector<real>& coords,
                                                FG_ELEMENT& value) const {
  assert(coords.size() == this->getDimension());
  // get the lowest-index point of the points
  // whose basis functions contribute to the interpolated value
  const auto& h = getGridSpacing();
  static thread_local IndexVector localIndexLowerNonzeroNeighborPoint;
  localIndexLowerNonzeroNeighborPoint.resize(dim_);
  for (DimType d = 0; d < dim_; ++d) {
#ifndef NDEBUG
    if (coords[d] < 0. || coords[d] > 1.) {
      std::cout << "coords " << coords << " out of bounds" << std::endl;
    }
    assert(coords[d] >= 0. && coords[d] <= 1.);
#endif  // ndef NDEBUG
    // this is the local index of the point that is lower than the coordinate
    // may also be negative if the coordinate is lower than this processes' coordinates
    IndexType localIndexLowerNonzeroNeighborIndexInThisDimension = static_cast<IndexType>(
        std::floor((coords[d] - this->getLowerBoundsCoord(d)) * this->getInverseGridSpacingIn(d)));

    // check if we even need to evaluate on this process
    if (localIndexLowerNonzeroNeighborIndexInThisDimension < -1) {
      // index too small
      value = 0.;
      return;
    } else if ((coords[d] >= 1.0 - h[d]) && this->hasBoundaryPoints_[d] == 1 &&
               this->getCartesianUtils().isOnLowerBoundaryInDimension(d)) {
      // if we have periodic boundary and this process is at the lower end of the dimension d
      // we need the periodic coordinate => don't return 0
    } else if (localIndexLowerNonzeroNeighborIndexInThisDimension > this->getLocalSizes()[d] - 1) {
      // index too high
      value = 0.;
      return;
    }
    localIndexLowerNonzeroNeighborPoint[d] = localIndexLowerNonzeroNeighborIndexInThisDimension;
  }
  // evaluate at those points and sum up according to the basis function
  value = evalIndexAndAllUpperNeighbors(localIndexLowerNonzeroNeighborPoint, coords);
}

template <typename FG_ELEMENT>
std::vector<FG_ELEMENT> DistributedFullGrid<FG_ELEMENT>::getInterpolatedValues(
    const std::vector<std::vector<real>>& interpolationCoords) const {
  auto numValues = interpolationCoords.size();
  std::vector<FG_ELEMENT> values;
  values.resize(numValues);
#pragma omp parallel for default(none) firstprivate(numValues) shared(values, interpolationCoords) \
    schedule(static)
  for (size_t i = 0; i < numValues; ++i) {
    this->evalLocal(interpolationCoords[i], values[i]);
  }
  MPI_Allreduce(MPI_IN_PLACE, values.data(), static_cast<int>(numValues), this->getMPIDatatype(),
                MPI_SUM, this->getCommunicator());
  return values;
}

template <typename FG_ELEMENT>
void DistributedFullGrid<FG_ELEMENT>::gatherFullGrid(FullGrid<FG_ELEMENT>& fg, RankType root) {
  int size = this->getCommunicatorSize();
  int rank = this->getRank();
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

      std::vector<int> csizes(sizes.begin(), sizes.end());
      std::vector<int> csubsizes(subsizes.begin(), subsizes.end());
      std::vector<int> cstarts(starts.begin(), starts.end());

      // create subarray view on data
      MPI_Datatype mysubarray;
      MPI_Type_create_subarray(static_cast<int>(this->getDimension()), &csizes[0], &csubsizes[0],
                               &cstarts[0], MPI_ORDER_FORTRAN, this->getMPIDatatype(), &mysubarray);
      MPI_Type_commit(&mysubarray);
      subarrayTypes.push_back(mysubarray);

      int src = r;
      MPI_Request req;
      MPI_Irecv(fg.getData(), 1, mysubarray, src, 0, this->getCommunicator(), &req);
      requests.push_back(req);
    }
  }

  MPI_Wait(&sendRequest, MPI_STATUS_IGNORE);

  if (rank == root) {
    MPI_Waitall(static_cast<int>(requests.size()), &requests[0], MPI_STATUSES_IGNORE);
  }

  for (size_t i = 0; i < subarrayTypes.size(); ++i) MPI_Type_free(&subarrayTypes[i]);
}

template <typename FG_ELEMENT>
std::vector<IndexType> DistributedFullGrid<FG_ELEMENT>::getFGPointsOfSubspace(
    const LevelVector& l) const {
  IndexVector subspaceIndices;
  IndexType numPointsOfSubspace = 1;
  static thread_local std::vector<IndexVector> oneDIndices;
  oneDIndices.resize(dim_);
  for (DimType d = 0; d < dim_; ++d) {
    if (l[d] > levels_[d]) {
      return subspaceIndices;
    }
    get1dIndicesLocal(d, l[d], oneDIndices[d]);
    numPointsOfSubspace *= static_cast<decltype(numPointsOfSubspace)>(oneDIndices[d].size());
  }
  if (numPointsOfSubspace > 0) {
    subspaceIndices.reserve(numPointsOfSubspace);

    IndexType localLinearIndexSum = 0;
    getFGPointsOfSubspaceRecursive(static_cast<DimType>(dim_ - 1), localLinearIndexSum, oneDIndices,
                                   subspaceIndices);
  }
  assert(static_cast<IndexType>(subspaceIndices.size()) == numPointsOfSubspace);
  return subspaceIndices;
}

template <typename FG_ELEMENT>
template <bool sparseGridFullyAllocated>
size_t DistributedFullGrid<FG_ELEMENT>::extractFromUniformSG(
    const DistributedSparseGridUniform<FG_ELEMENT>& dsg) {
  assert(dsg.isSubspaceDataCreated());

  // all the hierarchical subspaces contained in this full grid
  const auto downwardClosedSet = combigrid::getDownSet(levels_);

  // loop over all subspaces (-> somewhat linear access in the sg)
  size_t numCopied = 0;
  static thread_local IndexVector subspaceIndices;
  typename AnyDistributedSparseGrid::SubspaceIndexType sIndex = 0;
#pragma omp parallel for shared(dsg, downwardClosedSet) default(none) schedule(guided) \
    firstprivate(sIndex) reduction(+ : numCopied)
  for (const auto& level : downwardClosedSet) {
    sIndex = dsg.getIndexInRange(level, sIndex);
    bool shouldBeCopied = sIndex > -1 && dsg.getDataSize(sIndex) > 0;
    if constexpr (!sparseGridFullyAllocated) {
      shouldBeCopied = shouldBeCopied && dsg.isSubspaceCurrentlyAllocated(sIndex);
    }
    if (shouldBeCopied) {
      auto sPointer = dsg.getData(sIndex);
      subspaceIndices = std::move(this->getFGPointsOfSubspace(level));
      // #pragma omp simd linear(sPointer : 1) // no simd benefit
      for (size_t fIndex = 0; fIndex < subspaceIndices.size(); ++fIndex) {
        this->getData()[subspaceIndices[fIndex]] = *sPointer;
        ++sPointer;
      }
      assert(dsg.getDataSize(sIndex) == subspaceIndices.size());
      numCopied += subspaceIndices.size();
    }
  }
  return numCopied;
}

template <typename FG_ELEMENT>
void DistributedFullGrid<FG_ELEMENT>::get1dIndicesLocal(DimType d, LevelType l,
                                                        IndexVector& oneDIndices) const {
  assert(l > 0);
  assert(d < this->getDimension());
  if (l > levels_[d]) {
    oneDIndices.clear();
    return;
  }

  const IndexType strideForThisLevel = getStrideForThisLevel(l, d);
  IndexType localStart = getLocalStartForThisLevel(l, d, strideForThisLevel);
  const IndexType numPointsOnThisPartition =
      getNumPointsOnThisPartition(d, localStart, strideForThisLevel);
  oneDIndices.resize(numPointsOnThisPartition);

  std::generate(oneDIndices.begin(), oneDIndices.end(), [&localStart, &strideForThisLevel]() {
    auto localStartBefore = localStart;
    localStart += strideForThisLevel;
    return localStartBefore;
  });
}

template <typename FG_ELEMENT>
real DistributedFullGrid<FG_ELEMENT>::getLpNorm(int p) const {
  assert(p >= 0);
  // special case maximum norm
  MPI_Datatype dtype = abstraction::getMPIDatatype(abstraction::getabstractionDataType<real>());
  if (p == 0) {
    real max = 0.0;
    auto data = this->getData();
#pragma omp parallel for reduction(max : max) default(none) shared(data) schedule(static)
    for (IndexType i = 0; i < this->getNrLocalElements(); ++i) {
      max = std::max(max, std::abs(data[i]));
    }
    real globalMax(-1);
    MPI_Allreduce(&max, &globalMax, 1, dtype, MPI_MAX, getCommunicator());
    return globalMax;
  } else {
    real p_f = static_cast<real>(p);
    auto data = this->getData();
    real res = 0.0;

#pragma omp parallel for reduction(+ : res) default(none) shared(data, oneOverPowOfTwo) \
    firstprivate(p_f) schedule(static)
    for (IndexType i = 0; i < this->getNrLocalElements(); ++i) {
      auto isOnBoundary = this->isLocalLinearIndexOnBoundary(i);
      auto countBoundary = std::count(isOnBoundary.begin(), isOnBoundary.end(), true);
      double hatFcnIntegral = oneOverPowOfTwo[countBoundary];
      // std::cout << hatFcnIntegral << std::endl;
      // sumIntegral += hatFcnIntegral;
      real abs = std::abs(data[i]);
      res += std::pow(abs, p_f) * hatFcnIntegral;
    }
    // std::cout << " sum " << sumIntegral << std::endl;
    real globalRes(0.);
    MPI_Allreduce(&res, &globalRes, 1, dtype, MPI_SUM, getCommunicator());
    // TODO this is only correct for 1 norm
    // other norms have mixed terms and need additional boundary communication
    return std::pow(globalRes * std::pow(this->getInnerNodalBasisFunctionIntegral(), p_f),
                    1.0 / p_f);
  }
}

template <typename FG_ELEMENT>
void DistributedFullGrid<FG_ELEMENT>::writePlotFile(const char* filename) const {
  auto dim = getDimension();

  // create subarray data type
  auto [csizes, csubsizes, cstarts] = this->getSizesSubsizesStartsOfSubtensor();

  // create subarray view on data
  MPI_Datatype mysubarray;
  MPI_Type_create_subarray(static_cast<int>(getDimension()), &csizes[0], &csubsizes[0], &cstarts[0],
                           MPI_ORDER_FORTRAN, getMPIDatatype(), &mysubarray);
  MPI_Type_commit(&mysubarray);

  // open file
  MPI_File fh;
  MPI_File_open(getCommunicator(), filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

  MPI_Offset offset = 0;

  // .raw files can be read by paraview, in that case write the header separately
  std::string fn = filename;
  if (fn.find(".raw") != std::string::npos) {
    // rank 0 write human-readable header
    if (this->getRank() == 0) {
      auto headername = fn + "_header";
      std::ofstream ofs(headername);

      // first line: dimension
      ofs << "dimensionality " << getDimension() << std::endl;

      // grid points per dimension in the order
      // x_0, x_1, ... , x_d
      ofs << "Extents ";
      for (auto s : csizes) {
        ofs << s << " ";
      }
      ofs << std::endl;

      // data type
      ofs << "Data type size " << sizeof(FG_ELEMENT) << std::endl;
      // TODO (pollinta) write endianness and data spacing
    }
  } else {
    // rank 0 write dim and resolution (and data format?)
    if (this->getRank() == 0) {
      MPI_File_write(fh, &dim, 1, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);

      std::vector<int> res(csizes.begin(), csizes.end());
      MPI_File_write(fh, &res[0], 6, MPI_INT, MPI_STATUS_IGNORE);
    }

    // set file view to right offset (in bytes)
    offset = (1 + dim) * sizeof(int);
  }

  MPI_File_set_view(fh, offset, getMPIDatatype(), mysubarray, "native", MPI_INFO_NULL);

  // write subarray
  MPI_File_write_all(fh, getData(), static_cast<int>(getNrLocalElements()), getMPIDatatype(),
                     MPI_STATUS_IGNORE);
  // close file
  MPI_File_close(&fh);
  MPI_Type_free(&mysubarray);
}

template <typename FG_ELEMENT>
void DistributedFullGrid<FG_ELEMENT>::writePlotFileVTK(const char* filename) const {
  auto dim = getDimension();
  assert(dim < 4);  // vtk supports only up to 3D

  // create subarray data type
  auto [csizes, csubsizes, cstarts] = this->getSizesSubsizesStartsOfSubtensor();

  // create subarray view on data
  MPI_Datatype mysubarray;
  MPI_Type_create_subarray(static_cast<int>(getDimension()), &csizes[0], &csubsizes[0], &cstarts[0],
                           MPI_ORDER_FORTRAN, getMPIDatatype(), &mysubarray);
  MPI_Type_commit(&mysubarray);

  // open file
  MPI_File fh;
  MPI_File_open(getCommunicator(), filename, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

  std::stringstream vtk_header;
  // cf https://lorensen.github.io/VTKExamples/site/VTKFileFormats/ => structured points
  vtk_header << "# vtk DataFile Version 2.0.\n"
             << "This file contains the combination solution evaluated on a full grid\n"
             << "BINARY\n"
             << "DATASET STRUCTURED_POINTS\n";
  if (dim == 3) {  // TODO change for non-boundary grids using getGridSpacing
    vtk_header << "DIMENSIONS " << csizes[0] << " " << csizes[1] << " " << csizes[2] << "\n"
               << "ORIGIN 0 0 0\n"
               << "SPACING " << 1. / static_cast<double>(csizes[0] - 1) << " "
               << 1. / static_cast<double>(csizes[1] - 1) << " "
               << 1. / static_cast<double>(csizes[2] - 1) << "\n";
  } else if (dim == 2) {
    vtk_header << "DIMENSIONS " << csizes[0] << " " << csizes[1] << " 1\n"
               << "ORIGIN 0 0 0\n"
               << "SPACING " << 1. / static_cast<double>(csizes[0] - 1) << " "
               << 1. / static_cast<double>(csizes[1] - 1) << " 1\n";
  } else if (dim == 1) {
    vtk_header << "DIMENSIONS " << csizes[0] << " 1 1\n"
               << "ORIGIN 0 0 0\n"
               << "SPACING " << 1. / static_cast<double>(csizes[0] - 1) << " 1 1\n";
  } else {
    assert(false);
  }
  vtk_header << "POINT_DATA "
             << std::accumulate(csizes.begin(), csizes.end(), 1, std::multiplies<int>()) << "\n"
             << "SCALARS quantity double 1\n"
             << "LOOKUP_TABLE default\n";
  // TODO set the right data type from combidatatype, for now double by default
  [[maybe_unused]] bool rightDataType = std::is_same<CombiDataType, double>::value;
  assert(rightDataType);
  auto header_string = vtk_header.str();
  auto header_size = header_string.size();

  // rank 0 write header
  if (this->getRank() == 0) {
    MPI_File_write(fh, header_string.data(), static_cast<int>(header_size), MPI_CHAR,
                   MPI_STATUS_IGNORE);
  }

  // set file view to right offset (in bytes)
  MPI_Offset offset = header_size * sizeof(char);
  // external32 not supported in OpenMPI < 5. -> writes "native" endianness
  // might work with MPICH
  MPI_File_set_view(fh, offset, getMPIDatatype(), mysubarray, "external32", MPI_INFO_NULL);

  // write subarray
  MPI_File_write_all(fh, getData(), static_cast<int>(getNrLocalElements()), getMPIDatatype(),
                     MPI_STATUS_IGNORE);

  // close file
  MPI_File_close(&fh);
  MPI_Type_free(&mysubarray);
}

template <typename FG_ELEMENT>
void DistributedFullGrid<FG_ELEMENT>::getHighestAndLowestNeighbor(DimType d, int& highest,
                                                                  int& lowest) const {
  lowest = MPI_PROC_NULL;
  highest = MPI_PROC_NULL;

  auto d_reverse = this->getDimension() - d - 1;
  if (!reverseOrderingDFGPartitions) {
    d_reverse = d;
  }

  MPI_Cart_shift(this->getCommunicator(), static_cast<int>(d_reverse),
                 static_cast<int>(getParallelization()[d] - 1), &lowest, &highest);

  // assert only boundaries have those neighbors (remove in case of periodicity)
  // this assumes no periodicity!
  if (!this->getLowerBounds()[d] == 0) {
    assert(highest < 0);
  }
  if (!this->getUpperBounds()[d] == this->getGlobalSizes()[d]) {
    assert(lowest < 0);
  }
}

template <typename FG_ELEMENT>
void DistributedFullGrid<FG_ELEMENT>::writeLowerBoundaryToUpperBoundary(DimType d) {
  assert(hasBoundaryPoints_[d] == 2);

  // create MPI datatypes
  auto downSubarray = getDownwardSubarray(d);
  auto upSubarray = getUpwardSubarray(d);

  // if I have the highest neighbor (i. e. I am the lowest rank), I need to send my lowest layer
  // in d to them, if I have the lowest neighbor (i. e. I am the highest rank), I can receive it
  int lower, higher;
  getHighestAndLowestNeighbor(d, higher, lower);

  // TODO asynchronous over d??
  [[maybe_unused]] auto success = MPI_Sendrecv(
      this->getData(), 1, downSubarray, higher, TRANSFER_GHOST_LAYER_TAG, this->getData(), 1,
      upSubarray, lower, TRANSFER_GHOST_LAYER_TAG, this->getCommunicator(), MPI_STATUS_IGNORE);
  assert(success == MPI_SUCCESS);
  MPI_Type_free(&downSubarray);
  MPI_Type_free(&upSubarray);
}

template <typename FG_ELEMENT>
void DistributedFullGrid<FG_ELEMENT>::exchangeGhostLayerUpward(DimType d,
                                                               std::vector<int>& subarrayExtents,
                                                               std::vector<FG_ELEMENT>& recvbuffer,
                                                               MPI_Request* recvRequest) {
  subarrayExtents.resize(this->getDimension());
  subarrayExtents.assign(this->getLocalSizes().begin(), this->getLocalSizes().end());
  subarrayExtents[d] = 1;

  // create MPI datatype
  auto subarray = getUpwardSubarray(d);

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
  MPI_Cart_shift(this->getCommunicator(), static_cast<int>(d_reverse), 1, &lower, &higher);

  if (this->returnBoundaryFlags()[d] == 1) {
    // if one-sided boundary, assert that everyone has neighbors
    assert(lower >= 0);
    assert(higher >= 0);
  } else if (this->returnBoundaryFlags()[d] == 2) {
    // assert that boundaries have no neighbors
    if (this->getLowerBounds()[d] == 0) {
      assert(lower < 0);
    }
    if (this->getUpperBounds()[d] == this->getGlobalSizes()[d]) {
      assert(higher < 0);
    }
  }

  // create recvbuffer
  auto numElements = std::accumulate(subarrayExtents.begin(), subarrayExtents.end(), 1,
                                     std::multiplies<IndexType>());
  if (lower < 0) {
    numElements = 0;
    std::memset(subarrayExtents.data(), 0, subarrayExtents.size() * sizeof(int));
  }
  recvbuffer.resize(numElements);
  std::memset(recvbuffer.data(), 0, recvbuffer.size() * sizeof(FG_ELEMENT));

  if (recvRequest != nullptr) {
    [[maybe_unused]] auto success =
        MPI_Irecv(recvbuffer.data(), numElements, this->getMPIDatatype(), lower,
                  TRANSFER_GHOST_LAYER_TAG, this->getCommunicator(), recvRequest);
    assert(success == MPI_SUCCESS);
    success = MPI_Send(this->getData(), 1, subarray, higher, TRANSFER_GHOST_LAYER_TAG,
                       this->getCommunicator());
    assert(success == MPI_SUCCESS);
  } else {
    [[maybe_unused]] auto success =
        MPI_Sendrecv(this->getData(), 1, subarray, higher, TRANSFER_GHOST_LAYER_TAG,
                     recvbuffer.data(), numElements, this->getMPIDatatype(), lower,
                     TRANSFER_GHOST_LAYER_TAG, this->getCommunicator(), MPI_STATUS_IGNORE);
    assert(success == MPI_SUCCESS);
  }
  MPI_Type_free(&subarray);
}

template <typename FG_ELEMENT>
std::vector<FG_ELEMENT> DistributedFullGrid<FG_ELEMENT>::exchangeGhostLayerUpward(
    DimType d, std::vector<int>& subarrayExtents) {
  std::vector<FG_ELEMENT> recvbuffer{};
  exchangeGhostLayerUpward(d, subarrayExtents, recvbuffer);
  return recvbuffer;
}

template <typename FG_ELEMENT>
std::vector<bool> DistributedFullGrid<FG_ELEMENT>::isGlobalLinearIndexOnBoundary(
    IndexType globalLinearIndex) const {
  // this could likely be done way more efficiently, but it's
  // currently not in any performance critical spot
  std::vector<bool> isOnBoundary(this->getDimension(), false);

  // convert to global vector index
  IndexVector globalAxisIndex(dim_);
  getGlobalVectorIndex(globalLinearIndex, globalAxisIndex);
  for (DimType d = 0; d < this->getDimension(); ++d) {
    if (this->returnBoundaryFlags()[d] == 2) {
      if (globalAxisIndex[d] == 0 ||
          globalAxisIndex[d] == this->globalNumPointsInDimension(d) - 1) {
        isOnBoundary[d] = true;
      }
    }
  }
  return isOnBoundary;
}

template <typename FG_ELEMENT>
std::vector<bool> DistributedFullGrid<FG_ELEMENT>::isLocalLinearIndexOnBoundary(
    IndexType localLinearIndex) const {
  return isGlobalLinearIndexOnBoundary(getGlobalLinearIndex(localLinearIndex));
}
template <typename FG_ELEMENT>
std::vector<IndexVector> DistributedFullGrid<FG_ELEMENT>::getCornersGlobalVectorIndicesRecursive(
    std::vector<IndexVector> indicesSoFar, DimType dim) const {
  if (dim < this->getDimension()) {
    std::vector<IndexVector> newIndicesSoFar{};
    for (const auto& indexVec : indicesSoFar) {
      newIndicesSoFar.push_back(indexVec);
      newIndicesSoFar.back().push_back(0);
      newIndicesSoFar.push_back(indexVec);
      newIndicesSoFar.back().push_back(this->globalNumPointsInDimension(dim) - 1);
    }
    assert(newIndicesSoFar.size() == 2 * indicesSoFar.size());
    return getCornersGlobalVectorIndicesRecursive(newIndicesSoFar, static_cast<DimType>(dim + 1));
  } else {
    return indicesSoFar;
  }
}

template <typename FG_ELEMENT>
std::vector<IndexVector> DistributedFullGrid<FG_ELEMENT>::getCornersGlobalVectorIndices() const {
  auto emptyVectorInVector = std::vector<IndexVector>(1);
  assert(emptyVectorInVector.size() == 1);
  assert(emptyVectorInVector[0].size() == 0);
  std::vector<IndexVector> cornersVectors =
      getCornersGlobalVectorIndicesRecursive(emptyVectorInVector, 0);
  assert(cornersVectors.size() == static_cast<size_t>(powerOfTwo[this->getDimension()]));
  return cornersVectors;
}

template <typename FG_ELEMENT>
std::vector<FG_ELEMENT> DistributedFullGrid<FG_ELEMENT>::getCornersValues() const {
  std::vector<FG_ELEMENT> values(powerOfTwo[this->getDimension()]);
  auto corners = getCornersGlobalVectorIndices();
  for (size_t cornerNo = 0; cornerNo < corners.size(); ++cornerNo) {
    if (this->isGlobalIndexHere(corners[cornerNo])) {
      // convert to local vector index, then to linear index
      IndexVector locAxisIndex(this->getDimension());
      [[maybe_unused]] bool present = getLocalVectorIndex(corners[cornerNo], locAxisIndex);
      assert(present);
      auto index = getLocalLinearIndex(locAxisIndex);
      values[cornerNo] = this->getData()[index];
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, values.data(), static_cast<int>(values.size()),
                this->getMPIDatatype(), MPI_SUM, this->getCommunicator());
  return values;
}

template <typename FG_ELEMENT>
void DistributedFullGrid<FG_ELEMENT>::InitMPI(MPI_Comm comm, const std::vector<int>& procs) {
  // check if communicator is already cartesian
  int status;
  MPI_Topo_test(comm, &status);

  if (status == MPI_CART) {
    if (comm != cartesianUtils_.getComm()) {
      assert(uniformDecomposition);
      cartesianUtils_ = MPICartesianUtils(comm);
    }
    if (procs != cartesianUtils_.getCartesianDimensions()) {
      std::cerr << "unpart" << std::endl;
      throw std::runtime_error("The given communicator is not partitioned as desired");
    }
  } else {
    // MPI_Comm_dup(comm, &communicator_);
    // cf. https://www.researchgate.net/publication/220439585_MPI_on_millions_of_cores
    // "Figure 3 shows the memory consumption in all these cases after 32 calls to MPI Comm dup"
    std::cerr << "undup" << std::endl;
    throw std::runtime_error(
        "Currently testing to not duplicate communicator (to save memory), \
                          if you do want to use this code please take care that \
                          MPI_Comm_free(&communicator_) will be called at some point");
  }
}

template <typename FG_ELEMENT>
void DistributedFullGrid<FG_ELEMENT>::setDecomposition(
    const std::vector<IndexVector>& decomposition) {
#ifndef NDEBUG
  assert(decomposition.size() == dim_);
  for (DimType i = 0; i < dim_; ++i)
    assert(decomposition[i].size() ==
           static_cast<size_t>(this->getCartesianUtils().getCartesianDimensions()[i]));

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
#endif  // not def NDEBUG

  decomposition_ = decomposition;
}

template <typename FG_ELEMENT>
void DistributedFullGrid<FG_ELEMENT>::getFGPointsOfSubspaceRecursive(
    DimType d, IndexType localLinearIndexSum, std::vector<IndexVector>& oneDIndices,
    std::vector<IndexType>& subspaceIndices) const {
  assert(d < dim_);
  assert(oneDIndices.size() == dim_);
  assert(!oneDIndices.empty());

  for (const auto idx : oneDIndices[d]) {
    auto updatedLocalIndexSum = localLinearIndexSum;
    updatedLocalIndexSum += this->getLocalOffsets()[d] * idx;
    if (d > 0) {
      getFGPointsOfSubspaceRecursive(static_cast<DimType>(d - 1), updatedLocalIndexSum, oneDIndices,
                                     subspaceIndices);
    } else {
      subspaceIndices.emplace_back(updatedLocalIndexSum);
    }
  }
}

template <typename FG_ELEMENT>
IndexType DistributedFullGrid<FG_ELEMENT>::getStrideForThisLevel(LevelType l, DimType d) const {
  assert(d < this->getDimension());
  // special treatment for level 1 suspaces with boundary
  return (l == 1 && hasBoundaryPoints_[d] > 0)
             ? combigrid::powerOfTwoByBitshift(static_cast<LevelType>(levels_[d] - 1))
             : combigrid::powerOfTwoByBitshift(static_cast<LevelType>(levels_[d] - l + 1));
}

template <typename FG_ELEMENT>
IndexType DistributedFullGrid<FG_ELEMENT>::getLocalStartForThisLevel(
    LevelType l, DimType d, IndexType strideForThisLevel) const {
  const auto firstGlobal1dIdx = getFirstGlobal1dIndex(d);

  // get global offset to find indices of this level
  // this is the first global index that has level l in dimension d
  IndexType offsetForThisLevel;
  if (hasBoundaryPoints_[d] > 0) {
    if (l == 1) {
      offsetForThisLevel = 0;
    } else {
      // offsetForThisLevel = strideForThisLevel / 2;
      offsetForThisLevel = combigrid::powerOfTwoByBitshift(static_cast<LevelType>(levels_[d] - l));
    }
  } else {
    // offsetForThisLevel = strideForThisLevel / 2 - 1;
    offsetForThisLevel =
        combigrid::powerOfTwoByBitshift(static_cast<LevelType>(levels_[d] - l)) - 1;
  }
  assert(offsetForThisLevel > -1);
  assert(strideForThisLevel >= offsetForThisLevel);

  // first level-index that is on our partition
  const IndexType firstGlobalIndexOnThisPartition =
      (firstGlobal1dIdx < offsetForThisLevel)
          ? 0
          : (firstGlobal1dIdx - offsetForThisLevel + strideForThisLevel - 1) / strideForThisLevel;
  // get global index of first local index which has level l
  const IndexType globalStart =
      (firstGlobalIndexOnThisPartition)*strideForThisLevel + offsetForThisLevel;
  assert(getLevel(d, globalStart) == l ||
         (getLevel(d, globalStart) == 0 && hasBoundaryPoints_[d] > 0 && l == 1));
  assert(globalStart >= firstGlobal1dIdx);

  return globalStart - firstGlobal1dIdx;
}

template <typename FG_ELEMENT>
IndexType DistributedFullGrid<FG_ELEMENT>::getNumPointsOnThisPartition(
    DimType d, IndexType localStart, IndexType strideForThisLevel) const {
  assert(d < this->getDimension());
  return (this->getLocalSizes()[d] - 1 < localStart)
             ? 0
             : (this->getLocalSizes()[d] - 1 - localStart) / strideForThisLevel + 1;
}

template <typename FG_ELEMENT>
IndexType DistributedFullGrid<FG_ELEMENT>::getNumPointsOnThisPartition(LevelType l,
                                                                       DimType d) const {
  assert(!(l > levels_[d]));
  const auto strideForThisLevel = getStrideForThisLevel(l, d);
  return getNumPointsOnThisPartition(d, getLocalStartForThisLevel(l, d, strideForThisLevel),
                                     strideForThisLevel);
}

template <typename FG_ELEMENT>
MPI_Datatype DistributedFullGrid<FG_ELEMENT>::getUpwardSubarray(DimType d) {
  // do index calculations
  // set lower bounds of subarray
  auto subarrayLowerBounds = this->getLowerBounds();
  auto subarrayUpperBounds = this->getUpperBounds();
  subarrayLowerBounds[d] += this->getLocalSizes()[d] - 1;

  auto subarrayExtents = subarrayUpperBounds - subarrayLowerBounds;
  assert(subarrayExtents[d] == 1);
  auto subarrayStarts = subarrayLowerBounds - this->getLowerBounds();

  // create MPI datatype
  std::vector<int> sizes(this->getLocalSizes().begin(), this->getLocalSizes().end());
  std::vector<int> subsizes(subarrayExtents.begin(), subarrayExtents.end());
  // the starts are local indices
  std::vector<int> starts(subarrayStarts.begin(), subarrayStarts.end());

  // create subarray view on data
  MPI_Datatype mysubarray;
  MPI_Type_create_subarray(static_cast<int>(this->getDimension()), sizes.data(), subsizes.data(),
                           starts.data(), MPI_ORDER_FORTRAN, this->getMPIDatatype(), &mysubarray);
  MPI_Type_commit(&mysubarray);
  return mysubarray;
}

template <typename FG_ELEMENT>
std::vector<MPI_Datatype> DistributedFullGrid<FG_ELEMENT>::getUpwardSubarrays() {
  // initialize upwardSubarrays_ only once
  if (upwardSubarrays_.size() == 0) {
    upwardSubarrays_.resize(this->getDimension());
    for (DimType d = 0; d < this->getDimension(); ++d) {
      upwardSubarrays_[d] = getUpwardSubarray(d);
    }
  }
  return upwardSubarrays_;
}

template <typename FG_ELEMENT>
MPI_Datatype DistributedFullGrid<FG_ELEMENT>::getDownwardSubarray(DimType d) {
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
  std::vector<int> sizes(this->getLocalSizes().begin(), this->getLocalSizes().end());
  std::vector<int> subsizes(subarrayExtents.begin(), subarrayExtents.end());
  // the starts are local indices
  std::vector<int> starts(subarrayStarts.begin(), subarrayStarts.end());

  // create subarray view on data
  MPI_Datatype mysubarray;
  MPI_Type_create_subarray(static_cast<int>(this->getDimension()), sizes.data(), subsizes.data(),
                           starts.data(), MPI_ORDER_FORTRAN, this->getMPIDatatype(), &mysubarray);
  MPI_Type_commit(&mysubarray);
  return mysubarray;
}

template <typename FG_ELEMENT>
std::vector<MPI_Datatype> DistributedFullGrid<FG_ELEMENT>::getDownwardSubarrays() {
  // initialize downwardSubarrays_ only once
  if (downwardSubarrays_.size() == 0) {
    downwardSubarrays_.resize(this->getDimension());
    for (DimType d = 0; d < this->getDimension(); ++d) {
      downwardSubarrays_[d] = getDownwardSubarray(d);
    }
  }
  return downwardSubarrays_;
}

template <typename FG_ELEMENT>
MPICartesianUtils DistributedFullGrid<FG_ELEMENT>::cartesianUtils_;

// output operator
template <typename FG_ELEMENT>
inline std::ostream& operator<<(std::ostream& os, const DistributedFullGrid<FG_ELEMENT>& dfg) {
  dfg.print();

  return os;
}

/**
 * @brief OwningDistributedFullGrid : a DistributedFullGrid with ownership of the data
 *
 * subclass of DistributedFullGrid that allocates a data vector of type FG_ELEMENT
 */
template <typename FG_ELEMENT>
class OwningDistributedFullGrid : public DistributedFullGrid<FG_ELEMENT> {
 public:
  OwningDistributedFullGrid() = default;

  explicit OwningDistributedFullGrid(
      DimType dim, const LevelVector& levels, CommunicatorType const& comm,
      const std::vector<BoundaryType>& hasBdrPoints, const std::vector<int>& procs,
      bool forwardDecomposition = true,
      const std::vector<IndexVector>& decomposition = std::vector<IndexVector>())
      : DistributedFullGrid<FG_ELEMENT>(dim, levels, comm, hasBdrPoints, ownedDataVector_.data(),
                                        procs, forwardDecomposition, decomposition) {
    ownedDataVector_.resize(this->getNrLocalElements());
    this->setData(ownedDataVector_.data());
  }

  const std::vector<FG_ELEMENT>& getDataVector() const { return ownedDataVector_; }

  void setDataVector(std::vector<FG_ELEMENT>&& otherVector) {
    assert(otherVector.size() == ownedDataVector_.size());
    ownedDataVector_ = std::move(otherVector);
    this->setData(ownedDataVector_.data());
  }

  void swapDataVector(std::vector<FG_ELEMENT>& otherVector) {
    assert(otherVector.size() == ownedDataVector_.size());
    ownedDataVector_.swap(otherVector);
    this->setData(ownedDataVector_.data());
  }

 private:
  std::vector<FG_ELEMENT> ownedDataVector_{};
};

}  // namespace combigrid

#endif /* DISTRIBUTEDCOMBIFULLGRID_HPP_ */
