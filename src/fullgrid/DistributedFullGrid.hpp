#ifndef DISTRIBUTEDCOMBIFULLGRID_HPP_
#define DISTRIBUTEDCOMBIFULLGRID_HPP_

#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>
#include <string>

#include "fullgrid/FullGrid.hpp"
#include "fullgrid/SliceIterator.hpp"
#include "mpi/MPICartesianUtils.hpp"
#include "mpi/MPISystem.hpp"
#include "sparsegrid/DistributedSparseGridUniform.hpp"
#include "utils/IndexVector.hpp"
#include "utils/LevelSetUtils.hpp"
#include "utils/LevelVector.hpp"
#include "utils/PowerOfTwo.hpp"
#include "utils/Stats.hpp"
#include "utils/Types.hpp"

namespace combigrid {

/* a regular (equidistant) domain decompositioning for an even number of processes
 * leads to grid points on the (geometrical) process boundaries.
 * with the forwardDecomposition flag it can be decided if the grid points on
 * the process boundaries belong to the process on the right-hand side (true)
 * of the process boundary, or to the one on the left-hand side (false).
 */
static inline std::vector<IndexVector> getDefaultDecomposition(
    const IndexVector& globalNumPointsPerDimension,
    const std::vector<int>& cartesianProcsPerDimension, bool forwardDecomposition) {
  auto dim = static_cast<DimType>(globalNumPointsPerDimension.size());
  assert(cartesianProcsPerDimension.size() == dim);

  // create decomposition vectors
  std::vector<IndexVector> decomposition(dim);
  for (DimType i = 0; i < dim; ++i) {
    if (cartesianProcsPerDimension[i] > globalNumPointsPerDimension[i]) {
      throw std::invalid_argument(
          "change to coarser parallelization! currently not all processes can have points.");
    }
    IndexVector& llbnd = decomposition[i];
    llbnd.resize(cartesianProcsPerDimension[i]);

    for (int j = 0; j < cartesianProcsPerDimension[i]; ++j) {
      double tmp = static_cast<double>(globalNumPointsPerDimension[i]) * static_cast<double>(j) /
                   static_cast<double>(cartesianProcsPerDimension[i]);

      if (forwardDecomposition)
        llbnd[j] = static_cast<IndexType>(std::ceil(tmp));
      else
        llbnd[j] = static_cast<IndexType>(std::floor(tmp));
    }
  }
  return decomposition;
}

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
                      const std::vector<BoundaryType>& hasBdrPoints, const std::vector<int>& procs,
                      bool forwardDecomposition = true,
                      const std::vector<IndexVector>& decomposition = std::vector<IndexVector>(),
                      const BasisFunctionBasis* basis = NULL)
      : dim_(dim), levels_(levels), hasBoundaryPoints_(hasBdrPoints) {
    assert(levels_.size() == dim);
    assert(hasBoundaryPoints_.size() == dim);
    assert(procs.size() == dim);

    InitMPI(comm, procs);  // will also check grids per dim

    // set global num of elements and offsets
    nrElements_ = 1;
    offsets_.resize(dim_);
    nrPoints_.resize(dim_);
    gridSpacing_.resize(dim_);

    for (DimType j = 0; j < dim_; j++) {
      nrPoints_[j] = combigrid::getNumDofNodal(levels_[j], hasBoundaryPoints_[j]);
      offsets_[j] = nrElements_;
      nrElements_ = nrElements_ * nrPoints_[j];
      if (hasBoundaryPoints_[j] == 1) {
        assert(!decomposition.empty() || !forwardDecomposition);
      }
      assert(levels_[j] < 30);
      gridSpacing_[j] = oneOverPowOfTwo[levels_[j]];
    }

    if (decomposition.size() == 0) {
      decomposition_ = getDefaultDecomposition(
          nrPoints_, this->getCartesianUtils().getCartesianDimensions(), forwardDecomposition);
    } else {
      setDecomposition(decomposition);
    }
    myPartitionsLowerBounds_ = getLowerBounds(rank_);
    myPartitionsUpperBounds_ = getUpperBounds(rank_);

    // set local elements and local offsets
    nrLocalPoints_ = getUpperBounds() - getLowerBounds();

    nrLocalElements_ = 1;
    localOffsets_.resize(dim);
    // cf. https://en.wikipedia.org/wiki/Row-_and_column-major_order#Address_calculation_in_general
    // -> column-major order
    for (DimType j = 0; j < dim_; ++j) {
      localOffsets_[j] = nrLocalElements_;
      nrLocalElements_ *= nrLocalPoints_[j];
    }

    // in contrast to serial implementation we directly create the grid
    fullgridVector_.resize(nrLocalElements_);
  }

  // explicit DistributedFullGrid(const DistributedFullGrid& other) {
  //   this = DistributedFullGrid(other.getDimension(), other.getLevels(), other.getCommunicator(),
  //                              other.returnBoundaryFlags(), other.getParallelization(), true,
  //                              other.getDecomposition());
  //   for (IndexType li = 0; li < other.getNrLocalElements(); ++li) {
  //     this->getData()[li] = other.getData()[li];
  //   }
  // }

  // copy construction would need to duplicate communicator_
  DistributedFullGrid(const DistributedFullGrid& other) = delete;
  DistributedFullGrid& operator=(const DistributedFullGrid&) = delete;
  DistributedFullGrid(DistributedFullGrid&& other) = default;
  DistributedFullGrid& operator=(DistributedFullGrid&& other) = default;

  virtual ~DistributedFullGrid() {
    // MPI_Comm_free(&communicator_);
    for (size_t i = 0; i < upwardSubarrays_.size(); ++i) {
      MPI_Type_free(&upwardSubarrays_[i]);
    }
    for (size_t i = 0; i < downwardSubarrays_.size(); ++i) {
      MPI_Type_free(&downwardSubarrays_[i]);
    }
  }


class SubarrayIterator : public SliceIterator<FG_ELEMENT> {
 public:
  SubarrayIterator(const MPI_Datatype& subarrayType, DistributedFullGrid& dfg) : SliceIterator<FG_ELEMENT>() {
    this->offsets_ = dfg.getLocalOffsets();
    this->dataPointer_ = &(dfg.getElementVector());

    // read out the subarray informaiton
    MPI_Datatype tmp;
    DimType dim = dfg.getDimension();
    std::vector<int> MPITypeIntegers(3 * dim + 2);
    auto success = MPI_Type_get_contents(subarrayType, static_cast<int>(MPITypeIntegers.size()), 0,
                                         1, MPITypeIntegers.data(), nullptr, &tmp);
    assert(success == MPI_SUCCESS);
    // check datatype
    assert(tmp == dfg.getMPIDatatype());
    // check order
    assert(MPITypeIntegers.back() == MPI_ORDER_FORTRAN);

    auto sizes = std::vector<int>(MPITypeIntegers.begin() + 1, MPITypeIntegers.begin() + dim + 1);
    this->subsizes_ =
        std::vector<int>(MPITypeIntegers.begin() + dim + 1, MPITypeIntegers.begin() + 2 * dim + 1);
    this->starts_ = std::vector<int>(MPITypeIntegers.begin() + 2 * dim + 1,
                                     MPITypeIntegers.begin() + 3 * dim + 1);
    this->currentLocalIndex_ = this->starts_;

    assert(dim == this->getDimension());
    // check that sizes match
    for (DimType d = 0; d < this->getDimension() - 1; ++d) {
      assert(sizes[d] == static_cast<int>(dfg.getLocalSizes()[d]));
      assert(this->subsizes_[d] + this->starts_[d] <= sizes[d]);
    }
    assert(std::accumulate(sizes.begin(), sizes.end(), 1, std::multiplies<int>()) ==
           this->dataPointer_->size());

    this->setEndIndex();
    this->validateSizes();
  }
};
  FG_ELEMENT evalLocalIndexOn(const IndexVector& localIndex,
                              const std::vector<real>& coords) const {
    const auto& lowerBounds = this->getLowerBounds();

    // get product of 1D hat functions on coords
    const auto& h = getGridSpacing();
    real phi_c = 1.;  // value of product of basis function on coords
    for (DimType d = 0; d < dim_; ++d) {
      // get distance between coords and point
      // this should be the same as
      // coordDistance = this->getCoordsLocal(localIndex) - coords;
      // cf. getLowerBoundsCoords()
      auto coordDistance =
          static_cast<double>(lowerBounds[d] + (hasBoundaryPoints_[d] > 0 ? 0 : 1) +
                              localIndex[d]) *
              h[d] -
          coords[d];
#ifndef NDEBUG
      if (std::abs(coordDistance) > h[d]) {
        std::cout << "assert bounds " << coordDistance << coords << h << static_cast<int>(d)
                  << localIndex << std::endl;
        assert(false &&
               "should only be called for coordinates within the support of this point's basis "
               "function");
      }
#endif  // ndef NDEBUG
      phi_c *= 1. - std::abs(coordDistance / h[d]);
    }

    auto localLinearIndex = getLocalLinearIndex(localIndex);
    // std::cout << "coords " <<  localIndex << coords << localLinearIndex << h << std::endl;
    // std::cout << "phi_c " << phi_c << this->getElementVector()[localLinearIndex] << std::endl;
    assert(phi_c >= 0.);
    return phi_c * this->getElementVector()[localLinearIndex];
  }

  /**
   * @brief recursive call to evaluate all neighbor points' contributions to the coordinate (on this
   * part of the grid)
   *
   * @param localIndex the (in-all-dimensions lower) neighbor of coords
   * @param dim the current dimension to split on (start with 0)
   * @param coords the coordinate to interpolate on
   * @return FG_ELEMENT the interpolated value at coords
   */
  FG_ELEMENT evalMultiindexRecursively(const IndexVector& localIndex, DimType dim,
                                       const std::vector<real>& coords) const {
    assert(!(dim > this->getDimension()));
    if (dim == this->getDimension()) {
      return evalLocalIndexOn(localIndex, coords);
    } else {
      FG_ELEMENT sum = 0.;
      auto lastIndexInDim = this->getLocalSizes()[dim] - 1;
      if (localIndex[dim] >= 0 && localIndex[dim] <= lastIndexInDim) {
        sum += evalMultiindexRecursively(localIndex, static_cast<DimType>(dim + 1), coords);
      } else if (localIndex[dim] > lastIndexInDim &&
                 !(this->hasBoundaryPoints_[dim] == 1 &&
                   this->getCartesianUtils().isOnLowerBoundaryInDimension(dim))) {
        throw std::runtime_error("localIndex[dim] > lastIndexInDim");
      }

      auto localIndexDimPlusOne = localIndex[dim] + 1;
      if (localIndexDimPlusOne >= 0) {
        if (this->hasBoundaryPoints_[dim] == 1 &&
            this->getCartesianUtils().isOnLowerBoundaryInDimension(dim) &&
            localIndexDimPlusOne == this->getGlobalSizes()[dim]) {
          // assume periodicity
          // if we are at the end of the dimension, wrap around
          localIndexDimPlusOne = 0;
          auto secondCoords = coords;
          secondCoords[dim] -= 1.;
          // these copies are used to find all the neighbor points of the current coordinate
          auto neighborIndex = localIndex;
          neighborIndex[dim] = localIndexDimPlusOne;
          sum +=
              evalMultiindexRecursively(neighborIndex, static_cast<DimType>(dim + 1), secondCoords);
        } else if (localIndexDimPlusOne <= lastIndexInDim) {
          // these copies are used to find all the neighbor points of the current coordinate
          auto neighborIndex = localIndex;
          neighborIndex[dim] = localIndexDimPlusOne;
          sum += evalMultiindexRecursively(neighborIndex, static_cast<DimType>(dim + 1),
                                            coords);
        }
      } else {
        throw std::runtime_error("localIndexDimPlusOne < 0");
      }
      return sum;
    }
  }

  /** evaluates the full grid on the specified coordinates
   * @param coords ND coordinates on the unit square [0,1]^D*/
  FG_ELEMENT evalLocal(const std::vector<real>& coords) const {
    FG_ELEMENT value;
    evalLocal(coords, value);
    return value;
  }

  void evalLocal(const std::vector<real>& coords, FG_ELEMENT& value) const {
    assert(coords.size() == this->getDimension());

    // get the lowest-index point of the points
    // whose basis functions contribute to the interpolated value
    auto lowerCoords = getLowerBoundsCoords();
    const auto& h = getGridSpacing();
    static IndexVector localIndexLowerNonzeroNeighborPoint;
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
      localIndexLowerNonzeroNeighborPoint[d] =
          static_cast<IndexType>(std::floor((coords[d] - lowerCoords[d]) / h[d]));

      // check if we even need to recursively evaluate on this process
      if (localIndexLowerNonzeroNeighborPoint[d] < -1) {
        // index too small
        value = 0.;
        return;
      } else if ((coords[d] >= 1.0 - h[d] && coords[d] <= 1.0) &&
                 this->hasBoundaryPoints_[d] == 1 &&
                 this->getCartesianUtils().isOnLowerBoundaryInDimension(d)) {
        // if we have periodic boundary and this process is at the lower end of the dimension d
        // we need the periodic coordinate => do nothing
      } else if (localIndexLowerNonzeroNeighborPoint[d] > this->getLocalSizes()[d] - 1) {
        // index too high
        value = 0.;
        return;
      }
    }
    // evaluate at those points and sum up according to the basis function
    // needs to be recursive in order to be dimensionally adaptive
    value = evalMultiindexRecursively(localIndexLowerNonzeroNeighborPoint, 0, coords);
  }

  /** evaluates the full grid on the specified coordinates
   * @param interpolationCoords vector of ND coordinates on the unit square [0,1]^D*/
  std::vector<FG_ELEMENT> getInterpolatedValues(
      const std::vector<std::vector<real>>& interpolationCoords) const {
    auto numValues = interpolationCoords.size();
    std::vector<FG_ELEMENT> values;
    values.resize(numValues);
    for (size_t i = 0; i < numValues; ++i) {
      this->evalLocal(interpolationCoords[i], values[i]);
    }
    MPI_Allreduce(MPI_IN_PLACE, values.data(), static_cast<int>(numValues), this->getMPIDatatype(),
                  MPI_SUM, this->getCommunicator());
    return values;
  }

  /** return the coordinates on the unit square corresponding to global idx
   * @param globalIndex [IN] global linear index of the element i
   * @param coords [OUT] the vector must be resized already */
  inline void getCoordsGlobal(IndexType globalLinearIndex, std::vector<real>& coords) const {
    IndexType ind = 0;
    IndexType tmp_add = 0;

    coords.resize(dim_);

    for (DimType j = 0; j < dim_; j++) {
      ind = globalLinearIndex % nrPoints_[j];
      globalLinearIndex = globalLinearIndex / nrPoints_[j];
      // set the coordinate based on if we have boundary points
      tmp_add = (hasBoundaryPoints_[j] > 0) ? (0) : (1);
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
      startindex = (hasBoundaryPoints_[k] > 0) ? 0 : 1;
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

    for (auto i_shifted = dim_; i_shifted > 0; --i_shifted) {
      auto dim_i = static_cast<DimType>(i_shifted - 1);
      globAxisIndex[dim_i] = tmp / (this->getOffset(dim_i));
      tmp = tmp % this->getOffset(dim_i);
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

    for (int i = static_cast<int>(dim_) - 1; i >= 0; --i) {
      const auto a = tmp / localOffsets_[i];
      const auto t = tmp % localOffsets_[i];

      locAxisIndex[i] = a;
      tmp = t;
    }
  }

  /** the local vector index corresponding to a global vector index
   if global index vector not contained in local domain false will be
   returned */
  inline bool getLocalVectorIndex(const IndexVector& globAxisIndex,
                                  IndexVector& locAxisIndex) const {
    assert(globAxisIndex.size() == dim_);

    if (globAxisIndex >= getLowerBounds() && globAxisIndex < getUpperBounds()) {
      locAxisIndex.assign(globAxisIndex.begin(), globAxisIndex.end());
      std::transform(locAxisIndex.begin(), locAxisIndex.end(), this->getLowerBounds().begin(),
                     locAxisIndex.begin(), std::minus<IndexType>());
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
    static IndexVector locAxisIndex(dim_);
    locAxisIndex.resize(dim_);
    getLocalVectorIndex(locLinIndex, locAxisIndex);

    // convert to global vector index
    static IndexVector globAxisIndex(dim_);
    globAxisIndex.resize(dim_);
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
    static IndexVector globAxisIndex(dim_);
    globAxisIndex.resize(dim_);
    getGlobalVectorIndex(globLinIndex, globAxisIndex);

    // convert to local vector index
    static IndexVector locAxisIndex(dim_);
    locAxisIndex.resize(dim_);

    if (getLocalVectorIndex(globAxisIndex, locAxisIndex)) {
      // convert to local linear index
      return getLocalLinearIndex(locAxisIndex);
    } else {
      return -1;
    }
  }

  inline bool isGlobalIndexHere(IndexType globLinIndex) const {
    return getLocalLinearIndex(globLinIndex) > -1;
  }

  inline bool isGlobalIndexHere(IndexVector globLinIndex) const {
    for (DimType d = 0; d < this->getDimension(); ++d) {
      if (globLinIndex[d] < this->getLowerBounds()[d] ||
          globLinIndex[d] >= this->getUpperBounds()[d]) {
        return false;
      }
    }
    return true;
  }

  /** returns the dimension of the full grid */
  inline DimType getDimension() const { return dim_; }

  /** the getters for the full grid vector */
  inline std::vector<FG_ELEMENT>& getElementVector() { return fullgridVector_; }

  inline const std::vector<FG_ELEMENT>& getElementVector() const { return fullgridVector_; }

  /** return the offset in the full grid vector of the dimension */
  inline IndexType getOffset(DimType i) const { return offsets_[i]; }

  inline const IndexVector& getOffsets() const { return offsets_; }

  inline const IndexVector& getLocalOffsets() const {
    assert(localOffsets_.size() == dim_);
    assert(localOffsets_[0] == 1);
    return localOffsets_;
  }

  /** return the level vector */
  inline const LevelVector& getLevels() const { return levels_; }

  /** returns the grid spacing (sometimes called h) */
  inline const std::vector<double>& getGridSpacing() const { return gridSpacing_; }

  inline std::vector<double> getInverseGridSpacing() const {
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
                   std::bind1st(std::divides<double>(), 1));
    for (DimType j = 0; j < dim_; j++) {
      assert(std::abs(oneOverHByDivision[j] < oneOverH[j]) < 1e-10);
      assert(oneOverHByDivision[j] == oneOverH[j]);
    }
#endif
    return oneOverH;
  }

  double getInnerNodalBasisFunctionIntegral() const {
    const auto& h = this->getGridSpacing();
    return std::accumulate(h.begin(), h.end(), 1., std::multiplies<double>());
  }

  /** returns the number of elements in the full grid */
  inline IndexType getNrElements() const { return nrElements_; }

  /** number of elements in the local partition */
  inline IndexType getNrLocalElements() const { return nrLocalElements_; }

  /** number of points per dimension i */
  inline IndexType length(DimType i) const { return nrPoints_[i]; }

  /** vector of flags to show if the dimension has boundary points*/
  inline const std::vector<BoundaryType>& returnBoundaryFlags() const { return hasBoundaryPoints_; }

  inline void setZero() {
    std::fill(this->getElementVector().begin(), this->getElementVector().end(), 0.);
  }

  inline FG_ELEMENT* getData() { return &fullgridVector_[0]; }

  inline const FG_ELEMENT* getData() const { return &fullgridVector_[0]; }

  /** MPI Communicator*/
  inline CommunicatorType getCommunicator() const { return communicator_; }

  inline RankType getRank() const { return rank_; }

  inline int getCommunicatorSize() const { return size_; }

  /** lower Bounds of this process */
  inline const IndexVector& getLowerBounds() const {
    // return getLowerBounds(rank_);
    return myPartitionsLowerBounds_;
  }

  /** lower bounds of rank r */
  inline IndexVector getLowerBounds(RankType r) const {
    assert(r >= 0 && r < size_);
    // get coords of r in cart comm
    std::vector<int> coords(dim_);
    IndexVector lowerBounds(dim_);
    cartesianUtils_.getPartitionCoordsOfRank(r, coords);

    for (DimType i = 0; i < dim_; ++i) {
      lowerBounds[i] = getDecomposition()[i][coords[i]];
    }
    return lowerBounds;
  }

  inline IndexType getLowestGlobalIndexOfRank(RankType r) const {
    auto lowestIndex = getGlobalLinearIndex(getLowerBounds(r));
    return lowestIndex;
  }

  /** coordinates of this process' lower bounds */
  inline std::vector<real> getLowerBoundsCoords() const {
    const auto& lowerBounds = this->getLowerBounds();
    static std::vector<real> coords(dim_);
    coords.resize(dim_);
    for (DimType i = 0; i < dim_; ++i) {
      coords[i] = static_cast<double>(lowerBounds[i] + (hasBoundaryPoints_[i] > 0 ? 0 : 1)) *
                  this->getGridSpacing()[i];
    }
    return coords;
  }

  /** coordinates of rank r's lower bounds */
  inline std::vector<real> getLowerBoundsCoords(RankType r) const {
    assert(r >= 0 && r < size_);
    IndexType lbli = getLowestGlobalIndexOfRank(r);
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
    cartesianUtils_.getPartitionCoordsOfRank(r, coords);

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
    for (size_t j = 0; j < dim_; ++j)
      decompositionCoords[j].resize(this->getCartesianUtils().getCartesianDimensions()[j]);
    assert(false && "this is pretty much untested, please add test before using this");

    for (RankType r = 0; r < size_; ++r) {
      // get coords of r in cart comm
      std::vector<int> coords(dim_);
      cartesianUtils_.getPartitionCoordsOfRank(r, coords);

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
      tmp_add = (hasBoundaryPoints_[j] > 0) ? (0) : (1);
      coords[j] = static_cast<double>(ubvi[j] + tmp_add) * getGridSpacing()[j];
    }
    return coords;
  }

  /* returns the neighboring process (in the sense that the neighbor has the same
   * partion coordinates in all other dimensions than d) in dimension d which
   * contains the point with the one-dimensional index idx1d
   */
  RankType getNeighbor1dFromAxisIndex(DimType dim, IndexType idx1d) const {
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

  /** Number of Grids in every dimension*/
  inline const std::vector<int>& getParallelization() const {
    return this->getCartesianUtils().getCartesianDimensions();
  }

  /** MPI Rank */
  inline int getMpiRank() const { return rank_; }

  /** MPI Size */
  inline int getMpiSize() const { return size_; }

  /** returns the 1d global index of the first point in the local domain
   *
   */
  inline IndexType getFirstGlobal1dIndex(DimType d) const { return getLowerBounds()[d]; }

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
    LevelType l = static_cast<LevelType>(lmax - __builtin_ctzl(static_cast<unsigned long>(idx1d)));
    return l;
  }

  // get 1d index of LeftPredecessor of a point
  // returns negative number if point has no left predecessor
  inline IndexType getLeftPredecessor(DimType d, IndexType idx1d) const {
    LevelType l = getLevel(d, idx1d);

    // boundary points
    if (l == 0) return -1;

    LevelType ldiff = static_cast<LevelType>(levels_[d] - l);
    IndexType lpidx = idx1d - combigrid::powerOfTwoByBitshift(ldiff);

    return lpidx;
  }

  inline IndexType getRightPredecessor(DimType d, IndexType idx1d) const {
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

  // get coordinates of the partition which contains the point specified
  // by the global index vector
  inline void getPartitionCoords(IndexVector& globalAxisIndex, std::vector<int>& partitionCoords) const {
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

        std::vector<int> csizes(sizes.begin(), sizes.end());
        std::vector<int> csubsizes(subsizes.begin(), subsizes.end());
        std::vector<int> cstarts(starts.begin(), starts.end());

        // create subarray view on data
        MPI_Datatype mysubarray;
        MPI_Type_create_subarray(static_cast<int>(this->getDimension()), &csizes[0], &csubsizes[0],
                                 &cstarts[0], MPI_ORDER_FORTRAN, this->getMPIDatatype(),
                                 &mysubarray);
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

  inline void getFGPointsOfSubspaceRecursive(DimType d, IndexType localLinearIndexSum,
                                             std::vector<IndexVector>& oneDIndices,
                                             std::vector<IndexType>& subspaceIndices) const {
    assert(d < dim_);
    assert(oneDIndices.size() == dim_);
    assert(!oneDIndices.empty());

    for (const auto idx : oneDIndices[d]) {
      auto updatedLocalIndexSum = localLinearIndexSum;
      updatedLocalIndexSum += localOffsets_[d] * idx;
      if (d > 0) {
        getFGPointsOfSubspaceRecursive(static_cast<DimType>(d - 1), updatedLocalIndexSum,
                                       oneDIndices, subspaceIndices);
      } else {
        subspaceIndices.emplace_back(updatedLocalIndexSum);
      }
    }
  }

  /**
   * @brief Get the indices of the points of the subspace on this partition
   *
   * @param l level of hierarchical subspace
   * @return the indices of points on this partition
   */
  inline std::vector<IndexType> getFGPointsOfSubspace(const LevelVector& l) const {
    IndexVector subspaceIndices;
    IndexType numPointsOfSubspace = 1;
    static auto oneDIndices = std::vector<IndexVector>(dim_);
    oneDIndices.resize(dim_);
    for (DimType d = 0; d < dim_; ++d) {
      if (l[d] > levels_[d]) {
        return subspaceIndices;
      }
      get1dIndicesLocal(d, l[d], oneDIndices[d]);
      numPointsOfSubspace *= oneDIndices[d].size();
    }
    if (numPointsOfSubspace > 0) {
      subspaceIndices.reserve(numPointsOfSubspace);

      IndexType localLinearIndexSum = 0;
      getFGPointsOfSubspaceRecursive(static_cast<DimType>(dim_ - 1), localLinearIndexSum,
                                     oneDIndices, subspaceIndices);
    }
    assert(static_cast<IndexType>(subspaceIndices.size()) == numPointsOfSubspace);
    return subspaceIndices;
  }

  /**
   * @brief extracts the (hopefully) hierarchical coefficients from dsg
   *        to the full grid's data structure
   *
   * @param dsg the DSG to extract from
   */
  void extractFromUniformSG(const DistributedSparseGridUniform<FG_ELEMENT>& dsg) {
    assert(dsg.isSubspaceDataCreated());

    // all the hierarchical subspaces contained in this full grid
    const auto downwardClosedSet = combigrid::getDownSet(levels_);

    // loop over all subspaces (-> somewhat linear access in the sg)
    typename DistributedSparseGridUniform<FG_ELEMENT>::SubspaceIndexType sIndex = 0;
    static IndexVector subspaceIndices;
    for (const auto& level : downwardClosedSet) {
      sIndex = dsg.getIndexInRange(level, sIndex);
      if (sIndex > -1 && dsg.getDataSize(sIndex) > 0) {
        auto sPointer = dsg.getData(sIndex);
        subspaceIndices = getFGPointsOfSubspace(level);
        for (const auto& fIndex : subspaceIndices) {
          fullgridVector_[fIndex] = *sPointer;
          ++sPointer;
        }
      }
    }
  }

  inline IndexType getStrideForThisLevel(LevelType l, DimType d) const {
    assert(d < this->getDimension());
    // special treatment for level 1 suspaces with boundary
    return (l == 1 && hasBoundaryPoints_[d] > 0)
               ? combigrid::powerOfTwoByBitshift(static_cast<LevelType>(levels_[d] - 1))
               : combigrid::powerOfTwoByBitshift(static_cast<LevelType>(levels_[d] - l + 1));
  }

  inline IndexType getLocalStartForThisLevel(LevelType l, DimType d, IndexType strideForThisLevel) const {
    const auto firstGlobal1dIdx = getFirstGlobal1dIndex(d);

    // get global offset to find indices of this level
    // this is the first global index that has level l in dimension d
    IndexType offsetForThisLevel;
    if (hasBoundaryPoints_[d] > 0) {
      if (l == 1) {
        offsetForThisLevel = 0;
      } else {
        // offsetForThisLevel = strideForThisLevel / 2;
        offsetForThisLevel =
            combigrid::powerOfTwoByBitshift(static_cast<LevelType>(levels_[d] - l));
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

  inline IndexType getNumPointsOnThisPartition(DimType d, IndexType localStart,
                                               IndexType strideForThisLevel) const {
    assert(d < this->getDimension());
    return (nrLocalPoints_[d] - 1 < localStart)
               ? 0
               : (nrLocalPoints_[d] - 1 - localStart) / strideForThisLevel + 1;
  }

  inline IndexType getNumPointsOnThisPartition(LevelType l, DimType d) const {
    assert(!(l > levels_[d]));
    const auto strideForThisLevel = getStrideForThisLevel(l, d);
    return getNumPointsOnThisPartition(d, getLocalStartForThisLevel(l, d, strideForThisLevel),
                                       strideForThisLevel);
  }

  /**
   * @brief get the local 1d indices that correspond to a given level
   *
   * @param d the dimension in which we want the indices for the level
   * @param l the level
   * @param oneDIndices a list of the local vector indices
   */
  inline void get1dIndicesLocal(DimType d, LevelType l, IndexVector& oneDIndices) const {
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

  /**
   * @brief Get the Lp Norm of the data on the dfg:
   *      norm = (sum_i(abs(x_i)^p)*integral(basis_i))^1/p
   *
   * @param p : the (polynomial) parameter, 0 is interpreted as maximum norm
   * @return real : the norm
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

      // double sumIntegral = 0;
      for (size_t i = 0; i < data.size(); ++i) {
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

    // we store our data in fortran notation, with the
    // first index in indexvectors being the first dimension.
    std::vector<int> csizes(sizes.begin(), sizes.end());
    std::vector<int> csubsizes(subsizes.begin(), subsizes.end());
    std::vector<int> cstarts(starts.begin(), starts.end());

    // create subarray view on data
    MPI_Datatype mysubarray;
    MPI_Type_create_subarray(static_cast<int>(getDimension()), &csizes[0], &csubsizes[0],
                             &cstarts[0], MPI_ORDER_FORTRAN, getMPIDatatype(), &mysubarray);
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
        MPI_File_write(fh, &dim, 1, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);

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

    // we store our data in fortran notation, with the
    // first index in indexvectors being the first dimension.
    std::vector<int> csizes(sizes.begin(), sizes.end());
    std::vector<int> csubsizes(subsizes.begin(), subsizes.end());
    std::vector<int> cstarts(starts.begin(), starts.end());

    // create subarray view on data
    MPI_Datatype mysubarray;
    MPI_Type_create_subarray(static_cast<int>(getDimension()), &csizes[0], &csubsizes[0],
                             &cstarts[0], MPI_ORDER_FORTRAN, getMPIDatatype(), &mysubarray);
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
      MPI_File_write(fh, header_string.data(), static_cast<int>(header_size), MPI_CHAR, MPI_STATUS_IGNORE);
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

  MPI_Datatype getUpwardSubarray(DimType d) {
    // do index calculations
    // set lower bounds of subarray
    IndexVector subarrayLowerBounds = this->getLowerBounds();
    IndexVector subarrayUpperBounds = this->getUpperBounds();
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

  std::vector<MPI_Datatype> getUpwardSubarrays() {
    // initialize upwardSubarrays_ only once
    if (upwardSubarrays_.size() == 0) {
      upwardSubarrays_.resize(this->getDimension());
      for (DimType d = 0; d < this->getDimension(); ++d) {
        upwardSubarrays_[d] = getUpwardSubarray(d);
      }
    }
    return upwardSubarrays_;
  }

  MPI_Datatype getDownwardSubarray(DimType d) {
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

  std::vector<MPI_Datatype> getDownwardSubarrays() {
    // initialize downwardSubarrays_ only once
    if (downwardSubarrays_.size() == 0) {
      downwardSubarrays_.resize(this->getDimension());
      for (DimType d = 0; d < this->getDimension(); ++d) {
        downwardSubarrays_[d] = getDownwardSubarray(d);
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

  void writeUpperBoundaryToLowerBoundary(DimType d) {
    assert(hasBoundaryPoints_[d] == 2);

    // create MPI datatypes
    auto downSubarray = getDownwardSubarray(d);
    auto upSubarray = getUpwardSubarray(d);

    // if I have the lowest neighbor (i. e. I am the highest rank), I need to send my highest layer in d to them,
    // if I have the highest neighbor (i. e. I am the lowest rank), I can receive it
    int lower, higher;
    getHighestAndLowestNeighbor(d, higher, lower);

    auto success =
        MPI_Sendrecv(this->getData(), 1, upSubarray, lower, TRANSFER_GHOST_LAYER_TAG,
                     this->getData(), 1, downSubarray, higher,
                     TRANSFER_GHOST_LAYER_TAG, this->getCommunicator(), MPI_STATUS_IGNORE);
    MPI_Type_free(&downSubarray);
    MPI_Type_free(&upSubarray);
  }

  void writeLowerBoundaryToUpperBoundary(DimType d) {
    assert(hasBoundaryPoints_[d] == 2);

    // create MPI datatypes
    auto downSubarray = getDownwardSubarray(d);
    auto upSubarray = getUpwardSubarray(d);

    // if I have the highest neighbor (i. e. I am the lowest rank), I need to send my lowest layer in d to them,
    // if I have the lowest neighbor (i. e. I am the highest rank), I can receive it
    int lower, higher;
    getHighestAndLowestNeighbor(d, higher, lower);

    // TODO asynchronous over d??
    auto success =
        MPI_Sendrecv(this->getData(), 1, downSubarray, higher, TRANSFER_GHOST_LAYER_TAG,
                     this->getData(), 1, upSubarray, lower,
                     TRANSFER_GHOST_LAYER_TAG, this->getCommunicator(), MPI_STATUS_IGNORE);
    assert(success == MPI_SUCCESS);
    MPI_Type_free(&downSubarray);
    MPI_Type_free(&upSubarray);
  }

  void exchangeBoundaryLayers(DimType d, std::vector<FG_ELEMENT>& recvbufferFromUp,
                              std::vector<FG_ELEMENT>& recvbufferFromDown) {
    assert(hasBoundaryPoints_[d] == 2);
    auto subarrayExtents = this->getLocalSizes();
    subarrayExtents[d] = 1;
    auto numElements = std::accumulate(subarrayExtents.begin(), subarrayExtents.end(), 1,
                                       std::multiplies<IndexType>());

    // create MPI datatypes
    auto downSubarray = getDownwardSubarray(d);
    auto upSubarray = getUpwardSubarray(d);

    // if I have the highest neighbor (i. e. I am the lowest rank), I need to send my lowest layer
    // in d to them, if I have the lowest neighbor (i. e. I am the highest rank), I can receive it
    int lower, higher;
    getHighestAndLowestNeighbor(d, higher, lower);
    if (lower != MPI_PROC_NULL) {
      recvbufferFromDown.resize(numElements);
    } else {
      recvbufferFromDown.resize(0);
    }
    if (higher != MPI_PROC_NULL) {
      recvbufferFromUp.resize(numElements);
    } else {
      recvbufferFromUp.resize(0);
    }

    // TODO asynchronous??
    // send lower boundary values
    auto success =
        MPI_Sendrecv(this->getData(), 1, downSubarray, higher, TRANSFER_GHOST_LAYER_TAG,
                     recvbufferFromDown.data(), static_cast<int>(recvbufferFromDown.size()),
                     this->getMPIDatatype(), lower, TRANSFER_GHOST_LAYER_TAG,
                     this->getCommunicator(), MPI_STATUS_IGNORE);
    assert(success == MPI_SUCCESS);
    // send upper boundary values
    success = MPI_Sendrecv(this->getData(), 1, upSubarray, lower, TRANSFER_GHOST_LAYER_TAG,
                           recvbufferFromUp.data(), static_cast<int>(recvbufferFromUp.size()),
                           this->getMPIDatatype(), higher, TRANSFER_GHOST_LAYER_TAG,
                           this->getCommunicator(), MPI_STATUS_IGNORE);
    assert(success == MPI_SUCCESS);
    MPI_Type_free(&downSubarray);
    MPI_Type_free(&upSubarray);
  }

  // non-RVO dependent version of ghost layer exchange
  void exchangeGhostLayerUpward(DimType d, std::vector<int>& subarrayExtents,
                                std::vector<FG_ELEMENT>& recvbuffer) {
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
    MPI_Cart_shift( this->getCommunicator(), static_cast<int>(d_reverse), 1, &lower, &higher );

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
      std::fill(subarrayExtents.begin(), subarrayExtents.end(), 0);
    }
    recvbuffer.resize(numElements);
    std::fill(recvbuffer.begin(), recvbuffer.end(), 0.);

    // TODO asynchronous over d??
    auto success =
        MPI_Sendrecv(this->getData(), 1, subarray, higher, TRANSFER_GHOST_LAYER_TAG,
                     recvbuffer.data(), numElements, this->getMPIDatatype(), lower,
                     TRANSFER_GHOST_LAYER_TAG, this->getCommunicator(), MPI_STATUS_IGNORE);
    assert(success == MPI_SUCCESS);
    MPI_Type_free(&subarray);
  }

  std::vector<FG_ELEMENT> exchangeGhostLayerUpward(DimType d, std::vector<int>& subarrayExtents) {
    std::vector<FG_ELEMENT> recvbuffer{};
    exchangeGhostLayerUpward(d, subarrayExtents, recvbuffer);
    return recvbuffer;
  }

  void averageBoundaryValues(DimType d) {
    std::vector<FG_ELEMENT> recvbufferFromUp, recvbufferFromDown;
    // exchanges just like in the ghost layer exchange, but across boundary
    exchangeBoundaryLayers(d, recvbufferFromUp, recvbufferFromDown);
    MPI_Barrier(this->getCommunicator());
    // get MPI datatypes
    auto downSubarrays = getDownwardSubarrays();
    auto upSubarrays = getUpwardSubarrays();

    // if both lower and higher are set, it's only me and I can just average
    // not strictly necessary but good to test the iterators
    if (recvbufferFromUp.size() > 0 && recvbufferFromDown.size() > 0) {
      SubarrayIterator downIt(downSubarrays[d], *this);
      SubarrayIterator upIt(upSubarrays[d], *this);
      assert(downIt.size() == upIt.size());
#ifndef NDEBUG
      auto initialIndexDiff = upIt.getIndex() - downIt.getIndex();
      auto initialUpIndex = upIt.getIndex();
      int ctr = 0;
#endif  // ndef NDEBUG
      while (!downIt.isAtEnd()) {
        assert(initialIndexDiff == upIt.getIndex() - downIt.getIndex());
        assert(upIt.getIndex() >= initialUpIndex);
        auto avg = 0.5 * ((*upIt) + (*downIt));
        *upIt = avg;
        *downIt = avg;
        ++downIt;
        ++upIt;
        assert(++ctr);
      }
      assert(upIt.isAtEnd());
      assert(ctr == upIt.size());
    } else if (recvbufferFromUp.size() > 0) {
      SubarrayIterator downIt(downSubarrays[d], *this);
      auto upIt = recvbufferFromUp.begin();
      while (!downIt.isAtEnd()) {
        auto avg = 0.5 * ((*upIt) + (*downIt));
        *downIt = avg;
        ++downIt;
        ++upIt;
      }
      assert(upIt == recvbufferFromUp.end());
    } else if (recvbufferFromDown.size() > 0) {
      SubarrayIterator upIt(upSubarrays[d], *this);
      auto downIt = recvbufferFromDown.begin();
      while (!upIt.isAtEnd()) {
        auto avg = 0.5 * ((*upIt) + (*downIt));
        *upIt = avg;
        ++downIt;
        ++upIt;
      }
      assert(downIt == recvbufferFromDown.end());
    }
  }

  void averageBoundaryValues() {
    for (DimType d = 0; d < this->getDimension(); ++d) {
      if (this->hasBoundaryPoints_[d] == 2) {
        this->averageBoundaryValues(d);
      }
    }
  }

  /**
   * @brief check if given globalLinearIndex is on the boundary of this DistributedFullGrid
   */
  std::vector<bool> isGlobalLinearIndexOnBoundary(IndexType globalLinearIndex) const {
    // this could likely be done way more efficiently, but it's
    // currently not in any performance critical spot
    std::vector<bool> isOnBoundary(this->getDimension(), false);

    // convert to global vector index
    IndexVector globalAxisIndex(dim_);
    getGlobalVectorIndex(globalLinearIndex, globalAxisIndex);
    for (DimType d = 0; d < this->getDimension(); ++d) {
      if (this->returnBoundaryFlags()[d] == 2) {
        if (globalAxisIndex[d] == 0 || globalAxisIndex[d] == this->length(d) - 1) {
          isOnBoundary[d] = true;
        }
      }
    }
    return isOnBoundary;
  }

  /**
   * @brief check if given localLinearIndex is on the boundary of this DistributedFullGrid
   */
  std::vector<bool> isLocalLinearIndexOnBoundary(IndexType localLinearIndex) const {
    return isGlobalLinearIndexOnBoundary(getGlobalLinearIndex(localLinearIndex));
  }

  /**
   * @brief recursive helper function for getCornersGlobal*Indices
   *        (otherwise, it can't be dimension-independent)
   *
   */
  std::vector<IndexVector> getCornersGlobalVectorIndicesRecursive(
      std::vector<IndexVector> indicesSoFar, DimType dim) const {
    if (dim < this->getDimension()) {
      std::vector<IndexVector> newIndicesSoFar{};
      for (const auto& indexVec : indicesSoFar) {
        newIndicesSoFar.push_back(indexVec);
        newIndicesSoFar.back().push_back(0);
        newIndicesSoFar.push_back(indexVec);
        newIndicesSoFar.back().push_back(this->length(dim) - 1);
      }
      assert(newIndicesSoFar.size() == 2 * indicesSoFar.size());
      return getCornersGlobalVectorIndicesRecursive(newIndicesSoFar, static_cast<DimType>(dim + 1));
    } else {
      return indicesSoFar;
    }
  }

  /**
   * @brief get a vector containing the global vector indices of the 2^d corners of this dfg
   *
   */
  std::vector<IndexVector> getCornersGlobalVectorIndices() const {
    auto emptyVectorInVector = std::vector<IndexVector>(1);
    assert(emptyVectorInVector.size() == 1);
    assert(emptyVectorInVector[0].size() == 0);
    std::vector<IndexVector> cornersVectors =
        getCornersGlobalVectorIndicesRecursive(emptyVectorInVector, 0);
    assert(cornersVectors.size() == static_cast<size_t>(powerOfTwo[this->getDimension()]));
    return cornersVectors;
  }

  /**
   * @brief get a vector containing the global linear indices of the 2^d corners of this dfg
   *
   */
  std::vector<IndexType> getCornersGlobalLinearIndices() const {
    std::vector<IndexVector> cornersVectors = getCornersGlobalVectorIndices();
    std::vector<IndexType> corners;
    for (const auto& corner : cornersVectors) {
      corners.push_back(this->getGlobalLinearIndex(corner));
    }
    assert(corners.size() == powerOfTwo[this->getDimension()]);
    return corners;
  }

  std::vector<FG_ELEMENT> getCornersValues() const {
    std::vector<FG_ELEMENT> values(powerOfTwo[this->getDimension()]);
    auto corners = getCornersGlobalVectorIndices();
    for (size_t cornerNo = 0; cornerNo < corners.size(); ++cornerNo) {
      if (this->isGlobalIndexHere(corners[cornerNo])) {
        // convert to local vector index, then to linear index
        IndexVector locAxisIndex(this->getDimension());
        bool present = getLocalVectorIndex(corners[cornerNo], locAxisIndex);
        assert(present);
        auto index = getLocalLinearIndex(locAxisIndex);
        values[cornerNo] = this->getElementVector()[index];
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, values.data(), static_cast<int>(values.size()),
                  this->getMPIDatatype(), MPI_SUM, this->getCommunicator());
    return values;
  }

  const MPICartesianUtils& getCartesianUtils() const {
    return cartesianUtils_;
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

  /** the grid spacing h for each dimension */
  std::vector<double> gridSpacing_;

  /** number of points per axis boundary included*/
  IndexVector nrPoints_;

  /** flag to show if the dimension has boundary points*/
  std::vector<BoundaryType> hasBoundaryPoints_;

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

  /** utility to get info about cartesian communicator  */
  static MPICartesianUtils cartesianUtils_;

  /** mpi rank */
  RankType rank_;

  /** mpi size */
  int size_;

  // the MPI Datatypes representing the boundary layers of the MPI processes' subgrid
  std::vector<MPI_Datatype> downwardSubarrays_;
  std::vector<MPI_Datatype> upwardSubarrays_;

  /** number of local (in this grid cell) points per axis*/
  IndexVector nrLocalPoints_;

  /**
   * @brief sets the MPI-related members rank_, size_, communicator_, and cartesianUtils_
   *
   * @param comm the communicator to use (assumed to be cartesian)
   * @param procs the desired partition (for sanity checking)
   */
  void InitMPI(MPI_Comm comm, const std::vector<int>& procs) {
    MPI_Comm_rank(comm, &rank_);
    MPI_Comm_size(comm, &size_);

#ifndef NDEBUG
    auto numSubgrids =
        std::accumulate(procs.begin(), procs.end(), 1, std::multiplies<int>());

    ASSERT(size_ == static_cast<int>(numSubgrids),
           " size_: " << size_ << " numSubgrids: " << static_cast<int>(numSubgrids));
    assert(size_ == static_cast<int>(numSubgrids));
#endif

    // check if communicator is already cartesian
    int status;
    MPI_Topo_test(comm, &status);

    if (status == MPI_CART) {
#ifndef NDEBUG
      // check if process grid of comm uses the required ordering
      auto maxdims = static_cast<int>(procs.size());
      std::vector<int> cartdims(maxdims), periods(maxdims), coords(maxdims);
      MPI_Cart_get(comm, static_cast<int>(maxdims), &cartdims[0], &periods[0], &coords[0]);
      // important: note reverse ordering of dims!
      std::vector<int> dims(procs.size());
      if (reverseOrderingDFGPartitions) {
        dims.assign(procs.rbegin(), procs.rend());
      } else {
        dims.assign(procs.begin(), procs.end());
      }
      ASSERT(cartdims == dims, " cartdims: " << cartdims << " dims: " << dims);
      assert(cartdims == dims);
      // not sure if this should be asserted -> for now, check only in exchangeGhostLayer...
      // for (int d = 0; d < maxdims; ++d) {
      //   if (hasBoundaryPoints_[d] == 1) {
      //     // assert periodicity
      //     assert(periods[d] == 1);
      //   } else {
      //     // assert no periodicity
      //     assert(periods[d] == 0);
      //   }
      // }
#endif

      // MPI_Comm_dup(comm, &communicator_);
      // cf. https://www.researchgate.net/publication/220439585_MPI_on_millions_of_cores
      // "Figure 3 shows the memory consumption in all these cases after 32 calls to MPI Comm dup"
      communicator_ = comm;

    } else {
      // important: note reverse ordering of dims!
      std::vector<int> dims(procs.size());
      if (reverseOrderingDFGPartitions) {
        dims.assign(procs.rbegin(), procs.rend());
      } else {
        dims.assign(procs.begin(), procs.end());
      }
      // todo mh: think whether periodic bc will be useful
      std::vector<int> periods(dim_, 0);
      int reorder = 0;
      MPI_Cart_create(comm, static_cast<int>(dim_), &dims[0], &periods[0], reorder, &communicator_);
      throw std::runtime_error("Currently testing to not duplicate communicator (to save memory), \
                          if you do want to use this code please take care that \
                          MPI_Comm_free(&communicator_) will be called at some point");
    }
    {
      if (communicator_ != cartesianUtils_.getComm()) {
        assert(uniformDecomposition);
        cartesianUtils_ = MPICartesianUtils(communicator_);
      }
    }
  }

  void setDecomposition(const std::vector<IndexVector>& decomposition) {
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
#endif // not def NDEBUG

    decomposition_ = decomposition;
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

template <typename FG_ELEMENT>
MPICartesianUtils DistributedFullGrid<FG_ELEMENT>::cartesianUtils_;

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

static inline std::vector<IndexVector> downsampleDecomposition(
    const std::vector<IndexVector> decomposition, const LevelVector& referenceLevel,
    const LevelVector& newLevel, const std::vector<BoundaryType>& boundary) {
  auto newDecomposition = decomposition;
  if (decomposition.size() > 0) {
    for (DimType d = 0 ; d < static_cast<DimType>(referenceLevel.size()); ++ d) {
      // for now, assume that we never want to interpolate on a level finer than referenceLevel
      assert(referenceLevel[d] >= newLevel[d]);
      auto levelDiff = referenceLevel[d] - newLevel[d];
      auto stepFactor = oneOverPowOfTwo[levelDiff];
      if(boundary[d] > 0) {
        // all levels contain the boundary points -> point 0 is the same
        for (auto& dec: newDecomposition[d]) {
          dec = static_cast<IndexType>(std::ceil(static_cast<double>(dec)*stepFactor));
        }
      } else {
        // all levels do not contain the boundary points -> mid point is the same
        auto leftProtrusion = powerOfTwo[levelDiff] - 1;
        for (auto& dec: newDecomposition[d]) {
          // same as before, but subtract the "left" protrusion on the finer level
          dec = static_cast<IndexType>(std::max(0., std::ceil(static_cast<double>(dec - leftProtrusion)*stepFactor)));
        }
      }
    }
  }
  return newDecomposition;
}

}  // namespace combigrid

#endif /* DISTRIBUTEDCOMBIFULLGRID_HPP_ */
