#ifndef SRC_SGPP_COMBIGRID_SPARSEGRID_DISTRIBUTEDSPARSEGRIDUNIFORM_HPP_
#define SRC_SGPP_COMBIGRID_SPARSEGRID_DISTRIBUTEDSPARSEGRIDUNIFORM_HPP_
#include "sparsegrid/AnyDistributedSparseGrid.hpp"
#include "utils/LevelSetUtils.hpp"
#include "utils/Types.hpp"

#include <boost/iterator/counting_iterator.hpp>
#include <cassert>
#include <numeric>
#include <set>
#include <vector>


namespace combigrid {
// forward declarations
template <typename FG_ELEMENT>
class DistributedFullGrid;

template <typename FG_ELEMENT>
class DistributedSparseGridUniform;

template <typename FG_ELEMENT>
class DistributedSparseGridDataContainer {
 public:
  using SubspaceIndexType = AnyDistributedSparseGrid::SubspaceIndexType;

  explicit DistributedSparseGridDataContainer(DistributedSparseGridUniform<FG_ELEMENT>& dsgu)
      : dsgu_(dsgu) {
    subspaces_.resize(dsgu_.getNumSubspaces(), nullptr);
    kahanDataBegin_.resize(dsgu_.getNumSubspaces(), nullptr);
  }

  DistributedSparseGridDataContainer(const DistributedSparseGridDataContainer&) = delete;
  DistributedSparseGridDataContainer& operator=(const DistributedSparseGridDataContainer&) = delete;
  DistributedSparseGridDataContainer(DistributedSparseGridDataContainer&&) = default;
  DistributedSparseGridDataContainer& operator=(DistributedSparseGridDataContainer&&) = default;

  // allocates memory for subspace data and sets pointers to subspaces
  void createSubspaceData() {
    if (subspacesWithData_.empty()) {
      // create data for all subspaces from "iota"
      subspacesWithData_ = std::set<SubspaceIndexType>{
          boost::counting_iterator<SubspaceIndexType>(0),
          boost::counting_iterator<SubspaceIndexType>(dsgu_.getNumSubspaces())};
    }
    size_t numDataPoints = dsgu_.getAccumulatedDataSize(subspacesWithData_);
    assert(numDataPoints > 0 && "all subspaces in dsg have 0 size");
    subspacesData_.resize(numDataPoints);
    std::memset(subspacesData_.data(), 0, subspacesData_.size() * sizeof(FG_ELEMENT));
    std::memset(subspaces_.data(), 0, subspaces_.size() * sizeof(FG_ELEMENT*));

    // update pointers and sizes in subspaces
    SubspaceSizeType offset = 0;
    for (size_t i = 0; i < subspaces_.size(); i++) {
      subspaces_[i] = subspacesData_.data() + offset;
      if (subspacesWithData_.find(i) != subspacesWithData_.end()) {
        offset += dsgu_.getSubspaceDataSizes()[i];
      }
    }
    assert(offset <= subspacesData_.size() && "offset exceeds data size");
    assert(std::is_sorted(std::begin(subspaces_), std::end(subspaces_)));
    if (kahanData_.empty()) {
      // create kahan buffer implicitly only once,
      // needs to be called explicitly if relevant sizes change
      this->createKahanBuffer();
    }
  }

  // allocates memory for kahan term data and sets pointers for it
  void createKahanBuffer() {
    if (subspacesWithData_.empty()) {
      // create data for all subspaces from "iota"
      subspacesWithData_ = std::set<SubspaceIndexType>{
          boost::counting_iterator<SubspaceIndexType>(0),
          boost::counting_iterator<SubspaceIndexType>(dsgu_.getNumSubspaces())};
    }
    size_t numDataPoints = dsgu_.getAccumulatedDataSize(subspacesWithData_);
    kahanData_.resize(numDataPoints, 0.);
    kahanDataBegin_.resize(dsgu_.getSubspaceDataSizes().size());
    std::memset(kahanDataBegin_.data(), 0, kahanDataBegin_.size() * sizeof(FG_ELEMENT*));

    // update pointers for begin of subspacen in kahan buffer
    SubspaceSizeType offset = 0;
    for (size_t i = 0; i < kahanDataBegin_.size(); i++) {
      kahanDataBegin_[i] = kahanData_.data() + offset;
      if (subspacesWithData_.find(i) != subspacesWithData_.end()) {
        offset += dsgu_.getSubspaceDataSizes()[i];
      }
    }
  }

  // deletes memory for subspace data and invalidates pointers to subspaces
  void deleteSubspaceData() {
    subspacesData_.clear();
    kahanData_.clear();
    // update pointers in subspaces
    std::memset(subspaces_.data(), 0, subspaces_.size() * sizeof(FG_ELEMENT*));
    std::memset(kahanDataBegin_.data(), 0, kahanDataBegin_.size() * sizeof(FG_ELEMENT*));
  }

  // sets all data elements to value zero
  void setZero() {
    std::memset(subspacesData_.data(), 0, subspacesData_.size() * sizeof(FG_ELEMENT));
    std::memset(kahanData_.data(), 0, kahanData_.size() * sizeof(FG_ELEMENT));
  }

  // returns the number of allocated grid points == size of the raw data vector
  inline size_t getRawDataSize() const { return subspacesData_.size(); }

  // returns true if data for the subspaces has been created
  bool isSubspaceDataCreated() const { return !subspacesData_.empty(); }

  void swap(DistributedSparseGridDataContainer& other) {
    assert(&dsgu_ == &other.dsgu_ && "cannot swap data containers of different dsgs");
    subspacesWithData_.swap(other.subspacesWithData_);
    subspaces_.swap(other.subspaces_);
    subspacesData_.swap(other.subspacesData_);
    kahanDataBegin_.swap(other.kahanDataBegin_);
    kahanData_.swap(other.kahanData_);
  }

  void allocateDifferentSubspaces(std::set<SubspaceIndexType>&& subspaces) {
    subspacesWithData_ = std::move(subspaces);
    subspacesData_.clear();
    kahanData_.clear();
    createSubspaceData();
    createKahanBuffer();
    setZero();
  }

 private:
  friend class DistributedSparseGridUniform<FG_ELEMENT>;

  const DistributedSparseGridUniform<FG_ELEMENT>&
      dsgu_;  // a reference to the dsgu to whose subspaces it belongs

  std::set<SubspaceIndexType> subspacesWithData_;  // set of subspaces currently allocated

  std::vector<FG_ELEMENT*> subspaces_;  // pointers to subspaces of the dsg

  std::vector<FG_ELEMENT> subspacesData_;  // allows linear access to all subspaces data

  std::vector<FG_ELEMENT*> kahanDataBegin_;  // pointers to Kahan summation residual terms

  std::vector<FG_ELEMENT> kahanData_;  // Kahan summation residual terms
};

/* This class can store a distributed sparse grid with a uniform space
 * decomposition. During construction no data is created and the data size of
 * the subspaces is initialized to zero (data sizes are usually set the
 * first time by registering the dsg in a distributed fullgrid during local
 * reduce). The data can be explicitly created by calling the
 * createSubspaceData() method or is implicitly generated as soon as it is
 * accessed. By calling deleteSubspaceData() the data can be deallocated.
 */
template <typename FG_ELEMENT>
class DistributedSparseGridUniform : public AnyDistributedSparseGrid {
 public:
  using ElementType = FG_ELEMENT;
  // type used to index the subspaces
  // should be enough for the current scenario (cf. test_createTruncatedHierarchicalLevels_large)
  // using SubspaceIndexType = AnyDistributedSparseGrid::SubspaceIndexType;

  /** create sparse grid of dimension d and specify for each dimension the
   * maximum discretization level
   * No data is allocated and the sizes of the subspace data is initialized to 0.
   */
  explicit DistributedSparseGridUniform(DimType dim, const LevelVector& lmax,
                                        const LevelVector& lmin, CommunicatorType comm);

  /**
   * create an empty (no data) sparse grid with given subspaces.
   */
  explicit DistributedSparseGridUniform(DimType dim, const std::vector<LevelVector>& subspaces,
                                        CommunicatorType comm);

  virtual ~DistributedSparseGridUniform() = default;

  // cheap rule of 5
  DistributedSparseGridUniform() = delete;
  DistributedSparseGridUniform(const DistributedSparseGridUniform& other) = delete;
  DistributedSparseGridUniform& operator=(const DistributedSparseGridUniform&) = delete;
  DistributedSparseGridUniform(DistributedSparseGridUniform&& other) = delete;
  DistributedSparseGridUniform& operator=(DistributedSparseGridUniform&& other) = delete;

  void swapDataContainers(DistributedSparseGridDataContainer<FG_ELEMENT>& otherContainer);

  void allocateDifferentSubspaces(std::set<SubspaceIndexType>&& subspaces);

  const std::set<SubspaceIndexType>& getCurrentlyAllocatedSubspaces() const;

  bool isSubspaceCurrentlyAllocated(SubspaceIndexType subspaceIndex) const;

  void print(std::ostream& os) const;

  // allocates memory for subspace data and sets pointers to subspaces
  void createSubspaceData();

  // allocates memory for kahan term data and sets pointers for it
  void createKahanBuffer();

  // deletes memory for subspace data and invalids pointers to subspaces
  void deleteSubspaceData();

  // creates data if necessary and sets all data elements to zero
  void setZero();

  // return all level vectors
  inline const std::vector<LevelVector>& getAllLevelVectors() const;

  // return level vector of subspace i
  inline const LevelVector& getLevelVector(SubspaceIndexType i) const;

  inline SubspaceIndexType getIndexInRange(const LevelVector& l, IndexType lowerBound) const;

  // return index of subspace i
  inline SubspaceIndexType getIndex(const LevelVector& l) const;

  // returns a pointer to first element in subspace i
  inline FG_ELEMENT* getData(SubspaceIndexType i);

  // returns a const pointer to first element in subspace i
  inline const FG_ELEMENT* getData(SubspaceIndexType i) const;

  // allows a linear access to the whole subspace data stored in this dsg
  inline FG_ELEMENT* getRawData();

  inline const FG_ELEMENT* getRawData() const;

  inline DimType getDim() const;

  // clear the levels_ vector, it is only necessary for registration of full grids
  void resetLevels();

  inline void setDataSize(SubspaceIndexType i, SubspaceSizeType newSize) override;

  SubspaceSizeType getAllocatedDataSize(SubspaceIndexType i) const;

  inline void registerDistributedFullGrid(const DistributedFullGrid<FG_ELEMENT>& dfg);

  template <bool sparseGridFullyAllocated = true>
  inline void addDistributedFullGrid(const DistributedFullGrid<FG_ELEMENT>& dfg,
                                     combigrid::real coeff);

  // returns the number of allocated grid points == size of the raw data vector
  inline size_t getRawDataSize() const;

  // returns true if data for the subspaces has been created
  bool isSubspaceDataCreated() const;

  // max-reduce subspace sizes from another DSGU (which has the same subspaces, but they may be less
  // or more populated than in this DSGU)
  void maxReduceSubspaceSizes(const DistributedSparseGridUniform<FG_ELEMENT>& other);

  // copy data from another DSGU (which has the same subspaces, but they may be less or more
  // populated than in this DSGU)
  void copyDataFrom(const DistributedSparseGridUniform<FG_ELEMENT>& other);

  template <typename SubspaceIndexContainer>
  void copyDataFrom(const DistributedSparseGridUniform<FG_ELEMENT>& other,
                    const SubspaceIndexContainer& subspaceIndices);

 private:
  std::vector<LevelVector> createLevels(DimType dim, const LevelVector& nmax,
                                        const LevelVector& lmin) const;
  DimType dim_;

  std::vector<LevelVector> levels_;  // linear access to all subspaces; may be reset to save memory

  DistributedSparseGridDataContainer<FG_ELEMENT> subspacesDataContainer_;
};

}  // namespace combigrid

namespace combigrid {

// at construction create only levels, no data
template <typename FG_ELEMENT>
DistributedSparseGridUniform<FG_ELEMENT>::DistributedSparseGridUniform(DimType dim,
                                                                       const LevelVector& lmax,
                                                                       const LevelVector& lmin,
                                                                       CommunicatorType comm)
    : DistributedSparseGridUniform(dim, createLevels(dim, lmax, lmin), comm) {}

// at construction create only levels, no data
template <typename FG_ELEMENT>
DistributedSparseGridUniform<FG_ELEMENT>::DistributedSparseGridUniform(
    DimType dim, const std::vector<LevelVector>& subspaces, CommunicatorType comm)
    : AnyDistributedSparseGrid(subspaces.size(), comm),
      dim_(dim),
      levels_(subspaces),
      subspacesDataContainer_(*this) {
  assert(dim > 0);

  for (auto& l : subspaces) {
    assert(l.size() == dim);
    for (size_t i = 0; i < l.size(); ++i) {
      assert(l[i] > 0);
    }
  }
}

template <typename FG_ELEMENT>
bool DistributedSparseGridUniform<FG_ELEMENT>::isSubspaceDataCreated() const {
  return subspacesDataContainer_.isSubspaceDataCreated();
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::maxReduceSubspaceSizes(
    const DistributedSparseGridUniform<FG_ELEMENT>& other) {
  assert(this->getNumSubspaces() == other.getNumSubspaces());
#pragma omp parallel for default(none) shared(other) schedule(static)
  for (decltype(this->getNumSubspaces()) i = 0; i < this->getNumSubspaces(); ++i) {
    assert(other.getDataSize(i) == this->getDataSize(i) || this->getDataSize(i) == 0 ||
           other.getDataSize(i) == 0);
    this->setDataSize(i, std::max(this->getDataSize(i), other.getDataSize(i)));
  }
  this->deleteSubspaceData();
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::copyDataFrom(
    const DistributedSparseGridUniform<FG_ELEMENT>& other) {
  assert(this->getNumSubspaces() == other.getNumSubspaces());
  assert(this->isSubspaceDataCreated());
  assert(other.isSubspaceDataCreated());

#pragma omp parallel for default(none) shared(other) schedule(guided)
  for (decltype(this->getNumSubspaces()) i = 0; i < this->getNumSubspaces(); ++i) {
    assert(other.getDataSize(i) == this->getDataSize(i) || this->getDataSize(i) == 0 ||
           other.getDataSize(i) == 0);
    assert(other.getAllocatedDataSize(i) == this->getAllocatedDataSize(i) ||
           this->getAllocatedDataSize(i) == 0 || other.getAllocatedDataSize(i) == 0);
    auto numPointsToCopy = std::min(other.getAllocatedDataSize(i), this->getAllocatedDataSize(i));
    std::copy_n(other.getData(i), numPointsToCopy, this->getData(i));
  }
}

template <typename FG_ELEMENT>
template <typename SubspaceIndexContainer>
void DistributedSparseGridUniform<FG_ELEMENT>::copyDataFrom(
    const DistributedSparseGridUniform<FG_ELEMENT>& other,
    const SubspaceIndexContainer& subspaceIndices) {
  assert(this->getNumSubspaces() == other.getNumSubspaces());
  assert(this->isSubspaceDataCreated() && other.isSubspaceDataCreated());

  // TODO this doesnt work with std::set
  // #pragma omp parallel for default(none) shared(other,subspaceIndices) schedule(guided)
  // for (auto it = subspaceIndices.cbegin(); it != subspaceIndices.cend(); ++it) {
  for (auto i : subspaceIndices) {
    assert(other.getDataSize(i) == this->getDataSize(i));
    assert(other.getAllocatedDataSize(i) == this->getAllocatedDataSize(i));
    auto numPointsToCopy = std::min(other.getDataSize(i), this->getDataSize(i));
    std::copy_n(other.getData(i), numPointsToCopy, this->getData(i));
  }
}

/** Zero initializes the dsgu data in case no data is already present.
 * Otherwise nothing happens.
 */
template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::createSubspaceData() {
  subspacesDataContainer_.createSubspaceData();
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::createKahanBuffer() {
  subspacesDataContainer_.createKahanBuffer();
}

/** Deallocates the dsgu data.
 *  This affects the values stored at the grid points and pointers which address
 *  the subspaces data.
 */
template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::deleteSubspaceData() {
  subspacesDataContainer_.deleteSubspaceData();
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::setZero() {
  subspacesDataContainer_.setZero();
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::print(std::ostream& os) const {
  auto levelIterator = levels_.cbegin();
  for (size_t i = 0; i < this->subspacesDataContainer_.subspaces_.size(); ++i) {
    os << i << " " << *levelIterator << " " << this->subspacesDataSizes_[i]
       << " "
       //  << subspaces_[i].size_
       << std::endl;
    ++levelIterator;
  }
}

template <typename FG_ELEMENT>
std::ostream& operator<<(std::ostream& os, const DistributedSparseGridUniform<FG_ELEMENT>& sg) {
  sg.print(os);
  return os;
}

template <typename FG_ELEMENT>
std::vector<LevelVector> DistributedSparseGridUniform<FG_ELEMENT>::createLevels(
    DimType dim, const LevelVector& nmax, const LevelVector& lmin) const {
  std::vector<LevelVector> created{};
  combigrid::createTruncatedHierarchicalLevels(nmax, lmin, created);
  // std::sort(created.begin(), created.end());
  assert(std::is_sorted(created.begin(), created.end()));
  if (created.size() > std::numeric_limits<SubspaceIndexType>::max()) {
    throw std::runtime_error("number of subspaces exceeds the maximum value of SubspaceIndexType");
  }
  return created;
}

template <typename FG_ELEMENT>
inline const std::vector<LevelVector>&
DistributedSparseGridUniform<FG_ELEMENT>::getAllLevelVectors() const {
  return levels_;
}

template <typename FG_ELEMENT>
inline const LevelVector& DistributedSparseGridUniform<FG_ELEMENT>::getLevelVector(
    SubspaceIndexType i) const {
  auto levelIterator = levels_.cbegin();
  std::advance(levelIterator, i);
  return *levelIterator;
}

template <typename FG_ELEMENT>
typename DistributedSparseGridUniform<FG_ELEMENT>::SubspaceIndexType
DistributedSparseGridUniform<FG_ELEMENT>::getIndexInRange(const LevelVector& l,
                                                          IndexType lowerBound) const {
#ifndef NDEBUG
  for (const auto& l_i : l) {
    assert(l_i > 0);
  }
#endif  // NDEBUG
  auto start = levels_.cbegin();
  std::advance(start, lowerBound);
  auto found = std::lower_bound(start, levels_.end(), l);
  if (found != levels_.end() && *found == l) {
    return static_cast<SubspaceIndexType>(std::distance(levels_.cbegin(), found));
  } else {
    // assert(false && "space not found in levels_");
    return -1;
  }
}

/* get index of space with l. returns -1 if not included */
template <typename FG_ELEMENT>
typename DistributedSparseGridUniform<FG_ELEMENT>::SubspaceIndexType
DistributedSparseGridUniform<FG_ELEMENT>::getIndex(const LevelVector& l) const {
  return getIndexInRange(l, 0);
}

template <typename FG_ELEMENT>
inline FG_ELEMENT* DistributedSparseGridUniform<FG_ELEMENT>::getData(SubspaceIndexType i) {
#ifndef NDEBUG
  assert(isSubspaceDataCreated());
  assert(i < this->subspacesDataContainer_.subspaces_.size());
  assert(this->subspacesDataContainer_.subspaces_[i] <=
         &(*(this->subspacesDataContainer_.subspacesData_.end())));
  // if (this->subspacesDataContainer_.subspaces_[i] ==
  //     &(*(this->subspacesDataContainer_.subspacesData_.end()))) {
  //   assert(this->getDataSize(i) == 0);
  // }
#endif
  return this->subspacesDataContainer_.subspaces_[i];
}

template <typename FG_ELEMENT>
inline const FG_ELEMENT* DistributedSparseGridUniform<FG_ELEMENT>::getData(
    SubspaceIndexType i) const {
#ifndef NDEBUG
  assert(isSubspaceDataCreated());
  assert(i < this->subspacesDataContainer_.subspaces_.size());
  assert(this->subspacesDataContainer_.subspaces_[i] <=
         &(*(this->subspacesDataContainer_.subspacesData_.end())));
  if (this->subspacesDataContainer_.subspaces_[i] ==
      &(*(this->subspacesDataContainer_.subspacesData_.end()))) {
    assert(this->getDataSize(i) == 0);
  }
#endif
  return this->subspacesDataContainer_.subspaces_[i];
}

template <typename FG_ELEMENT>
inline FG_ELEMENT* DistributedSparseGridUniform<FG_ELEMENT>::getRawData() {
  return this->subspacesDataContainer_.subspacesData_.data();
}

template <typename FG_ELEMENT>
inline const FG_ELEMENT* DistributedSparseGridUniform<FG_ELEMENT>::getRawData() const {
  assert(isSubspaceDataCreated() && "subspace data not created");
  return this->subspacesDataContainer_.subspacesData_.data();
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::swapDataContainers(
    DistributedSparseGridDataContainer<FG_ELEMENT>& otherContainer) {
  subspacesDataContainer_.swap(otherContainer);
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::allocateDifferentSubspaces(
    std::set<typename AnyDistributedSparseGrid::SubspaceIndexType>&& subspaces) {
  subspacesDataContainer_.allocateDifferentSubspaces(std::move(subspaces));
}

template <typename FG_ELEMENT>
const std::set<typename AnyDistributedSparseGrid::SubspaceIndexType>&
DistributedSparseGridUniform<FG_ELEMENT>::getCurrentlyAllocatedSubspaces() const {
  return subspacesDataContainer_.subspacesWithData_;
}

template <typename FG_ELEMENT>
bool DistributedSparseGridUniform<FG_ELEMENT>::isSubspaceCurrentlyAllocated(
    SubspaceIndexType subspaceIndex) const {
  // try to find value in subspacesWithData_
  // if it is not found, the subspace is not allocated
  return subspacesDataContainer_.subspacesWithData_.find(subspaceIndex) !=
         subspacesDataContainer_.subspacesWithData_.end();
}

template <typename FG_ELEMENT>
inline DimType DistributedSparseGridUniform<FG_ELEMENT>::getDim() const {
  return dim_;
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::resetLevels() {
  levels_.clear();
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::setDataSize(SubspaceIndexType i,
                                                           SubspaceSizeType newSize) {
#ifndef NDEBUG
  assert(subspacesDataSizes_[i] == 0 || subspacesDataSizes_[i] == newSize);
  if (i >= getNumSubspaces()) {
    std::cout << "Index too large, no subspace with this index included in distributed sparse grid"
              << std::endl;
    assert(false);
  }
#endif  // NDEBUG
  if (newSize != subspacesDataSizes_[i]) {
    // invalidate the data vector
    this->deleteSubspaceData();
  }
  subspacesDataSizes_[i] = newSize;
}

template <typename FG_ELEMENT>
SubspaceSizeType DistributedSparseGridUniform<FG_ELEMENT>::getAllocatedDataSize(
    SubspaceIndexType i) const {
  if (this->isSubspaceCurrentlyAllocated(i)) {
    return subspacesDataSizes_[i];
  } else {
    return 0;
  }
}

/**
 * @brief "registers" the DistributedFullGrid with this DistributedSparseGridUniform:
 *         sets the dsg's subspaceSizes where they are not yet set to contain all of the DFG's
 *         subspaces.
 *        The size of a subspace in the dsg is chosen according to the corresponding
 *        subspace size in the dfg.
 *
 * @param dfg the DFG to register
 */
template <typename FG_ELEMENT>
inline void DistributedSparseGridUniform<FG_ELEMENT>::registerDistributedFullGrid(
    const DistributedFullGrid<FG_ELEMENT>& dfg) {
  assert(dfg.getDimension() == dim_);
  // all the hierarchical subspaces contained in the full grid
  const auto downwardClosedSet = combigrid::getDownSet(dfg.getLevels());

  SubspaceIndexType index = 0;
  // resize all common subspaces in dsg, if necessary
#pragma omp parallel for default(none) shared(downwardClosedSet, dfg, std::cout, std::cerr) \
    firstprivate(index) schedule(guided)
  for (const auto& level : downwardClosedSet) {
    index = this->getIndexInRange(level, index);
    if (index > -1) {
      IndexType numPointsOfSubspace = 1;
      for (DimType d = 0; d < dim_; ++d) {
        numPointsOfSubspace *= dfg.getNumPointsOnThisPartition(level[d], d);
      }
      const auto subSgDataSize = this->getDataSize(index);
      // resize DSG subspace if it has zero size
      if (subSgDataSize == 0) {
        this->setDataSize(index, numPointsOfSubspace);
      } else {
        ASSERT(subSgDataSize == numPointsOfSubspace,
               "subSgDataSize: " << subSgDataSize
                                 << ", numPointsOfSubspace: " << numPointsOfSubspace << " , level "
                                 << level << " , rank " << this->rank_ << std::endl);
      }
    }
  }
}

/**
 * @brief adds the (hopefully) hierarchical coefficients from the DFG
 *        to the DSG's data structure, multiplied by coeff
 *
 * @param dfg the DFG to add
 * @param coeff the coefficient that gets multiplied to all entries in DFG
 */
template <typename FG_ELEMENT>
template <bool sparseGridFullyAllocated>
inline void DistributedSparseGridUniform<FG_ELEMENT>::addDistributedFullGrid(
    const DistributedFullGrid<FG_ELEMENT>& dfg, combigrid::real coeff) {
  assert(this->isSubspaceDataCreated());
  if (this->subspacesDataContainer_.kahanData_.empty() ||
      this->subspacesDataContainer_.kahanDataBegin_.empty()) {
    throw std::runtime_error("Kahan data not initialized");
  }

  // all the hierarchical subspaces contained in this full grid
  const auto downwardClosedSet = combigrid::getDownSet(dfg.getLevels());

  static thread_local IndexVector subspaceIndices;
  SubspaceIndexType sIndex = 0;
// loop over all subspaces of the full grid
#pragma omp parallel for default(none) shared(downwardClosedSet, dfg, std::cout) \
    firstprivate(coeff, sIndex) schedule(guided)
  for (const auto& level : downwardClosedSet) {
    sIndex = this->getIndexInRange(level, sIndex);
    bool shouldBeCopied = sIndex > -1 && this->getDataSize(sIndex) > 0;
    if constexpr (!sparseGridFullyAllocated) {
      shouldBeCopied = shouldBeCopied && this->isSubspaceCurrentlyAllocated(sIndex);
    }
    if (shouldBeCopied) {
      auto sPointer = this->getData(sIndex);
      auto kPointer = this->subspacesDataContainer_.kahanDataBegin_[sIndex];
#ifndef NDEBUG
      if (sIndex < this->getNumSubspaces() - 1) {
        auto sDataSize = this->subspacesDataContainer_.subspaces_[sIndex + 1] - sPointer;
        assert(sDataSize == this->getDataSize(sIndex));
        assert(sDataSize == this->getAllocatedDataSize(sIndex));
        assert(std::distance(this->getData(0), sPointer) <
               this->subspacesDataContainer_.subspacesData_.size());
        auto kDataSize = this->subspacesDataContainer_.kahanDataBegin_[sIndex + 1] - kPointer;
        assert(kDataSize == this->getDataSize(sIndex));
      }
#endif  // NDEBUG
      subspaceIndices = std::move(dfg.getFGPointsOfSubspace(level));
#pragma omp simd linear(sPointer, kPointer : 1)
      for (size_t fIndex = 0; fIndex < subspaceIndices.size(); ++fIndex) {
        FG_ELEMENT summand = coeff * dfg.getData()[subspaceIndices[fIndex]];
        // cf. https://en.wikipedia.org/wiki/Kahan_summation_algorithm
        FG_ELEMENT y = summand - *kPointer;  // TODO check if these should be volatile
        FG_ELEMENT t = *sPointer + y;
        *kPointer = (t - *sPointer) - y;
        *sPointer = t;
        ++sPointer;
        ++kPointer;
      }
      assert(this->getDataSize(sIndex) == subspaceIndices.size());
    }
  }
}

template <typename FG_ELEMENT>
inline size_t DistributedSparseGridUniform<FG_ELEMENT>::getRawDataSize() const {
  return subspacesDataContainer_.getRawDataSize();
}

} /* namespace combigrid */

#endif /* SRC_SGPP_COMBIGRID_SPARSEGRID_DISTRIBUTEDSPARSEGRID_HPP_ */
