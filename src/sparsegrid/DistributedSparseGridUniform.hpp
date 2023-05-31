#ifndef SRC_SGPP_COMBIGRID_SPARSEGRID_DISTRIBUTEDSPARSEGRIDUNIFORM_HPP_
#define SRC_SGPP_COMBIGRID_SPARSEGRID_DISTRIBUTEDSPARSEGRIDUNIFORM_HPP_

#include <cassert>

#include "utils/Types.hpp"
#include "utils/LevelSetUtils.hpp"
#include "sparsegrid/AnyDistributedSparseGrid.hpp"
#include <numeric>

namespace combigrid {
// forward declaration
template <typename FG_ELEMENT>
class DistributedFullGrid;

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
  using SubspaceIndexType = int32_t;

  /** create sparse grid of dimension d and specify for each dimension the
   * maximum discretization level
   * No data is allocated and the sizes of the subspace data is initialized to 0.
   */
  explicit DistributedSparseGridUniform(DimType dim, const LevelVector& lmax, const LevelVector& lmin,
                               CommunicatorType comm);

  /**
   * create an empty (no data) sparse grid with given subspaces.
   */
  explicit DistributedSparseGridUniform(DimType dim, const std::vector<LevelVector>& subspaces,
                               CommunicatorType comm);

  virtual ~DistributedSparseGridUniform();

  // cheap rule of 5
  DistributedSparseGridUniform() = delete;
  DistributedSparseGridUniform(const DistributedSparseGridUniform& other) = delete;
  DistributedSparseGridUniform& operator=(const DistributedSparseGridUniform&) = delete;
  DistributedSparseGridUniform(DistributedSparseGridUniform&& other) = delete;
  DistributedSparseGridUniform& operator=(DistributedSparseGridUniform&& other) = delete;

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

  inline void registerDistributedFullGrid(const DistributedFullGrid<FG_ELEMENT>& dfg);

  inline void addDistributedFullGrid(const DistributedFullGrid<FG_ELEMENT>& dfg,
                                     combigrid::real coeff);

  // returns the number of allocated grid points == size of the raw data vector
  inline size_t getRawDataSize() const;

  // returns true if data for the subspaces has been created
  bool isSubspaceDataCreated() const;

  // copy data from another DSGU (which has the same subspaces, but they may be less or more
  // populated than in this DSGU)
  void copyDataFrom(const DistributedSparseGridUniform<FG_ELEMENT>& other);

  const std::vector<std::pair<CommunicatorType, MPI_Datatype>>& getDatatypesByComm() const;

 private:
  std::vector<LevelVector> createLevels(DimType dim, const LevelVector& nmax,
                                        const LevelVector& lmin) const;
  void setReductionDatatypes();

  DimType dim_;

  std::vector<LevelVector> levels_;  // linear access to all subspaces; may be reset to save memory

  std::vector<FG_ELEMENT*> subspaces_;  // pointers to subspaces of the dsg

  std::vector<FG_ELEMENT> subspacesData_;  // allows linear access to all subspaces data

  std::vector<FG_ELEMENT*> kahanDataBegin_;  // pointers to Kahan summation residual terms

  std::vector<FG_ELEMENT> kahanData_;  // Kahan summation residual terms

  std::vector<std::pair<CommunicatorType, MPI_Datatype>> datatypesByComm_;
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
      subspaces_(subspaces.size()) {
  assert(dim > 0);

  for (auto& l : subspaces) {
    assert(l.size() == dim);
    for (size_t i = 0; i < l.size(); ++i) {
      assert(l[i] > 0);
    }
  }
}

template <typename FG_ELEMENT>
DistributedSparseGridUniform<FG_ELEMENT>::~DistributedSparseGridUniform() {
  // free datatypes and comms
  for (auto& it : datatypesByComm_) {
    MPI_Type_free(&it.second);
    // MPI_Comm_free(&it.first);
  }
}

template <typename FG_ELEMENT>
bool DistributedSparseGridUniform<FG_ELEMENT>::isSubspaceDataCreated() const {
  return subspacesData_.size() != 0;
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::copyDataFrom(
    const DistributedSparseGridUniform<FG_ELEMENT>& other) {
  assert(this->isSubspaceDataCreated() && other.isSubspaceDataCreated());
  // #pragma omp parallel for
  for (decltype(this->getNumSubspaces()) i = 0; i < this->getNumSubspaces(); ++i) {
    assert(other.getDataSize(i) == this->getDataSize(i) || this->getDataSize(i) == 0 ||
           other.getDataSize(i) == 0);
    auto numPointsToCopy = std::min(other.getDataSize(i), this->getDataSize(i));
    std::copy_n(other.getData(i), numPointsToCopy, this->getData(i));
  }
}

/** Zero initializes the dsgu data in case no data is already present.
 * Otherwise nothing happens.
 */
template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::createSubspaceData() {
  if (not isSubspaceDataCreated()) {
    size_t numDataPoints = this->getAccumulatedDataSize();
    assert(numDataPoints > 0 && "all subspaces in dsg have 0 size");
    subspacesData_.resize(numDataPoints, 0.);

    // update pointers and sizes in subspaces
    SubspaceSizeType offset = 0;
    for (size_t i = 0; i < subspaces_.size(); i++) {
      subspaces_[i] = subspacesData_.data() + offset;
      offset += subspacesDataSizes_[i];
    }
    if (kahanData_.empty()) {
      // create kahan buffer implicitly only once,
      // needs to be called explicitly if relevant sizes change
      this->createKahanBuffer();
    }

    // if the subspacesByComm_ is set, we need to set the MPI data types for reduction
    if (!this->subspacesByComm_.empty()) {
      this->setReductionDatatypes();
    }
  }
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::createKahanBuffer() {
  size_t numDataPoints = this->getAccumulatedDataSize();
  kahanData_.resize(numDataPoints, 0.);
  kahanDataBegin_.resize(this->subspacesDataSizes_.size());

  // update pointers for begin of subspacen in kahan buffer
  SubspaceSizeType offset = 0;
  for (size_t i = 0; i < kahanDataBegin_.size(); i++) {
    kahanDataBegin_[i] = kahanData_.data() + offset;
    offset += this->subspacesDataSizes_[i];
  }
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::setReductionDatatypes() {
  // like for sparse grid reduce, allow only up to 16MiB per reduction
  //(when using double precision)
  auto chunkSize = 2097152;

  // iterate subspacesByComm_ and create MPI datatypes for each communicator
  FG_ELEMENT* rawDataStart = this->getRawData();
  for (auto& it : this->subspacesByComm_) {
    auto comm = it.first;
    const auto& subspaces = it.second;
    // get chunked subspaces for this data type
    {
      auto subspaceIt = subspaces.cbegin();
      while (subspaceIt != subspaces.cend()) {
        SubspaceSizeType chunkDataSize = 0;
        std::vector<SubspaceIndexType> chunkSubspaces;
        auto nextAddedDataSize = getDataSize(*subspaceIt);
        do {
          chunkDataSize += nextAddedDataSize;
          chunkSubspaces.push_back(*subspaceIt);
          ++subspaceIt;
        } while (subspaceIt != subspaces.cend() && (nextAddedDataSize = getDataSize(*subspaceIt)) &&
                 (chunkDataSize + nextAddedDataSize) < chunkSize);

        // create datatype for this chunk
        std::vector<int> arrayOfBlocklengths, arrayOfDisplacements;
        arrayOfBlocklengths.reserve(chunkSubspaces.size());
        arrayOfDisplacements.reserve(chunkSubspaces.size());
        for (const auto& ss : chunkSubspaces) {
          arrayOfBlocklengths.push_back(this->getDataSize(ss));
          arrayOfDisplacements.push_back(this->getData(ss) - rawDataStart);
        }

        MPI_Datatype myIndexedDatatype;
        MPI_Type_indexed(static_cast<int>(chunkSubspaces.size()), arrayOfBlocklengths.data(),
                         arrayOfDisplacements.data(),
                         getMPIDatatype(abstraction::getabstractionDataType<FG_ELEMENT>()),
                         &myIndexedDatatype);
        MPI_Type_commit(&myIndexedDatatype);
        datatypesByComm_.push_back(std::make_pair(comm, myIndexedDatatype));
      }
    }
  }
}

/** Deallocates the dsgu data.
 *  This affects the values stored at the grid points and pointers which address
 *  the subspaces data.
 */
template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::deleteSubspaceData() {
  if (isSubspaceDataCreated()) {
    subspacesData_.clear();

    // update pointers in subspaces
    for (auto& ss : subspaces_) {
      ss = nullptr;
    }
  }
  // free datatypes and ops
  for (auto& it : datatypesByComm_) {
    MPI_Type_free(&it.second);
    // MPI_Comm_free(&it.first);
  }
  datatypesByComm_.clear();
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::setZero() {
  if (isSubspaceDataCreated())
    std::memset(subspacesData_.data(), 0, subspacesData_.size() * sizeof(FG_ELEMENT));
  else
    createSubspaceData();
  std::memset(kahanData_.data(), 0, kahanData_.size() * sizeof(FG_ELEMENT));
  if (kahanData_.empty()) {
    // create kahan buffer implicitly only once,
    // needs to be called explicitly if relevant sizes change
    this->createKahanBuffer();
  }
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::print(std::ostream& os) const {
  auto levelIterator = levels_.cbegin();
  for (size_t i = 0; i < subspaces_.size(); ++i) {
    os << i << " " << *levelIterator << " " << subspacesDataSizes_[i] << " "
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
  assert(isSubspaceDataCreated());
  return subspaces_[i];
}

template <typename FG_ELEMENT>
inline const FG_ELEMENT* DistributedSparseGridUniform<FG_ELEMENT>::getData(
    SubspaceIndexType i) const {
  assert(isSubspaceDataCreated());
  return subspaces_[i];
}

template <typename FG_ELEMENT>
const std::vector<std::pair<CommunicatorType, MPI_Datatype>>&
DistributedSparseGridUniform<FG_ELEMENT>::getDatatypesByComm() const {
  return datatypesByComm_;
}

template <typename FG_ELEMENT>
inline FG_ELEMENT* DistributedSparseGridUniform<FG_ELEMENT>::getRawData() {
  return subspacesData_.data();
}

template <typename FG_ELEMENT>
inline const FG_ELEMENT* DistributedSparseGridUniform<FG_ELEMENT>::getRawData() const {
  assert(isSubspaceDataCreated() && "subspace data not created");
  return subspacesData_.data();
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
  IndexType numPointsOfSubspace = 1;
  // resize all common subspaces in dsg, if necessary
  for (const auto& level : downwardClosedSet) {
    index = this->getIndexInRange(level, index);
    if (index > -1) {
      numPointsOfSubspace = 1;
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
inline void DistributedSparseGridUniform<FG_ELEMENT>::addDistributedFullGrid(
    const DistributedFullGrid<FG_ELEMENT>& dfg, combigrid::real coeff) {
  assert(this->isSubspaceDataCreated());
  if (kahanData_.empty() || kahanDataBegin_.empty()) {
    throw std::runtime_error("Kahan data not initialized");
  }

  bool anythingWasAdded = false;

  // all the hierarchical subspaces contained in this full grid
  const auto downwardClosedSet = combigrid::getDownSet(dfg.getLevels());
  static IndexVector subspaceIndices;

  SubspaceIndexType sIndex = 0;
  // loop over all subspaces of the full grid
  for (const auto& level : downwardClosedSet) {
    sIndex = this->getIndexInRange(level, sIndex);
    if (sIndex > -1 && this->getDataSize(sIndex) > 0) {
      auto sPointer = this->getData(sIndex);
      auto kPointer = kahanDataBegin_[sIndex];
#ifndef NDEBUG
      if (sIndex < kahanDataBegin_.size() - 1) {
        auto kDataSize = kahanDataBegin_[sIndex + 1] - kPointer;
        assert(kDataSize == this->getDataSize(sIndex));
      }
#endif  // NDEBUG
      subspaceIndices = dfg.getFGPointsOfSubspace(level);
      for (const auto& fIndex : subspaceIndices) {
        FG_ELEMENT summand = coeff * dfg.getData()[fIndex];
        // cf. https://en.wikipedia.org/wiki/Kahan_summation_algorithm
        FG_ELEMENT y = summand - *kPointer;
        FG_ELEMENT t = *sPointer + y;
        *kPointer = (t - *sPointer) - y;
        *sPointer = t;
        ++sPointer;
        ++kPointer;
        anythingWasAdded = true;
      }
    }
  }

  // make sure that anything was added -- I can only think of weird setups
  // where that would not be the case
  assert(anythingWasAdded);
}

template <typename FG_ELEMENT>
inline size_t DistributedSparseGridUniform<FG_ELEMENT>::getRawDataSize() const {
  return subspacesData_.size();
}

} /* namespace combigrid */

#endif /* SRC_SGPP_COMBIGRID_SPARSEGRID_DISTRIBUTEDSPARSEGRID_HPP_ */
