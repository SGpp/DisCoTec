#ifndef SRC_SGPP_COMBIGRID_SPARSEGRID_DISTRIBUTEDSPARSEGRIDUNIFORM_HPP_
#define SRC_SGPP_COMBIGRID_SPARSEGRID_DISTRIBUTEDSPARSEGRIDUNIFORM_HPP_

#include <assert.h>

#include "utils/Types.hpp"
#include "utils/LevelSetUtils.hpp"
#include "manager/ProcessGroupSignals.hpp"
#include "mpi/MPITags.hpp"
#include "io/MPIInputOutput.hpp"
#include <numeric>

#include <boost/serialization/vector.hpp>

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
class DistributedSparseGridUniform {
 public:
  // type used to index the subspaces
  // should be enough for the current scenario (cf. test_createTruncatedHierarchicalLevels_large)
  using SubspaceIndexType = int32_t;

  /** create sparse grid of dimension d and specify for each dimension the
   * maximum discretization level
   * No data is allocated and the sizes of the subspace data is initialized to 0.
   */
  DistributedSparseGridUniform(DimType dim, const LevelVector& lmax, const LevelVector& lmin,
                               CommunicatorType comm);

  /**
   * create an empty (no data) sparse grid with given subspaces.
   */
  DistributedSparseGridUniform(DimType dim, const std::vector<LevelVector>& subspaces,
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

  // return the number of subspaces
  inline SubspaceIndexType getNumSubspaces() const;

  // check if a subspace with l is contained in the sparse grid
  // unlike getIndex this will not throw an assert in case l is not contained
  bool isContained(const LevelVector& l) const;

  // clear the levels_ vector, it is only necessary for registration of full grids
  void resetLevels();

  // data size of the subspace at index i
  inline SubspaceSizeType getDataSize(SubspaceIndexType i) const;

  // sum of all data sizes of all subspaces
  inline size_t getAccumulatedDataSize() const;

  // sets data size of subspace with index i to newSize
  inline void setDataSize(SubspaceIndexType i, SubspaceSizeType newSize);

  inline void registerDistributedFullGrid(const DistributedFullGrid<FG_ELEMENT>& dfg);

  inline void addDistributedFullGrid(const DistributedFullGrid<FG_ELEMENT>& dfg,
                                     combigrid::real coeff);

  // returns the number of allocated grid points == size of the raw data vector
  inline size_t getRawDataSize() const;

  inline CommunicatorType getCommunicator() const;

  inline int getCommunicatorSize() const;

  // allows linear access to the data sizes of all subspaces
  const std::vector<SubspaceSizeType>& getSubspaceDataSizes() const;

  // reduces the data sizes (between process groups) in-place
  void reduceSubspaceSizes(CommunicatorType comm);

  // broadcasts subspace sizes from one rank to all others in comm
  void broadcastDsgSizes(CommunicatorType comm, RankType sendingRank);

  void sendDsgSizesWithGather(CommunicatorType comm, RankType collectorRank);

  void receiveDsgSizesWithScatter(CommunicatorType comm, RankType collectorRank);

  // returns true if data for the subspaces has been created
  bool isSubspaceDataCreated() const;

  // copy data from another DSGU (which has the same subspaces, but they may be less or more populated than in this DSGU)
  void copyDataFrom(const DistributedSparseGridUniform<FG_ELEMENT>& other);

  void writeMinMaxCoefficents(const std::string& filename, size_t i) const;

  // naive read/write operations -- each rank writes their own data partition to a separate binary file
  void writeToDiskChunked(std::string filePrefix);

  void readFromDiskChunked(std::string filePrefix);

  // coordinated read/write to one single file containing the whole dsg data
  bool writeOneFile(std::string fileName) const;

  bool readOneFile(std::string fileName);

  bool readOneFileAndReduce(std::string fileName, int numberOfChunks = 1);

  bool writeSubspaceSizesToFile(std::string fileName) const;

  bool readSubspaceSizesFromFile(std::string fileName, bool withCollectiveBuffering = false);

  template <typename ReduceFunctionType>
  bool readReduceSubspaceSizesFromFile(std::string fileName, ReduceFunctionType reduceFunction,
                                       int numElementsToBuffer = 0,
                                       bool withCollectiveBuffering = false);

 private:
  std::vector<LevelVector> createLevels(DimType dim, const LevelVector& nmax,
                                        const LevelVector& lmin);

  DimType dim_;

  std::vector<LevelVector> levels_;  // linear access to all subspaces; may be reset to save memory

  CommunicatorType comm_;

  RankType rank_;

  int commSize_;

  std::vector<FG_ELEMENT*> subspaces_;  // pointers to subspaces of the dsg

  std::vector<FG_ELEMENT> subspacesData_;  // allows linear access to all subspaces data

  std::vector<SubspaceSizeType> subspacesDataSizes_;  // allocated data sizes of all subspaces

  std::vector<FG_ELEMENT*> kahanDataBegin_;  // pointers to Kahan summation residual terms

  std::vector<FG_ELEMENT> kahanData_;  // Kahan summation residual terms

  friend class boost::serialization::access;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version);
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
    : dim_(dim),
      levels_(subspaces),
      comm_(comm),
      subspaces_(subspaces.size()),
      subspacesDataSizes_(subspaces.size()) {
  assert(dim > 0);

  for (auto& l : subspaces) {
    assert(l.size() == dim);
    for (size_t i = 0; i < l.size(); ++i) {
      assert(l[i] > 0);
    }
  }
  MPI_Comm_rank(comm_, &rank_);
  MPI_Comm_size(comm_, &commSize_);
}

template <typename FG_ELEMENT>
bool DistributedSparseGridUniform<FG_ELEMENT>::isSubspaceDataCreated() const {
  return subspacesData_.size() != 0;
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::copyDataFrom(
    const DistributedSparseGridUniform<FG_ELEMENT>& other) {
  assert(this->isSubspaceDataCreated() && other.isSubspaceDataCreated());
  for (decltype(this->getNumSubspaces()) i = 0; i < this->getNumSubspaces(); ++i) {
    assert(other.getDataSize(i) == this->getDataSize(i) || this->getDataSize(i) == 0 || other.getDataSize(i) == 0);
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
  }
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::createKahanBuffer() {
  size_t numDataPoints = std::accumulate(subspacesDataSizes_.begin(), subspacesDataSizes_.end(),
                                         static_cast<size_t>(0));
  kahanData_.resize(numDataPoints, 0.);
  kahanDataBegin_.resize(subspacesDataSizes_.size());

  // update pointers for begin of subspacen in kahan buffer
  SubspaceSizeType offset = 0;
  for (size_t i = 0; i < kahanDataBegin_.size(); i++) {
    kahanDataBegin_[i] = kahanData_.data() + offset;
    offset += subspacesDataSizes_[i];
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
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::setZero() {
  if (isSubspaceDataCreated())
    std::fill(subspacesData_.begin(), subspacesData_.end(), 0.);
  else
    createSubspaceData();
  std::fill(kahanData_.begin(), kahanData_.end(), 0.);
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
DistributedSparseGridUniform<FG_ELEMENT>::~DistributedSparseGridUniform() {}

template <typename FG_ELEMENT>
std::vector<LevelVector> DistributedSparseGridUniform<FG_ELEMENT>::createLevels(
    DimType dim, const LevelVector& nmax, const LevelVector& lmin) {
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
inline typename DistributedSparseGridUniform<FG_ELEMENT>::SubspaceIndexType
DistributedSparseGridUniform<FG_ELEMENT>::getNumSubspaces() const {
  return static_cast<SubspaceIndexType>(subspaces_.size());
}

template <typename FG_ELEMENT>
bool DistributedSparseGridUniform<FG_ELEMENT>::isContained(const LevelVector& l) const {
  auto found = std::lower_bound(levels_.cbegin(), levels_.cend(), l);
  return (found != levels_.end() && *found == l);
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::resetLevels() {
  levels_.clear();
}

template <typename FG_ELEMENT>
SubspaceSizeType DistributedSparseGridUniform<FG_ELEMENT>::getDataSize(SubspaceIndexType i) const {
#ifndef NDEBUG
  if (i >= getNumSubspaces()) {
    std::cout << "Index too large, no subspace with this index included in distributed sparse grid"
              << std::endl;
    assert(false);
  }
#endif  // NDEBUG

  return subspacesDataSizes_[i];
}

template <typename FG_ELEMENT>
size_t DistributedSparseGridUniform<FG_ELEMENT>::getAccumulatedDataSize() const {
  return std::accumulate(subspacesDataSizes_.begin(), subspacesDataSizes_.end(),
                         static_cast<size_t>(0));
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
        FG_ELEMENT summand = coeff * dfg.getElementVector()[fIndex];
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

template <typename FG_ELEMENT>
CommunicatorType DistributedSparseGridUniform<FG_ELEMENT>::getCommunicator() const {
  return comm_;
}

template <typename FG_ELEMENT>
inline int DistributedSparseGridUniform<FG_ELEMENT>::getCommunicatorSize() const {
  return commSize_;
}

template <typename FG_ELEMENT>
inline const std::vector<SubspaceSizeType>&
DistributedSparseGridUniform<FG_ELEMENT>::getSubspaceDataSizes() const {
  return subspacesDataSizes_;
}

/** Performs a max allreduce in comm with subspace sizes of each dsg
 *
 * After calling, all workers which share the same spatial decomposition will
 * have the same subspace sizes and therefor. in the end have equally sized dsgs.
 */
template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::reduceSubspaceSizes(CommunicatorType comm) {
  assert(this->getNumSubspaces() > 0);

  // prepare for MPI call in globalReduceComm
  MPI_Datatype dtype = getMPIDatatype(abstraction::getabstractionDataType<SubspaceSizeType>());

  // perform allreduce
  assert(subspacesDataSizes_.size() <
         static_cast<SubspaceSizeType>(std::numeric_limits<int>::max()));
  MPI_Allreduce(MPI_IN_PLACE, subspacesDataSizes_.data(),
                static_cast<int>(subspacesDataSizes_.size()), dtype, MPI_MAX, comm);
  // assume that the sizes changed, the buffer might be the wrong size now
  this->deleteSubspaceData();
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::broadcastDsgSizes(CommunicatorType comm,
                                                                 RankType sendingRank) {
  assert(this->getNumSubspaces() > 0);
  MPI_Datatype dtype = getMPIDatatype(abstraction::getabstractionDataType<SubspaceSizeType>());

  // perform broadcast
  assert(subspacesDataSizes_.size() <
         static_cast<SubspaceSizeType>(std::numeric_limits<int>::max()));
  MPI_Bcast(subspacesDataSizes_.data(), static_cast<int>(subspacesDataSizes_.size()), dtype,
            sendingRank, comm);
  // assume that the sizes changed, the buffer might be the wrong size now
  this->deleteSubspaceData();
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::sendDsgSizesWithGather(CommunicatorType comm,
                                                                      RankType collectorRank) {
  auto numSubspaces = static_cast<int>(this->getNumSubspaces());
  assert(numSubspaces > 0);
  assert(numSubspaces == subspacesDataSizes_.size());
  MPI_Datatype dtype = getMPIDatatype(abstraction::getabstractionDataType<SubspaceSizeType>());

  // perform gather (towards tl-manager)
  // send size of buffer to manager
  MPI_Gather(&numSubspaces, 1, MPI_INT, nullptr, 0, MPI_INT, collectorRank, comm);

  // send subspace sizes to manager
  MPI_Gatherv(subspacesDataSizes_.data(), numSubspaces, dtype, nullptr, nullptr, nullptr, dtype,
              collectorRank, comm);
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::receiveDsgSizesWithScatter(CommunicatorType comm,
                                                                          RankType collectorRank) {
  auto numSubspaces = static_cast<int>(this->getNumSubspaces());
  assert(numSubspaces > 0);
  assert(numSubspaces == subspacesDataSizes_.size());
  MPI_Datatype dtype = getMPIDatatype(abstraction::getabstractionDataType<SubspaceSizeType>());

  // receive updated sizes from manager
  MPI_Scatterv(nullptr, 0, nullptr, dtype, subspacesDataSizes_.data(), numSubspaces, dtype,
               collectorRank, comm);
  // assume that the sizes changed, the buffer might be the wrong size now
  this->deleteSubspaceData();
}

template <typename FG_ELEMENT>
inline void DistributedSparseGridUniform<FG_ELEMENT>::writeMinMaxCoefficents(
    const std::string& filename, size_t outputIndex) const {
  assert(this->isSubspaceDataCreated());
  bool writerProcess = false;
  std::ofstream ofs;
  if (this->rank_ == 0) {
    writerProcess = true;
    ofs = std::ofstream(filename + "_" + std::to_string(outputIndex) + ".txt");
    // std::cout << *this << std::endl;
  }
  // iterate subspaces
  assert(levels_.size() > 0);
  MPI_Datatype dataType = getMPIDatatype(abstraction::getabstractionDataType<combigrid::real>());
  auto realmax = std::numeric_limits<combigrid::real>::max();
  auto realmin = std::numeric_limits<combigrid::real>::min();
  auto smaller_real = [](const FG_ELEMENT& one, const FG_ELEMENT& two) {
    return std::real(one) < std::real(two);
  };

  for (SubspaceIndexType i = 0; i < static_cast<SubspaceIndexType>(levels_.size()); ++i) {
    auto minimumValue = realmax;
    auto maximumValue = realmin;
    if (subspacesDataSizes_[i] > 0) {
      // auto first = subspacesData_.begin();
      auto first = subspaces_[i];
      auto last = first + subspacesDataSizes_[i];
      auto it = std::min_element(first, last, smaller_real);
      minimumValue = std::real(*it);
      first = subspaces_[i];
      it = std::max_element(first, last, smaller_real);
      maximumValue = std::real(*it);
    }
    // allreduce the minimum and maximum values
    MPI_Allreduce(MPI_IN_PLACE, &minimumValue, 1, dataType, MPI_MIN, getCommunicator());
    MPI_Allreduce(MPI_IN_PLACE, &maximumValue, 1, dataType, MPI_MAX, getCommunicator());

    // if on zero process, write them out to file
    if (writerProcess) {
      const auto& level = getLevelVector(i);
      if (minimumValue < realmax)
        ofs << level << " : " << minimumValue << ", " << maximumValue << std::endl;
    }
  }
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::writeToDiskChunked(std::string filePrefix) {
  std::string myFilename = filePrefix + std::to_string(this->rank_);
  std::ofstream ofp(myFilename, std::ios::out | std::ios::binary);
  ofp.write(reinterpret_cast<const char*>(this->getRawData()), this->getRawDataSize() * sizeof(FG_ELEMENT));
  ofp.close();
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::readFromDiskChunked(std::string filePrefix){
  std::string myFilename = filePrefix + std::to_string(this->rank_);
  std::ifstream ifp(myFilename, std::ios::in | std::ios::binary);
  ifp.read(reinterpret_cast<char*>(this->getRawData()), this->getRawDataSize() * sizeof(FG_ELEMENT));
  ifp.close();
}

template <typename FG_ELEMENT>
bool DistributedSparseGridUniform<FG_ELEMENT>::writeOneFile(std::string fileName) const {
  auto comm = this->getCommunicator();

  MPI_Offset len = this->getRawDataSize();
  auto data = this->getRawData();
  bool success = mpiio::writeValuesConsecutive<FG_ELEMENT>(data, len, fileName, comm);
  return success;
}

template <typename FG_ELEMENT>
bool DistributedSparseGridUniform<FG_ELEMENT>::readOneFile(std::string fileName) {
  auto comm = this->getCommunicator();

  // get offset in file
  MPI_Offset len = this->getRawDataSize();
  auto data = this->getRawData();
  bool success = mpiio::readValuesConsecutive<FG_ELEMENT>(data, len, fileName, comm);
  return success;
}

template <typename FG_ELEMENT>
bool DistributedSparseGridUniform<FG_ELEMENT>::readOneFileAndReduce(std::string fileName,
                                                                    int numberOfChunks) {
  auto comm = this->getCommunicator();

  const int numElementsInChunk = this->getRawDataSize() / numberOfChunks;
  const int remainder = this->getRawDataSize() % numberOfChunks;
  const int numElementsToBuffer = numElementsInChunk + (remainder == 0 ? 0 : 1);

  // get offset in file
  const MPI_Offset len = this->getRawDataSize();
  auto data = this->getRawData();
  bool success = mpiio::readReduceValuesConsecutive<FG_ELEMENT>(
      data, len, fileName, comm, numElementsToBuffer, std::plus<FG_ELEMENT>{});

  return success;
}

template <typename FG_ELEMENT>
bool DistributedSparseGridUniform<FG_ELEMENT>::writeSubspaceSizesToFile(
    std::string fileName) const {
  auto comm = this->getCommunicator();
  MPI_Offset len = this->getNumSubspaces();
  bool success = mpiio::writeValuesConsecutive<SubspaceSizeType>(
      this->getSubspaceDataSizes().data(), len, fileName, comm);
  return success;
}

template <typename FG_ELEMENT>
bool DistributedSparseGridUniform<FG_ELEMENT>::readSubspaceSizesFromFile(
    std::string fileName, bool withCollectiveBuffering) {
  auto comm = this->getCommunicator();
  MPI_Offset len = this->getNumSubspaces();
  bool success = mpiio::readValuesConsecutive<SubspaceSizeType>(
      this->subspacesDataSizes_.data(), len, fileName, comm, withCollectiveBuffering);
  return success;
}

template <typename FG_ELEMENT>
template <typename ReduceFunctionType>
bool DistributedSparseGridUniform<FG_ELEMENT>::readReduceSubspaceSizesFromFile(
    std::string fileName, ReduceFunctionType reduceFunction, int numElementsToBuffer,
    bool withCollectiveBuffering) {
  auto comm = this->getCommunicator();
  MPI_Offset len = this->getNumSubspaces();
  if (numElementsToBuffer == 0) {
    numElementsToBuffer = len;
  }

  bool success = mpiio::readReduceValuesConsecutive<SubspaceSizeType>(
      this->subspacesDataSizes_.data(), len, fileName, comm, numElementsToBuffer, reduceFunction,
      withCollectiveBuffering);

  return success;
}

/**
 * Sends the raw dsg data to the destination process in communicator comm.
 */
template <typename FG_ELEMENT>
static void sendDsgData(DistributedSparseGridUniform<FG_ELEMENT>* dsgu, RankType dest,
                        CommunicatorType comm) {
  FG_ELEMENT* data = dsgu->getRawData();
  auto dataSize = dsgu->getRawDataSize();
  MPI_Datatype dataType = getMPIDatatype(abstraction::getabstractionDataType<FG_ELEMENT>());

  size_t sentRecvd = 0;
  while ((dataSize - sentRecvd) / INT_MAX > 0) {
    MPI_Send(data + sentRecvd, (int)INT_MAX, dataType, dest, TRANSFER_DSGU_DATA_TAG, comm);
    sentRecvd += INT_MAX;
  }
  MPI_Send(data + sentRecvd, (int)(dataSize - sentRecvd), dataType, dest, TRANSFER_DSGU_DATA_TAG,
           comm);
}

/**
* Recvs the raw dsg data from the source process in communicator comm.
*/
template <typename FG_ELEMENT>
static void recvDsgData(DistributedSparseGridUniform<FG_ELEMENT> * dsgu,
                          RankType source, CommunicatorType comm) {
  FG_ELEMENT* data = dsgu->getRawData();
  auto dataSize = dsgu->getRawDataSize();
  MPI_Datatype dataType = getMPIDatatype(abstraction::getabstractionDataType<FG_ELEMENT>());

  size_t sentRecvd = 0;
  while ((dataSize - sentRecvd) / INT_MAX > 0) {
    MPI_Recv(data + sentRecvd, (int)INT_MAX, dataType, source, TRANSFER_DSGU_DATA_TAG, comm,
             MPI_STATUS_IGNORE);
    sentRecvd += INT_MAX;
  }
  MPI_Recv(data + sentRecvd, (int)(dataSize - sentRecvd), dataType, source, TRANSFER_DSGU_DATA_TAG,
           comm, MPI_STATUS_IGNORE);
}

/**
 * Asynchronous Bcast of the raw dsg data in the communicator comm.
 */
template <typename FG_ELEMENT>
static MPI_Request asyncBcastDsgData(DistributedSparseGridUniform<FG_ELEMENT>* dsgu, RankType root,
                                     CommunicatorType comm) {
  if (dsgu->getRawDataSize() >= INT_MAX) {
    throw std::runtime_error(
        "asyncBcastDsgData: Dsg is too large and can not be "
        "transferred in a single MPI Call (not "
        "supported yet) try a more refined"
        "decomposition");
  }

  FG_ELEMENT* data = dsgu->getRawData();
  int dataSize  = static_cast<int>(dsgu->getRawDataSize());
  MPI_Datatype dataType = getMPIDatatype(abstraction::getabstractionDataType<FG_ELEMENT>());
  MPI_Request request = MPI_REQUEST_NULL;

  auto success = MPI_Ibcast(data, dataSize, dataType, root, comm, &request);
  assert(success == MPI_SUCCESS);
  return request;
}

/**
* Sends all subspace data sizes to the receiver in communicator comm.
*/
template <typename FG_ELEMENT>
static void sendSubspaceDataSizes(DistributedSparseGridUniform<FG_ELEMENT> * dsgu,
                          RankType dest, CommunicatorType comm) {
  assert(dsgu->getNumSubspaces() > 0);

  const std::vector<int>& subspacesDataSizes = dsgu->getSubspaceDataSizes();
  MPI_Send(subspacesDataSizes.data(), subspacesDataSizes.size(), MPI_INT, dest, TRANSFER_SUBSPACE_DATA_SIZES_TAG, comm);
}

/**
* Receives reduced subspace data sizes from the sender in communicator recvComm
* and concurrently distributes them inside bcastComm.
*/
template <typename FG_ELEMENT>
static MPI_Request recvAndBcastSubspaceDataSizes(DistributedSparseGridUniform<FG_ELEMENT> * dsgu,
                                         RankType recvSrc,
                                         CommunicatorType recvComm,
                                         RankType bcastRoot,
                                         CommunicatorType bcastComm) {
  assert(dsgu->getNumSubspaces() > 0);
  const std::vector<int>& subspacesDataSizes = dsgu->getSubspaceDataSizes();
  std::vector<int> buf(subspacesDataSizes.size());

  // receive subspace data sizes from manager
  MPI_Status status;
  MPI_Recv(buf.data(), buf.size(), MPI_INT, recvSrc, TRANSFER_SUBSPACE_DATA_SIZES_TAG, recvComm, &status);

  // distribute subspace sizes asynchronously
  MPI_Request request;
  MPI_Ibcast(buf.data(), buf.size(), MPI_INT, bcastRoot, bcastComm, &request);

  // update subspace data sizes of dsgu
  for (size_t i = 0; i < subspacesDataSizes.size(); i++) {
    dsgu->setDataSize(i, buf[i]);
  }
  return request;
}

} /* namespace combigrid */

#endif /* SRC_SGPP_COMBIGRID_SPARSEGRID_DISTRIBUTEDSPARSEGRID_HPP_ */
