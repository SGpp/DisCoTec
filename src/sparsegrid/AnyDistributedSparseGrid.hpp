#pragma once

#include <cassert>
#include <vector>

#include "utils/Types.hpp"

namespace combigrid {

/* sparse grid base class; agnostic of the dimensionality and
 * data type used*/
class AnyDistributedSparseGrid {
 public:
  // type used to index the subspaces
  // should be enough for the current scenario (cf. test_createTruncatedHierarchicalLevels_large)
  using SubspaceIndexType = int32_t;

  explicit inline AnyDistributedSparseGrid(size_t numSubspaces, CommunicatorType comm);

  virtual ~AnyDistributedSparseGrid() = default;

  // cheap rule of 5
  AnyDistributedSparseGrid() = delete;
  AnyDistributedSparseGrid(const AnyDistributedSparseGrid& other) = delete;
  AnyDistributedSparseGrid& operator=(const AnyDistributedSparseGrid&) = delete;
  AnyDistributedSparseGrid(AnyDistributedSparseGrid&& other) = delete;
  AnyDistributedSparseGrid& operator=(AnyDistributedSparseGrid&& other) = delete;

  // sum of all data sizes of all subspaces
  inline size_t getAccumulatedDataSize() const;

  inline CommunicatorType getCommunicator() const;

  // data size of the subspace at index i
  inline SubspaceSizeType getDataSize(SubspaceIndexType i) const;

  // return the number of subspaces
  inline SubspaceIndexType getNumSubspaces() const;

  inline RankType getRank() const;

  // allows linear access to the data sizes of all subspaces
  inline const std::vector<SubspaceSizeType>& getSubspaceDataSizes() const;

  inline std::vector<SubspaceSizeType>& getSubspaceDataSizes();

  // sets data size of subspace with index i to newSize
  inline virtual void setDataSize(SubspaceIndexType i, SubspaceSizeType newSize);

 protected:
  CommunicatorType comm_;

  RankType rank_;

  std::vector<SubspaceSizeType> subspacesDataSizes_;  // data sizes of all subspaces
};

}  // namespace combigrid

namespace combigrid {

inline AnyDistributedSparseGrid::AnyDistributedSparseGrid(size_t numSubspaces,
                                                          CommunicatorType comm)
    : comm_(comm), subspacesDataSizes_(numSubspaces) {
  // make sure numSubspaces fits into SubspaceIndexType
  assert(numSubspaces <= std::numeric_limits<SubspaceIndexType>::max());
  MPI_Comm_rank(comm_, &rank_);
}
size_t AnyDistributedSparseGrid::getAccumulatedDataSize() const {
  return std::accumulate(subspacesDataSizes_.begin(), subspacesDataSizes_.end(),
                         static_cast<size_t>(0));
}

CommunicatorType AnyDistributedSparseGrid::getCommunicator() const { return comm_; }

SubspaceSizeType AnyDistributedSparseGrid::getDataSize(SubspaceIndexType i) const {
#ifndef NDEBUG
  if (i >= getNumSubspaces()) {
    std::cout << "Index too large, no subspace with this index included in distributed sparse grid"
              << std::endl;
    assert(false);
  }
#endif  // NDEBUG

  return subspacesDataSizes_[i];
}

inline typename AnyDistributedSparseGrid::SubspaceIndexType
AnyDistributedSparseGrid::getNumSubspaces() const {
  return static_cast<SubspaceIndexType>(subspacesDataSizes_.size());
}

inline RankType AnyDistributedSparseGrid::getRank() const { return rank_; }

inline const std::vector<SubspaceSizeType>& AnyDistributedSparseGrid::getSubspaceDataSizes() const {
  return subspacesDataSizes_;
}

inline std::vector<SubspaceSizeType>& AnyDistributedSparseGrid::getSubspaceDataSizes() {
  return subspacesDataSizes_;
}

void AnyDistributedSparseGrid::setDataSize(SubspaceIndexType i, SubspaceSizeType newSize) {
#ifndef NDEBUG
  assert(subspacesDataSizes_[i] == 0 || subspacesDataSizes_[i] == newSize);
  if (i >= getNumSubspaces()) {
    std::cout << "Index too large, no subspace with this index included in distributed sparse grid"
              << std::endl;
    assert(false);
  }
#endif  // NDEBUG
  subspacesDataSizes_[i] = newSize;
}

} /* namespace combigrid */
