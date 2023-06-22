#pragma once

#include <cassert>
#include <map>
#include <set>
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

  explicit AnyDistributedSparseGrid(size_t numSubspaces, CommunicatorType comm);

  virtual ~AnyDistributedSparseGrid();

  // cheap rule of 5
  AnyDistributedSparseGrid() = delete;
  AnyDistributedSparseGrid(const AnyDistributedSparseGrid& other) = delete;
  AnyDistributedSparseGrid& operator=(const AnyDistributedSparseGrid&) = delete;
  AnyDistributedSparseGrid(AnyDistributedSparseGrid&& other) = delete;
  AnyDistributedSparseGrid& operator=(AnyDistributedSparseGrid&& other) = delete;

  void clearSubspaceCommunicators();

  // sum of all data sizes of all subspaces
  size_t getAccumulatedDataSize() const;

  // sum of all data sizes of passed subspaces
  size_t getAccumulatedDataSize(const std::set<SubspaceIndexType>& subsetOfSubspaces) const;

  CommunicatorType getCommunicator() const;

  // data size of the subspace at index i
  SubspaceSizeType getDataSize(SubspaceIndexType i) const;

  // return the number of subspaces
  SubspaceIndexType getNumSubspaces() const;

  RankType getRank() const;

  // allows linear access to the data sizes of all subspaces
  const std::vector<SubspaceSizeType>& getSubspaceDataSizes() const;

  std::vector<SubspaceSizeType>& getSubspaceDataSizes();

  const std::vector<std::pair<CommunicatorType, std::vector<SubspaceIndexType>>>&
  getSubspacesByCommunicator() const;

  // sets data size of subspace with index i to newSize
  virtual void setDataSize(SubspaceIndexType i, SubspaceSizeType newSize);

  void setOutgroupCommunicator(CommunicatorType comm, RankType rankInComm);

  // sets the communicators for subspaces (required for subspace reduce)
  void setSubspaceCommunicators(CommunicatorType comm, RankType rankInComm);

 protected:
  CommunicatorType comm_;

  RankType rank_;

  std::vector<SubspaceSizeType> subspacesDataSizes_;  // data sizes of all subspaces

  std::vector<std::pair<CommunicatorType, std::vector<SubspaceIndexType>>> subspacesByComm_;

  bool myOwnSubspaceCommunicators_ = false;
};

}  // namespace combigrid