#pragma once

#include <cassert>
#include <map>
#include <set>
#include <vector>

#include "../utils/Types.hpp"

namespace combigrid {

/**
 * @brief sparse grid base class
 *
 * this class is agnostic of the dimensionality and data type of the grid
 */
class AnyDistributedSparseGrid {
 public:
  // type used to index the subspaces
  // should be enough for the current scenario (cf. test_createTruncatedHierarchicalLevels_large)
  using SubspaceIndexType = int32_t;

  explicit AnyDistributedSparseGrid(size_t numSubspaces, CommunicatorType comm);

  virtual ~AnyDistributedSparseGrid();

  // cheap rule of 5
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  AnyDistributedSparseGrid() = delete;
  AnyDistributedSparseGrid(const AnyDistributedSparseGrid& other) = delete;
  AnyDistributedSparseGrid& operator=(const AnyDistributedSparseGrid&) = delete;
  AnyDistributedSparseGrid(AnyDistributedSparseGrid&& other) = delete;
  AnyDistributedSparseGrid& operator=(AnyDistributedSparseGrid&& other) = delete;
#endif  // DOXYGEN_SHOULD_SKIP_THIS

  void clearSubspaceCommunicators();

  /**
   * @brief get sum of all data sizes of all subspaces of this grid on this rank
   *
   * @return the number of numbers stored; multiply by sizeof(FG_ELEMENT) to get bytes
   */
  size_t getAccumulatedDataSize() const;

  /**
   * @brief get sum of data sizes of select subspaces of this grid on this rank
   *
   * @param subsetOfSubspaces the indices of the subspaces to consider
   * @return the number of numbers stored; multiply by sizeof(FG_ELEMENT) to get bytes
   */
  size_t getAccumulatedDataSize(const std::set<SubspaceIndexType>& subsetOfSubspaces) const;

  /**
   * @brief get the communicator of this grid
   *
   * usally the same as theMPISystem()->getLocalComm()
   */
  CommunicatorType getCommunicator() const;

  /**
   * @brief get the number of numbers in subspace i
   */
  SubspaceSizeType getDataSize(SubspaceIndexType i) const;
  SubspaceSizeType getDataSize(size_t i) const;

  /**
   * @brief get the subspaces that are only in this process group
   */
  std::set<typename AnyDistributedSparseGrid::SubspaceIndexType>& getIngroupSubspaces() const;

  /**
   * @brief get the number of subspaces
   */
  SubspaceIndexType getNumSubspaces() const;

  /**
   * @brief get the rank of this process in the communicator
   *
   * usually the same as theMPISystem()->getLocalRank()
   */
  RankType getRank() const;

  /**
   * @brief get the data sizes of all subspaces
   */
  const std::vector<SubspaceSizeType>& getSubspaceDataSizes() const;
  std::vector<SubspaceSizeType>& getSubspaceDataSizes();

  /**
   * @brief get the subspace communicators and the subspaces they are responsible for
   *
   * used for subspace reduce only
   */
  const std::vector<std::pair<CommunicatorType, std::vector<SubspaceIndexType>>>&
  getSubspacesByCommunicator() const;

  /**
   * @brief set the data size of a subspace
   */
  virtual void setDataSize(SubspaceIndexType i, SubspaceSizeType newSize);

  /**
   * @brief set the communicator for communication with the other process groups
   *
   * usually perpendicular to the local communicator: all ranks in this communicator have an
   * AnyDistributedSparseGrid that represents the same part of the domain
   *
   * used only for sparse grid reduce and outgroup reduce
   *
   * @param comm the communicator to use
   */
  void setOutgroupCommunicator(CommunicatorType comm, RankType rankInComm);

  /**
   * @brief generate the communicators for communication with the other process groups
   *
   * usually perpendicular to the local communicator: all ranks in this communicator have an
   * AnyDistributedSparseGrid that represents the same part of the domain
   *
   * used only for subspace reduce
   *
   * @param comm the communicator to generate sub-communicators from
   * @param rankInComm the rank in \p comm
   */
  void setSubspaceCommunicators(CommunicatorType comm, RankType rankInComm);

 protected:
  CommunicatorType comm_;

  RankType rank_;

  std::vector<SubspaceSizeType> subspacesDataSizes_;  // data sizes of all subspaces

  std::vector<std::pair<CommunicatorType, std::vector<SubspaceIndexType>>> subspacesByComm_;

  bool myOwnSubspaceCommunicators_ = false;
};

}  // namespace combigrid