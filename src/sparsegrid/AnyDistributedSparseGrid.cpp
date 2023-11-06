#include "sparsegrid/AnyDistributedSparseGrid.hpp"

#include <algorithm>  // std::find
#include <boost/multiprecision/cpp_int.hpp>
#include <cassert>
#include <iostream>
#include <map>
#include <numeric>  // std::accumulate

namespace combigrid {

AnyDistributedSparseGrid::AnyDistributedSparseGrid(size_t numSubspaces, CommunicatorType comm)
    : comm_(comm), subspacesDataSizes_(numSubspaces) {
  // make sure numSubspaces fits into SubspaceIndexType
  assert(numSubspaces <= static_cast<size_t>(std::numeric_limits<SubspaceIndexType>::max()));
  MPI_Comm_rank(comm_, &rank_);
}
size_t AnyDistributedSparseGrid::getAccumulatedDataSize() const {
  return std::accumulate(subspacesDataSizes_.begin(), subspacesDataSizes_.end(),
                         static_cast<size_t>(0));
}

size_t AnyDistributedSparseGrid::getAccumulatedDataSize(
    const std::set<SubspaceIndexType>& subsetOfSubspaces) const {
  // make filtered accumulation
  size_t accumulatedDataSize = 0;
  for (auto i : subsetOfSubspaces) {
    accumulatedDataSize += subspacesDataSizes_[i];
  }
  return accumulatedDataSize;
}

AnyDistributedSparseGrid::~AnyDistributedSparseGrid() { clearSubspaceCommunicators(); }

void AnyDistributedSparseGrid::clearSubspaceCommunicators() {
  if (this->myOwnSubspaceCommunicators_) {
    // free all my subspace communicators
    for (auto& pair : subspacesByComm_) {
      MPI_Comm_free(&pair.first);
    }
  }
  subspacesByComm_.clear();
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

SubspaceSizeType AnyDistributedSparseGrid::getDataSize(size_t i) const {
  return this->getDataSize(static_cast<SubspaceIndexType>(i));
}

std::set<typename AnyDistributedSparseGrid::SubspaceIndexType>&
AnyDistributedSparseGrid::getIngroupSubspaces() const {
  assert(this->getSubspacesByCommunicator().size() < 2);
  static thread_local std::set<SubspaceIndexType> ingroupSubspaces{};
  ingroupSubspaces.clear();
  if (this->getSubspacesByCommunicator().empty()) {
    for (SubspaceIndexType i = 0; i < this->getNumSubspaces(); ++i) {
      // select those with a positive data size
      if (this->getDataSize(i) > 0) {
        ingroupSubspaces.insert(i);
      }
    }
  } else {
    for (SubspaceIndexType i = 0; i < this->getNumSubspaces(); ++i) {
      // select those with a positive data size and that are not associated to outgroup communicator
      bool notForOutputComm = std::find(this->getSubspacesByCommunicator()[0].second.cbegin(),
                                        this->getSubspacesByCommunicator()[0].second.cend(),
                                        i) == this->getSubspacesByCommunicator()[0].second.cend();

      if (this->getDataSize(i) > 0 && notForOutputComm) {
        ingroupSubspaces.insert(i);
      }
    }
  }
  return ingroupSubspaces;
}

typename AnyDistributedSparseGrid::SubspaceIndexType AnyDistributedSparseGrid::getNumSubspaces()
    const {
  return static_cast<SubspaceIndexType>(subspacesDataSizes_.size());
}

RankType AnyDistributedSparseGrid::getRank() const { return rank_; }

const std::vector<SubspaceSizeType>& AnyDistributedSparseGrid::getSubspaceDataSizes() const {
  return subspacesDataSizes_;
}

std::vector<SubspaceSizeType>& AnyDistributedSparseGrid::getSubspaceDataSizes() {
  return subspacesDataSizes_;
}

const std::vector<
    std::pair<CommunicatorType, std::vector<typename AnyDistributedSparseGrid::SubspaceIndexType>>>&
AnyDistributedSparseGrid::getSubspacesByCommunicator() const {
  return subspacesByComm_;
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

using UIntForGroupReductionType = boost::multiprecision::uint256_t;

std::vector<UIntForGroupReductionType> getSubspaceVote(
    CommunicatorType comm, RankType rankInComm,
    const std::vector<SubspaceSizeType>& subspacesDataSizes) {
  RankType maxNumRanks = std::numeric_limits<UIntForGroupReductionType>::digits;
  if (rankInComm >= maxNumRanks) {
    throw std::runtime_error("Rank in communicator is too large; try fewer process groups? " +
                             std::to_string(maxNumRanks) + " allowed");
    // or implement things with bitsets or a veeery long string altogether
  }

  UIntForGroupReductionType mySummand = static_cast<UIntForGroupReductionType>(1)
                                               << rankInComm;

  // allocate vector of long long
  std::vector<UIntForGroupReductionType> subspaceVote(subspacesDataSizes.size(), 0);
  // set to mySummand if we have data in this subspace
  for (AnyDistributedSparseGrid::SubspaceIndexType i = 0; i < subspacesDataSizes.size(); ++i) {
    if (subspacesDataSizes[i] > 0) {
      subspaceVote[i] = mySummand;
    }
  }

  // vote by binary or
  MPI_Allreduce(MPI_IN_PLACE, subspaceVote.data(),
                static_cast<int>(subspaceVote.size()) * sizeof(UIntForGroupReductionType),
                MPI_CHAR, MPI_BOR, comm);
  return subspaceVote;
}

void AnyDistributedSparseGrid::setOutgroupCommunicator(CommunicatorType comm, RankType rankInComm) {
  assert(this->getNumSubspaces() > 0);
  assert(this->subspacesByComm_.size() < 2);
  // free previous communicator, if any
  clearSubspaceCommunicators();

  RankType maxNumRanks = std::numeric_limits<UIntForGroupReductionType>::digits;

  std::vector<UIntForGroupReductionType> subspaceVote =
      getSubspaceVote(comm, rankInComm, subspacesDataSizes_);

  // use this data to get the subspaces by the reduced values
  UIntForGroupReductionType manyGroups = 0;
  std::vector<SubspaceIndexType> subspacesForMany;
  for (SubspaceIndexType i = 0; i < this->getNumSubspaces(); ++i) {
    // only add if the vote has > 1 bit set (i.e. more than one rank has data in this subspace)
    if (subspaceVote[i] != 0 && (subspaceVote[i] & (subspaceVote[i] - 1)) != 0) {
      subspacesForMany.push_back(i);
      manyGroups |= subspaceVote[i];
    }
  }

  // use this data to set subspacesByComm_
  std::vector<RankType> ranks;
  // find all ranks that have voted
  for (RankType r = 0; r < maxNumRanks; ++r) {
    if (manyGroups & (static_cast<UIntForGroupReductionType>(1) << r)) {
      ranks.push_back(r);
    }
  }

  // assert that the given comm has the right size and all ranks are contributing
  int commSize;
  MPI_Comm_size(comm, &commSize);
  assert(commSize == static_cast<int>(ranks.size()) || subspacesForMany.empty());
  MPI_Comm subspaceComm = comm;
  myOwnSubspaceCommunicators_ = false;  // now we also need to reduce the data sizes (like for
                                        // sparse grid reduce, but only for the
                                        // subspaces to be exchanged)
  std::vector<SubspaceSizeType> subspaceDataSizesAlmostCopy(this->subspacesDataSizes_.size(), 0);
  for (const auto& subspace : subspacesForMany) {
    subspaceDataSizesAlmostCopy[subspace] = this->subspacesDataSizes_[subspace];
  }
  if (!subspacesForMany.empty()) {
    MPI_Allreduce(MPI_IN_PLACE, subspaceDataSizesAlmostCopy.data(),
                  subspaceDataSizesAlmostCopy.size(),
                  getMPIDatatype(abstraction::getabstractionDataType<SubspaceSizeType>()), MPI_MAX,
                  subspaceComm);
    for (const auto& subspace : subspacesForMany) {
      this->subspacesDataSizes_[subspace] = subspaceDataSizesAlmostCopy[subspace];
    }
  }
  if (std::find(ranks.begin(), ranks.end(), rankInComm) != ranks.end()) {
    subspacesByComm_.push_back(std::make_pair(subspaceComm, std::move(subspacesForMany)));
  } else {
    assert(subspacesForMany.empty());
  }
#ifndef NDEBUG
  // max-reduce the number of communicators created on each rank
  size_t maxNumComms = subspacesByComm_.size();
  MPI_Allreduce(MPI_IN_PLACE, &maxNumComms, 1,
                getMPIDatatype(abstraction::getabstractionDataType<size_t>()), MPI_MAX, comm);
  assert(maxNumComms <= 1);
#endif  // NDEBUG
  assert(subspacesByComm_.size() < 2);
}

void AnyDistributedSparseGrid::setSubspaceCommunicators(CommunicatorType comm,
                                                        RankType rankInComm) {
  clearSubspaceCommunicators();
  assert(this->getNumSubspaces() > 0);
  RankType maxNumRanks = std::numeric_limits<UIntForGroupReductionType>::digits;

  std::vector<UIntForGroupReductionType> subspaceVote =
      getSubspaceVote(comm, rankInComm, subspacesDataSizes_);

  // use this data to get the subspaces by the reduced values
  std::map<UIntForGroupReductionType, std::vector<SubspaceIndexType>> subspacesByVote;
  for (SubspaceIndexType i = 0; i < this->getNumSubspaces(); ++i) {
    // only add if the vote has > 1 bit set (i.e. more than one rank has data in this subspace)
    if (subspaceVote[i] != 0 && (subspaceVote[i] & (subspaceVote[i] - 1)) != 0) {
      subspacesByVote[subspaceVote[i]].push_back(i);
    }
  }

  MPI_Group wholeGroup;
  MPI_Comm_group(comm, &wholeGroup);
  // use this data to set subspacesByComm_
  for (const auto& kv : subspacesByVote) {
    // get a new group communicator that includes all that have voted
    // for this subspace
    UIntForGroupReductionType vote = kv.first;
    std::vector<RankType> ranks;
    // find all ranks that have voted
    for (RankType r = 0; r < maxNumRanks; ++r) {
      if (vote & (static_cast<UIntForGroupReductionType>(1) << r)) {
        ranks.push_back(r);
      }
    }
    // make those a group
    MPI_Group subspaceGroup;
    MPI_Group_incl(wholeGroup, int(ranks.size()), ranks.data(), &subspaceGroup);
    MPI_Comm subspaceComm;
    MPI_Comm_create(comm, subspaceGroup, &subspaceComm);
    MPI_Group_free(&subspaceGroup);
    // if I am one of the ranks, store the subspaces and the communicator
    if (std::find(ranks.begin(), ranks.end(), rankInComm) != ranks.end()) {
      subspacesByComm_.push_back(std::make_pair(subspaceComm, std::move(kv.second)));
    }
  }
  myOwnSubspaceCommunicators_ = true;
  MPI_Group_free(&wholeGroup);
  // max-reduce the number of communicators created on each rank
  size_t maxNumComms = subspacesByComm_.size();
  MPI_Allreduce(MPI_IN_PLACE, &maxNumComms, 1,
                getMPIDatatype(abstraction::getabstractionDataType<size_t>()), MPI_MAX, comm);
  if (rankInComm == 0 && this->getRank() == 0) {
    // TODO remove
    std::cout << "Found max. " << maxNumComms << " different communicators for subspaces"
              << std::endl;
  }
}

} /* namespace combigrid */
