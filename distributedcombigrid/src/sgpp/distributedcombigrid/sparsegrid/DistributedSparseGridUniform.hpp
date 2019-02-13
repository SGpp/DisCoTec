/*
 * DistributedSparseGrid.h
 *
 *  Created on: Oct 19, 2015
 *      Author: heenemo
 */

#ifndef SRC_SGPP_COMBIGRID_SPARSEGRID_DISTRIBUTEDSPARSEGRIDUNIFORM_HPP_
#define SRC_SGPP_COMBIGRID_SPARSEGRID_DISTRIBUTEDSPARSEGRIDUNIFORM_HPP_

#include <assert.h>

#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/combischeme/DimAdaptiveCombiScheme.hpp"

using namespace combigrid;

/*
 * Instead of having private static functions, I put these functions in an
 * unnamed namespace. So, they are not accessible from outside the file, as well.
 * In the general case, this would have the advantage, that we can change
 * the declaration of these functions without touching the declaration of the
 * class. So we avoid recompilation of all files that use the class.
 */
namespace {

template<typename FG_ELEMENT>
struct SubspaceSGU {
  LevelVector level_;

  IndexVector sizes_;

  size_t dataSize_;

  std::vector<FG_ELEMENT> data_;
};

/*
 template<typename FG_ELEMENT>
 bool mycomp( const SubspaceSGU<FG_ELEMENT>& ss1, const SubspaceSGU<FG_ELEMENT>& ss2 ){
 return ( ss1.dataSize_ > ss2.dataSize_ );
 }*/

} // end anonymous namespace

namespace combigrid {

template<typename FG_ELEMENT>
class DistributedSparseGridUniform {
public:

  DistributedSparseGridUniform(const DimAdaptiveCombiScheme& scheme,
      const std::vector<bool>& boundary, CommunicatorType comm,
      size_t procsPerNode = 0);

  virtual ~DistributedSparseGridUniform();

  void print(std::ostream& os) const;

  // return level vector of subspace i
  inline const LevelVector& getLevelVector(size_t i) const;

  // return index of subspace i
  inline IndexType getIndex(const LevelVector& l) const;

  inline const std::vector<bool>& getBoundaryVector() const;

  // get pointer to first element in subspace with l
  inline FG_ELEMENT* getData(const LevelVector& l);

  // get pointer to first element in subspace i
  inline FG_ELEMENT* getData(size_t i);

  // get reference to data vector of subspace i.
  inline std::vector<FG_ELEMENT>& getDataVector(size_t i);

  // get reference to data vector of subspace with l.
  inline std::vector<FG_ELEMENT>& getDataVector(const LevelVector& l);

  inline size_t getDim() const;

  // return the number of subspaces
  inline size_t getNumSubspaces() const;

  // return the sizes for each dimension for subspace i
  inline const IndexVector& getSubspaceSizes(size_t i) const;

  // return the sizes for each dimension for subspace with l
  inline const IndexVector& getSubspaceSizes(const LevelVector& l) const;

  // return the number of elements of subspace i.
  // this number is independent of whether the subspace is initialized on this
  // process or not.
  inline size_t getSubspaceSize(size_t i) const;

  // return the number of elements of subspace i.
  // this number is independent of whether the subspace is initialized on this
  // process or not.
  inline size_t getSubspaceSize(const LevelVector& l) const;

  // check if a subspace with l is contained in the sparse grid
  // unlike getIndex this will not throw an assert in case l is not contained
  bool isContained(const LevelVector& l) const;

  inline size_t getDataSize(size_t i) const;

  inline size_t getDataSize(const LevelVector& l) const;

  inline CommunicatorType getCommunicator() const;

  inline int getCommunicatorSize() const;

private:
  void checkBasicInvariant();

  void setSizes();

  void calcProcAssignment(int procsPerNode);

  DimType dim_;

  std::vector<LevelVector> levels_;

  std::vector<bool> boundary_;

  CommunicatorType comm_;

  std::vector<RankType> subspaceToProc_;

  RankType rank_;

  int commSize_;

  std::vector<SubspaceSGU<FG_ELEMENT> > subspaces_;
};

// at construction create only levels, no data
template <typename FG_ELEMENT>
DistributedSparseGridUniform<FG_ELEMENT>::DistributedSparseGridUniform(
    const DimAdaptiveCombiScheme& scheme, const std::vector<bool>& boundary,
    CommunicatorType comm, size_t procsPerNode) : dim_(scheme.dim()) {

  MPI_Comm_rank(comm, &rank_);
  MPI_Comm_size(comm, &commSize_);
  comm_ = comm;

  boundary_ = boundary;

  levels_ = scheme.getLevels();

  subspaces_.resize(levels_.size());

  for (size_t i = 0; i < levels_.size(); ++i) subspaces_.at(i).level_ = levels_.at(i);

  setSizes();
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::print(std::ostream& os) const {
  for (size_t i = 0; i < subspaces_.size(); ++i) {
    os << i << " " << subspaces_[i].level_ << " " << subspaces_[i].sizes_ << " "
       << subspaces_[i].dataSize_ << std::endl;
  }
}


template<typename FG_ELEMENT>
std::ostream& operator<<(std::ostream& os, const DistributedSparseGridUniform<FG_ELEMENT>& sg) {
  sg.print(os);
  return os;
}

template <typename FG_ELEMENT>
DistributedSparseGridUniform<FG_ELEMENT>::~DistributedSparseGridUniform() {}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::setSizes() {
  for (size_t i = 0; i < subspaces_.size(); ++i) {
    IndexVector& sizes = subspaces_[i].sizes_;
    sizes.resize(dim_);
    const LevelVector& l = subspaces_[i].level_;

    // loop over all dimensions
    for (size_t j = 0; j < dim_; ++j) {
      if (l[j] == 1 && boundary_[j]) {
        sizes[j] = (size_t(std::pow(2.0, real(l[j] - 1))) + size_t(2));
      } else {
        sizes[j] = size_t(std::pow(2.0, real(l[j] - 1)));
      }
    }

    IndexType tmp(1);

    for (auto s : sizes) tmp *= s;

    subspaces_[i].dataSize_ = size_t(tmp);
  }
}

template <typename FG_ELEMENT>
inline const LevelVector& DistributedSparseGridUniform<FG_ELEMENT>::getLevelVector(size_t i) const {
  return levels_[i];
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::calcProcAssignment(int procsPerNode) {
  if (procsPerNode == 0) {
    for (size_t i = 0; i < subspaces_.size(); ++i) subspaceToProc_[i] = int(i) % commSize_;
  }
  /*
     else{
     // check if commsize a multiple of procs per node
     assert( (commSize_ % procsPerNode) == 0
     && "number of procs in comm must be multiple of procsPerNode" );
     int numNodes = commSize_ / procsPerNode;
     for( int i=0; i<int( subspaces_.size() ); ++i ){
     int nodeID = i % numNodes;
     int procInNodeID = ( i / numNodes ) % procsPerNode;
     subspaceToProc_[i] = nodeID * procsPerNode + procInNodeID;
     }
     }*/

  // use this to assign subspaces to first proc in group
  else {
    // todo: ceil
    int numNodes = static_cast<int>(std::ceil(real(commSize_) / real(procsPerNode)));

    for (int i = 0; i < int(subspaces_.size()); ++i) {
      int nodeID = i % numNodes;
      subspaceToProc_[i] = nodeID * procsPerNode;
    }
  }
}

/* get index of space with l. returns -1 if not included */
template <typename FG_ELEMENT>
IndexType DistributedSparseGridUniform<FG_ELEMENT>::getIndex(const LevelVector& l) const {
  // get index of l
  for (IndexType i = 0; i < IndexType(levels_.size()); ++i) {
    if (levels_[i] == l) {
      return i;
    }
  }

  return -1;
}

template <typename FG_ELEMENT>
inline const std::vector<bool>& DistributedSparseGridUniform<FG_ELEMENT>::getBoundaryVector() const {
  return boundary_;
}

template <typename FG_ELEMENT>
inline FG_ELEMENT* DistributedSparseGridUniform<FG_ELEMENT>::getData(const LevelVector& l) {
  IndexType i = getIndex(l);

  if (i < 0) {
    std::cout << "l = " << l << " not included in distributed sparse grid" << std::endl;
    assert(false);
  }

  return &subspaces_[i].data_[0];
}

template <typename FG_ELEMENT>
inline FG_ELEMENT* DistributedSparseGridUniform<FG_ELEMENT>::getData(size_t i) {
  return &subspaces_[i].data_[0];
}

template <typename FG_ELEMENT>
inline std::vector<FG_ELEMENT>& DistributedSparseGridUniform<FG_ELEMENT>::getDataVector(size_t i) {
  return subspaces_[i].data_;
}

template <typename FG_ELEMENT>
inline std::vector<FG_ELEMENT>& DistributedSparseGridUniform<FG_ELEMENT>::getDataVector(
    const LevelVector& l){
  IndexType i = getIndex(l);

  if (i < 0) {
    std::cout << "l = " << l << " not included in distributed sparse grid" << std::endl;
    assert(false);
  }

  return subspaces_[i].data_;
}

template <typename FG_ELEMENT>
inline DimType DistributedSparseGridUniform<FG_ELEMENT>::getDim() const {
  return dim_;
}

template <typename FG_ELEMENT>
inline size_t DistributedSparseGridUniform<FG_ELEMENT>::getNumSubspaces() const {
  return subspaces_.size();
}

template <typename FG_ELEMENT>
bool DistributedSparseGridUniform<FG_ELEMENT>::isContained(const LevelVector& l) const {
  // get index of l
  bool found = false;

  for (size_t i = 0; i < levels_.size(); ++i) {
    if (levels_[i] == l) {
      found = true;
      break;
    }
  }

  return found;
}

template <typename FG_ELEMENT>
const IndexVector& DistributedSparseGridUniform<FG_ELEMENT>::getSubspaceSizes(size_t i) const {
  return subspaces_[i].sizes_;
}

template <typename FG_ELEMENT>
const IndexVector& DistributedSparseGridUniform<FG_ELEMENT>::getSubspaceSizes(
    const LevelVector& l) const {
  IndexType i = getIndex(l);

  if (i < 0) {
    std::cout << "l = " << l << " not included in distributed sparse grid" << std::endl;
    assert(false);
  }

  return subspaces_[i].sizes_;
}

template <typename FG_ELEMENT>
size_t DistributedSparseGridUniform<FG_ELEMENT>::getSubspaceSize(size_t i) const {
  return subspaces_[i].dataSize_;
}

template <typename FG_ELEMENT>
size_t DistributedSparseGridUniform<FG_ELEMENT>::getSubspaceSize(const LevelVector& l) const {
  IndexType i = getIndex(l);

  if (i < 0) {
    std::cout << "l = " << l << " not included in distributed sparse grid" << std::endl;
    assert(false);
  }

  return subspaces_[i].dataSize_;
}

template <typename FG_ELEMENT>
size_t DistributedSparseGridUniform<FG_ELEMENT>::getDataSize(size_t i) const {
  return subspaces_[i].data_.size();
}

template <typename FG_ELEMENT>
size_t DistributedSparseGridUniform<FG_ELEMENT>::getDataSize(const LevelVector& l) const {
  IndexType i = getIndex(l);

  if (i < 0) {
    std::cout << "l = " << l << " not included in distributed sparse grid" << std::endl;
    assert(false);
  }

  return subspaces_[i].data_.size();
}

template <typename FG_ELEMENT>
CommunicatorType DistributedSparseGridUniform<FG_ELEMENT>::getCommunicator() const {
  return comm_;
}

template <typename FG_ELEMENT>
inline int DistributedSparseGridUniform<FG_ELEMENT>::getCommunicatorSize() const {
  return commSize_;
}

} /* namespace combigrid */

#endif /* SRC_SGPP_COMBIGRID_SPARSEGRID_DISTRIBUTEDSPARSEGRID_HPP_ */
