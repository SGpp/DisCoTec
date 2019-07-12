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
#include "sgpp/distributedcombigrid/manager/ProcessGroupSignals.hpp"

#include <boost/serialization/vector.hpp>

using namespace combigrid;

/*
 * Instead of having private static functions, I put these functions in an
 * unnamed namespace. So, they are not accessible from outside the file, as well.
 * In the general case, this would have the advantage, that we can change
 * the declaration of these functions without touching the declaration of the
 * class. So we avoid recompilation of all files that use the class.
 */
namespace {

template <typename FG_ELEMENT>
struct SubspaceSGU {
  LevelVector level_;

  IndexVector sizes_;

  size_t dataSize_;

  std::vector<FG_ELEMENT> data_;

  SubspaceSGU& operator+=(const SubspaceSGU& rhs) 
  {
    assert(this->level_ == rhs.level_);
    assert(data_.size() == rhs.data_.size());

    for(size_t i = 0; i < data_.size(); ++i){
      this->data_[i] += rhs.data_[i];
    }

    return *this;
  }

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& level_;
    ar& sizes_;
    ar& dataSize_;
    ar& data_;
  }
};

/*
 template<typename FG_ELEMENT>
 bool mycomp( const SubspaceSGU<FG_ELEMENT>& ss1, const SubspaceSGU<FG_ELEMENT>& ss2 ){
 return ( ss1.dataSize_ > ss2.dataSize_ );
 }*/

}  // end anonymous namespace

namespace combigrid {

template <typename FG_ELEMENT>
class DistributedSparseGridUniform {
 public:
  /** create sparse grid of dimension d and specify for each dimension the
   *  the maximum discretization level and whether there is a boundary or not
   */
  DistributedSparseGridUniform(DimType dim, const LevelVector& lmax, const LevelVector& lmin,
                               const std::vector<bool>& boundary, CommunicatorType comm,
                               size_t procsPerNode = 0);
  DistributedSparseGridUniform(){}

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

  inline const LevelVector& getNMax() const;

  inline const LevelVector& getNMin() const;

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

  void recvAndAddDSGUniform(RankType src, CommunicatorType comm);

 private:
  void createLevels(DimType dim, const LevelVector& nmax, const LevelVector& lmin);

  void createLevelsRec(size_t dim, size_t n, size_t d, LevelVector& l, const LevelVector& nmax);

  void setSizes();

  DimType dim_;

  LevelVector nmax_;

  LevelVector lmin_;

  std::vector<LevelVector> levels_;

  std::vector<bool> boundary_;

  CommunicatorType comm_;

  std::vector<RankType> subspaceToProc_;

  RankType rank_;

  int commSize_;

  std::vector<SubspaceSGU<FG_ELEMENT> > subspaces_;

  friend class boost::serialization::access;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version);
};

}  // namespace

namespace combigrid {

// at construction create only levels, no data
template <typename FG_ELEMENT>
DistributedSparseGridUniform<FG_ELEMENT>::DistributedSparseGridUniform(
    DimType dim, const LevelVector& lmax, const LevelVector& lmin,
    const std::vector<bool>& boundary, CommunicatorType comm, size_t procsPerNode)
    : dim_(dim) {
  assert(dim > 0);

  assert(lmax.size() == dim);

  for (size_t i = 0; i < lmax.size(); ++i) assert(lmax[i] > 0);

  assert(lmin.size() == dim);

  for (size_t i = 0; i < lmin.size(); ++i) assert(lmin[i] > 0);

  assert(boundary.size() == dim);

  MPI_Comm_rank(comm, &rank_);
  MPI_Comm_size(comm, &commSize_);
  comm_ = comm;

  nmax_ = lmax;
  lmin_ = lmin;
  boundary_ = boundary;

  createLevels(dim, nmax_, lmin_);

  subspaces_.resize(levels_.size());

  for (size_t i = 0; i < levels_.size(); ++i) subspaces_[i].level_ = levels_[i];

  setSizes();
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::print(std::ostream& os) const {
  for (size_t i = 0; i < subspaces_.size(); ++i) {
    os << i << " " << subspaces_[i].level_ << " " << subspaces_[i].sizes_ << " "
       << subspaces_[i].dataSize_ << std::endl;
  }
}

template <typename FG_ELEMENT>
std::ostream& operator<<(std::ostream& os, const DistributedSparseGridUniform<FG_ELEMENT>& sg) {
  sg.print(os);
  return os;
}

template <typename FG_ELEMENT>
DistributedSparseGridUniform<FG_ELEMENT>::~DistributedSparseGridUniform() {}

// start recursion by setting dim=d=dimensionality of the vector space
template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::createLevelsRec(size_t dim, size_t n, size_t d,
                                                               LevelVector& l,
                                                               const LevelVector& nmax) {
  // sum rightmost entries of level vector
  LevelType lsum(0);

  for (size_t i = dim; i < l.size(); ++i) lsum += l[i];

  for (LevelType ldim = 1; ldim <= LevelType(n) + LevelType(d) - 1 - lsum; ++ldim) {
    l[dim - 1] = ldim;

    if (dim == 1) {
      if (l <= nmax) {
        levels_.push_back(l);
        // std::cout << l << std::endl;
      }
    } else {
      createLevelsRec(dim - 1, n, d, l, nmax);
    }
  }
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::createLevels(DimType dim, const LevelVector& nmax,
                                                            const LevelVector& lmin) {
  assert(nmax.size() == dim);
  assert(lmin.size() == dim);

  // compute c which fulfills nmax - c*1  >= lmin

  LevelVector ltmp(nmax);
  LevelType c = 0;

  while (ltmp > lmin) {
    ++c;

    for (size_t i = 0; i < dim; ++i) {
      ltmp[i] = nmax[i] - c;

      if (ltmp[i] < 1) ltmp[i] = 1;
    }
  }

  LevelVector rlmin(dim);

  for (size_t i = 0; i < rlmin.size(); ++i) {
    rlmin[i] = nmax[i] - c;
  }

  LevelType n = sum(rlmin) + c - dim + 1;

  LevelVector l(dim);
  createLevelsRec(dim, n, dim, l, nmax);
}

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

/* get index of space with l. returns -1 if not included */
template <typename FG_ELEMENT>
IndexType DistributedSparseGridUniform<FG_ELEMENT>::getIndex(const LevelVector& l) const {
  for (const auto& l_i : l){
    #ifdef DEBUG_OUTPUT
    // std::cerr << "getIndex()"<< std::endl;
    #endif
    assert(l_i > 0);
  }
  // std::cout << "get index of "<< toString(l) <<" before"<< std::endl;
  for (IndexType i = 0; i < IndexType(levels_.size()); ++i) {
    // std::cout << "...vs "<< toString(levels_[i]) << std::endl;
    if (levels_[i] == l) {
      return i;
    }
  }
  // assert (false && "space not found in levels_");
  return -1;
}

template <typename FG_ELEMENT>
inline const std::vector<bool>& DistributedSparseGridUniform<FG_ELEMENT>::getBoundaryVector()
    const {
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
    const LevelVector& l) {
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
inline const LevelVector& DistributedSparseGridUniform<FG_ELEMENT>::getNMax() const {
  return nmax_;
}

template <typename FG_ELEMENT>
inline const LevelVector& DistributedSparseGridUniform<FG_ELEMENT>::getNMin() const {
  return lmin_;
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


/**
* Root broadcasts all subspaces to the other ranks in communicator comm.
*
* TODO: For now, we broadcast all subspaces sequentially instead of sending the whole
* dsgu in a single call. (Sending the whole dsgu takes approx. double the communication volume)
* Sending the whole dsgu could be faster.
*/
template <typename FG_ELEMENT>
static void broadcastSubspaces(DistributedSparseGridUniform<FG_ELEMENT> * dsgu, 
                               const std::vector<LevelVector>& subspaces,
                               RankType root, CommunicatorType comm) {
  assert(dsgu->getNumSubspaces() > 0);
  assert(dsgu->getDataVector(dsgu->getNumSubspaces() - 1).size() >= 0);

  for (const LevelVector& ss : subspaces) {
    std::vector<FG_ELEMENT>& ssData = dsgu->getDataVector(ss);
    //assert(ssData.size() > 0 && "Subspace not initialized");
    MPI_Datatype dataType = getMPIDatatype(abstraction::getabstractionDataType<FG_ELEMENT>());
    MPI_Request request;
    MPI_Ibcast(ssData.data(), ssData.size(), dataType, static_cast<int>(root), comm, &request);
    MPI_Wait(&request, MPI_STATUS_IGNORE);
  }
}


/**
* Sends all subspaces to the receiver in communicator comm.
*/
template <typename FG_ELEMENT>
static void sendSubspaces(DistributedSparseGridUniform<FG_ELEMENT> * dsgu,
                          const std::vector<LevelVector>& subspaces,
                          RankType dest, CommunicatorType comm) {
  assert(dsgu->getNumSubspaces() > 0);
  assert(dsgu->getDataVector(dsgu->getNumSubspaces() - 1).size() >= 0);

  for (const LevelVector& ss : subspaces) {
    std::vector<FG_ELEMENT>& ssData = dsgu->getDataVector(ss);
    //assert(ssData.size() > 0 && "Subspace not initialized");

    MPI_Datatype dataType = getMPIDatatype(abstraction::getabstractionDataType<FG_ELEMENT>());
    MPI_Send(ssData.data(), ssData.size(), dataType, dest, sendSSTag, comm);
    //std::cout << "Worker sent ss with size " << ssData.size() << " to manager" << std::endl;
  }
}

/**
* Receives subspaces from the sender in communicator recvComm and concurrently
* distributes them inside bcastComm.
*/
template <typename FG_ELEMENT>
static std::vector<MPI_Request> recvAndBcastSubspaces(DistributedSparseGridUniform<FG_ELEMENT> * dsgu, 
                                                      const std::vector<LevelVector>& subspaces,
                                                      RankType recvSrc,
                                                      CommunicatorType recvComm,
                                                      RankType bcastRoot,
                                                      CommunicatorType bcastComm) {
  assert(dsgu->getNumSubspaces() > 0);
  assert(dsgu->getDataVector(dsgu->getNumSubspaces() - 1).size() >= 0);

  std::vector<MPI_Request> requests;
  requests.reserve(subspaces.size());

  for (const LevelVector& ss : subspaces) {
    std::vector<FG_ELEMENT>& ssData = dsgu->getDataVector(ss);
    //assert(ssData.size() > 0 && "Subspace not initialized");

    // receive subspace from manager
    MPI_Datatype dataType = getMPIDatatype(abstraction::getabstractionDataType<FG_ELEMENT>());
    MPI_Status status;
    MPI_Recv(ssData.data(), ssData.size(), dataType, recvSrc, sendSSTag, recvComm, &status);
    //std::cout << "Worker recieved ss with size " << ssData.size() << " from manager" << std::endl;

    // distribute subspace asynchronously
    MPI_Request request;
    MPI_Ibcast(ssData.data(), ssData.size(), dataType, bcastRoot, bcastComm, &request);
    requests.push_back(request);
  }
  return requests;
}

/**
* Sequentially adds all received subspaces to the dsgu and concurrently
* broadcasts them inside bcastComm.
*/
template <typename FG_ELEMENT>
static std::vector<MPI_Request> recvAndAddAndBcastSubspaces(DistributedSparseGridUniform<FG_ELEMENT> * dsgu,
                                                            const std::vector<LevelVector>& subspaces,
                                                            RankType recvSrc,
                                                            CommunicatorType recvComm,
                                                            RankType bcastRoot,
                                                            CommunicatorType bcastComm) {
  assert(dsgu->getNumSubspaces() > 0);
  assert(dsgu->getDataVector(dsgu->getNumSubspaces() - 1).size() >= 0);

  std::vector<MPI_Request> requests;
  requests.reserve(subspaces.size());

  for (const LevelVector& ss : subspaces) {
    std::vector<FG_ELEMENT>& ssData = dsgu->getDataVector(ss);
    std::unique_ptr<FG_ELEMENT[]> buf(new FG_ELEMENT[ssData.size()]);
    //assert(ssData.size() > 0 && "Subspace not initialized");

    // receive subspace from manager
    MPI_Datatype dataType = getMPIDatatype(abstraction::getabstractionDataType<FG_ELEMENT>());
    MPI_Status status;
    MPI_Recv(buf.get(), ssData.size(), dataType, recvSrc, sendSSTag, recvComm, &status);
    //std::cout << "Worker recieved ss with size " << ssData.size() << " from manager" << std::endl;

    // add subspace
    MPI_Reduce_local(buf.get(), ssData.data(), ssData.size(), dataType, MPI_SUM);
    buf.release();

    // distribute subspace asynchronously
    MPI_Request request;
    MPI_Ibcast(ssData.data(), ssData.size(), dataType, bcastRoot, bcastComm, &request);
    requests.push_back(request);
  }
  return requests;
}

} /* namespace combigrid */

#endif /* SRC_SGPP_COMBIGRID_SPARSEGRID_DISTRIBUTEDSPARSEGRID_HPP_ */
