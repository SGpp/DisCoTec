#ifndef SRC_SGPP_COMBIGRID_SPARSEGRID_DISTRIBUTEDSPARSEGRIDUNIFORM_HPP_
#define SRC_SGPP_COMBIGRID_SPARSEGRID_DISTRIBUTEDSPARSEGRIDUNIFORM_HPP_

#include <assert.h>

#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupSignals.hpp"
#include "sgpp/distributedcombigrid/mpi/MPITags.hpp"
#include <numeric>

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
  // commented most members, since they were virtually unused
  // if needed in the future, they may be used again

  // combigrid::LevelVector level_; // level of the subspace

  // combigrid::IndexVector sizes_; // contains the number of points per dim of the whole ss

  // size_t size_; // contains the number of Points of the whole ss

  FG_ELEMENT * data_; // contains the values at the data points (Attention: Due to the decomposition, only part of the full ss may be stored)

  // size_t dataSize_; // contains the number of values stored in data_ == size of the ss part.
};

}  // end anonymous namespace

namespace combigrid {

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
  /** create sparse grid of dimension d and specify for each dimension the
   * maximum discretization level and whether there is a boundary or not
   * No data is allocated and the sizes of the subspace data is initialized to 0.
   */
  DistributedSparseGridUniform(DimType dim, const LevelVector& lmax, const LevelVector& lmin,
                               const std::vector<bool>& boundary, CommunicatorType comm,
                               size_t procsPerNode = 0);

  /**
   * create an empty (no data) sparse grid with given subspaces.
   */
  DistributedSparseGridUniform(DimType dim,
                               const std::vector<LevelVector>& subspaces,
                               const std::vector<bool>& boundary,
                               CommunicatorType comm,
                               size_t procsPerNode = 0);

  DistributedSparseGridUniform(){}

  virtual ~DistributedSparseGridUniform();

  void print(std::ostream& os) const;

  // allocates memory for subspace data and sets pointers to subspaces
  void createSubspaceData();

  // deletes memory for subspace data and invalids pointers to subspaces
  void deleteSubspaceData();

  // creates data if necessary and sets all data elements to zero
  void setZero();

  // return level vector of subspace i
  inline const LevelVector& getLevelVector(size_t i) const;

  // return index of subspace i
  inline IndexType getIndex(const LevelVector& l) const;

  inline const std::vector<bool>& getBoundaryVector() const;

  // returns a pointer to first element in subspace with l
  inline FG_ELEMENT* getData(const LevelVector& l);

  // returns a pointer to first element in subspace i
  inline FG_ELEMENT* getData(size_t i);

  // allows a linear access to the whole subspace data stored in this dsg
  inline FG_ELEMENT* getRawData();

  inline size_t getDim() const;

  // return the number of subspaces
  inline size_t getNumSubspaces() const;

  // return the sizes for each dimension for subspace i
  inline const IndexVector& getSubspaceSizes(size_t i) const;

  // return the sizes for each dimension for subspace with l
  inline const IndexVector& getSubspaceSizes(const LevelVector& l) const;

  // // return the number of elements of subspace i.
  // // this number is independent of whether the subspace is initialized on this
  // // process or not.
  // inline size_t getSubspaceSize(size_t i) const;

  // // return the number of elements of subspace i.
  // // this number is independent of whether the subspace is initialized on this
  // // process or not.
  // inline size_t getSubspaceSize(const LevelVector& l) const;

  // check if a subspace with l is contained in the sparse grid
  // unlike getIndex this will not throw an assert in case l is not contained
  bool isContained(const LevelVector& l) const;

  // data size of the subspace at index i
  inline size_t getDataSize(size_t i) const;

  // data size of the subspace with level l
  inline size_t getDataSize(const LevelVector& l) const;

  // sets data size of subspace with index i to newSize
  inline void setDataSize(size_t i, size_t newSize);

  // sets data size of subspace with level l to newSize
  inline void setDataSize(const LevelVector& l, size_t newSize);

  // returns the number of allocated grid points == size of the raw data vector
  inline size_t getRawDataSize() const;

  inline CommunicatorType getCommunicator() const;

  inline int getCommunicatorSize() const;

  void zeroOutSubspaces();

  // allows linear access to the data sizes of all subspaces
  const std::vector<size_t>& getSubspaceDataSizes() const;

  // returns true if data for the subspaces has been created
  bool isSubspaceDataCreated() const;

  void writeMinMaxCoefficents(const std::string& filename, size_t i) const;

 private:
  std::vector<LevelVector> createLevels(DimType dim, const LevelVector& nmax, const LevelVector& lmin);

  void createLevelsRec(size_t dim, size_t n, size_t d, LevelVector& l, const LevelVector& nmax, std::vector<LevelVector>& created) const;

  // void setSizes();

  DimType dim_;

  std::vector<LevelVector> levels_; // linear access to all subspaces

  std::vector<bool> boundary_;

  CommunicatorType comm_;

  std::vector<RankType> subspaceToProc_;

  RankType rank_;

  int commSize_;

  std::vector<SubspaceSGU<FG_ELEMENT> > subspaces_; // subspaces of the dsg

  std::vector<FG_ELEMENT> subspacesData_; // allows linear access to all subspaces data

  std::vector<size_t> subspacesDataSizes_; // allocated data sizes of all subspaces

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
    : DistributedSparseGridUniform(dim,createLevels(dim, lmax, lmin),boundary, comm, procsPerNode)
    {}

// at construction create only levels, no data
template <typename FG_ELEMENT>
DistributedSparseGridUniform<FG_ELEMENT>::DistributedSparseGridUniform(
    DimType dim,
    const std::vector<LevelVector>& subspaces,
    const std::vector<bool>& boundary,
    CommunicatorType comm,
    size_t procsPerNode /*= 0*/)
    : boundary_(boundary),
      comm_(comm),
      dim_(dim),
      levels_(subspaces),
      subspacesDataSizes_(subspaces.size())
{
  assert(dim > 0);
  assert(boundary.size() == dim);

  for (auto& l : subspaces){
    assert(l.size() == dim);
    for (size_t i = 0; i < l.size(); ++i) {
      assert(l[i] > 0);
    }
  }

  // init subspaces
  subspaces_.reserve(subspaces.size());
  for (size_t i = 0; i < subspaces.size(); i++) {
   subspaces_.push_back({ ///*.level_    =*/ subspaces[i],
                          // /*.sizes_    =*/ IndexVector(0),
                          // /*.size_     =*/ 0,
                          /*.data_     =*/ nullptr
                          });
  }
  // setSizes();

  MPI_Comm_rank(comm_, &rank_);
  MPI_Comm_size(comm_, &commSize_);
}

template <typename FG_ELEMENT>
bool DistributedSparseGridUniform<FG_ELEMENT>::isSubspaceDataCreated() const {
  return subspacesData_.size() != 0;
}

/** Zero initializes the dsgu data in case no data is already present.
 * Otherwise nothing happens.
 */
template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::createSubspaceData() {
  if (not isSubspaceDataCreated()) {
    size_t numDataPoints = std::accumulate(subspacesDataSizes_.begin(), subspacesDataSizes_.end(), 0);
    assert(numDataPoints > 0 && "all subspaces in dsg have 0 size");

    if (numDataPoints > 0) {
      subspacesData_ = std::vector<FG_ELEMENT>(numDataPoints);

      // update pointers and sizes in subspaces
      size_t offset = 0;
      for (size_t i = 0; i < subspaces_.size(); i++) {
        subspaces_[i].data_ = subspacesData_.data() + offset;
        offset += subspacesDataSizes_[i];
      }
    }
  }
}

/** Deallocates the dsgu data.
 *  This affects the values stored at the grid points, pointers which address
 *  the subspaces data and the subspaces data sizes. No meta data which
 *  characterizes the subspaces is affected.
 */
template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::deleteSubspaceData() {
  // assert(false &&
  //        "due to the way that DFGs register the DSG only once (currently), you should think of a "
  //        "way to reset the localFGIndexToLocalSGPointerList_ member, since it will become "
  //        "invalidated by the next line, or to have it use the subspaces_[*].data_ structures");
  if (isSubspaceDataCreated()) {
    subspacesData_ = std::vector<FG_ELEMENT>();

    // update pointers in subspaces
    for (auto& ss : subspaces_) {
      ss.data_ = nullptr;
    }
  }
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::setZero() {
  if (isSubspaceDataCreated())
    std::fill(subspacesData_.begin(), subspacesData_.end(), 0.);
  else
    createSubspaceData();
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::print(std::ostream& os) const {
  for (size_t i = 0; i < subspaces_.size(); ++i) {
    os << i << " " << levels_[i] << " " << subspacesDataSizes_[i] << " "
      //  << subspaces_[i].size_
       << std::endl;
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
                                                               const LevelVector& nmax, std::vector<LevelVector>& created) const {
  // sum rightmost entries of level vector
  LevelType lsum(0);

  for (size_t i = dim; i < l.size(); ++i) lsum += l[i];

  for (LevelType ldim = 1; ldim <= LevelType(n) + LevelType(d) - 1 - lsum; ++ldim) {
    l[dim - 1] = ldim;

    if (dim == 1) {
      if (l <= nmax) {
        created.push_back(l);
        // std::cout << l << std::endl;
      }
    } else {
      createLevelsRec(dim - 1, n, d, l, nmax, created);
    }
  }
}

template <typename FG_ELEMENT>
std::vector<LevelVector> DistributedSparseGridUniform<FG_ELEMENT>::createLevels(DimType dim, const LevelVector& nmax,
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
  std::vector<LevelVector> created {};
  createLevelsRec(dim, n, dim, l, nmax, created);
  return created;
}

// template <typename FG_ELEMENT>
// void DistributedSparseGridUniform<FG_ELEMENT>::setSizes() {
//   for (size_t i = 0; i < subspaces_.size(); ++i) {
//     // IndexVector& sizes = subspaces_[i].sizes_;
//     // sizes.resize(dim_);
//     IndexVector sizes(dim_);
//     // const LevelVector& l = subspaces_[i].level_;
//     const LevelVector& l = levels_[i];

//     // loop over all dimensions
//     for (size_t j = 0; j < dim_; ++j) {
//       if (l[j] == 1 && boundary_[j]) {
//         sizes[j] = (size_t(std::pow(2.0, real(l[j] - 1))) + size_t(2));
//       } else {
//         sizes[j] = size_t(std::pow(2.0, real(l[j] - 1)));
//       }
//     }

//     IndexType tmp(1);

//     for (auto s : sizes) tmp *= s;

//     subspaces_[i].size_ = size_t(tmp);
//   }
// }

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::zeroOutSubspaces() {
  #pragma omp parallel for schedule(dynamic,4)
  for (size_t i = 0; i < subspaces_.size(); ++i) {
    for (size_t j = 0; j < subspaces_[i].data_.size(); j++)
    {
      subspaces_[i].data_[j]=0;
    }
    
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
  createSubspaceData();
  IndexType i = getIndex(l);

  if (i < 0) {
    std::cout << "l = " << l << " not included in distributed sparse grid" << std::endl;
    assert(false);
  }

  return subspaces_[i].data_;
}

template <typename FG_ELEMENT>
inline FG_ELEMENT* DistributedSparseGridUniform<FG_ELEMENT>::getData(size_t i) {
  createSubspaceData();
  return subspaces_[i].data_;
}

template <typename FG_ELEMENT>
inline FG_ELEMENT* DistributedSparseGridUniform<FG_ELEMENT>::getRawData() {
  createSubspaceData();
  return subspacesData_.data();
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

// template <typename FG_ELEMENT>
// const IndexVector& DistributedSparseGridUniform<FG_ELEMENT>::getSubspaceSizes(size_t i) const {
//   return subspaces_[i].sizes_;
// }

// template <typename FG_ELEMENT>
// const IndexVector& DistributedSparseGridUniform<FG_ELEMENT>::getSubspaceSizes(
//     const LevelVector& l) const {
//   IndexType i = getIndex(l);

//   if (i < 0) {
//     std::cout << "l = " << l << " not included in distributed sparse grid" << std::endl;
//     assert(false);
//   }

//   return subspaces_[i].sizes_;
// }

// template <typename FG_ELEMENT>
// size_t DistributedSparseGridUniform<FG_ELEMENT>::getSubspaceSize(size_t i) const {
//   return subspaces_[i].size_;
// }

// template <typename FG_ELEMENT>
// size_t DistributedSparseGridUniform<FG_ELEMENT>::getSubspaceSize(const LevelVector& l) const {
//   IndexType i = getIndex(l);

//   if (i < 0) {
//     std::cout << "l = " << l << " not included in distributed sparse grid" << std::endl;
//     assert(false);
//   }

//   return subspaces_[i].size_;
// }

template <typename FG_ELEMENT>
size_t DistributedSparseGridUniform<FG_ELEMENT>::getDataSize(size_t i) const {
  if (i >= getNumSubspaces()) {
    std::cout << "Index too large, no subspace with this index included in distributed sparse grid" << std::endl;
    assert(false);
  }

  return subspacesDataSizes_[i];
}

template <typename FG_ELEMENT>
size_t DistributedSparseGridUniform<FG_ELEMENT>::getDataSize(const LevelVector& l) const {
  IndexType i = getIndex(l);

  if (i < 0) {
    std::cout << "l = " << l << " not included in distributed sparse grid" << std::endl;
    assert(false);
  }

  return subspacesDataSizes_[i];
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::setDataSize(size_t i, size_t newSize) {
  if (i >= getNumSubspaces()) {
    std::cout << "Index too large, no subspace with this index included in distributed sparse grid" << std::endl;
    assert(false);
  }

  subspacesDataSizes_[i] = newSize;
}

template <typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::setDataSize(const LevelVector& l, size_t newSize) {
  IndexType i = getIndex(l);

  if (i < 0) {
    std::cout << "l = " << l << " not included in distributed sparse grid" << std::endl;
    assert(false);
  }

  subspacesDataSizes_[i] = newSize;
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
inline const std::vector<size_t>& DistributedSparseGridUniform<FG_ELEMENT>::getSubspaceDataSizes() const {
  return subspacesDataSizes_;
}

template <typename FG_ELEMENT>
inline void DistributedSparseGridUniform<FG_ELEMENT>::writeMinMaxCoefficents(
    const std::string& filename, size_t i) const {
  bool writerProcess = false;
  std::ofstream ofs;
  MASTER_EXCLUSIVE_SECTION {
    writerProcess = true;
    ofs = std::ofstream(filename + "_" + std::to_string(i) + ".txt");
    // std::cout << *this << std::endl;
  }
  // iterate subspaces
  assert(levels_.size() > 0);
  MPI_Datatype dataType = getMPIDatatype(abstraction::getabstractionDataType<FG_ELEMENT>());
  for (size_t i = 0; i < levels_.size(); ++i) {
    // auto first = subspacesData_.begin();
    auto first = subspaces_[i].data_;
    auto last = first + subspacesDataSizes_[i];
    // allreduce the minimum and maxiumum values
    auto it = std::min_element(first, last);
    auto minimumValue = *it;
    MPI_Allreduce(MPI_IN_PLACE, &minimumValue, 1, dataType, MPI_MIN, getCommunicator());
    first = subspaces_[i].data_;
    it = std::max_element(first, last);
    auto maximumValue = *it;
    MPI_Allreduce(MPI_IN_PLACE, &maximumValue, 1, dataType, MPI_MAX, getCommunicator());

    // if on zero process, write them out to file
    if (writerProcess) {
      const auto& level = getLevelVector(i);
      ofs << level << " : " << minimumValue << ", " << maximumValue << std::endl;
    }
  }
}

/**
* Sends the raw dsg data to the destination process in communicator comm.
*/
template <typename FG_ELEMENT>
static void sendDsgData(DistributedSparseGridUniform<FG_ELEMENT> * dsgu,
                          RankType dest, CommunicatorType comm) {
  assert(dsgu->getRawDataSize() < INT_MAX && "Dsg is too large and can not be "
                                            "transferred in a single MPI Call (not "
                                            "supported yet) try a more refined"
                                            "decomposition");

  FG_ELEMENT* data = dsgu->getRawData();
  int dataSize  = static_cast<int>(dsgu->getRawDataSize());
  MPI_Datatype dataType = getMPIDatatype(abstraction::getabstractionDataType<FG_ELEMENT>());

  MPI_Send(data, dataSize, dataType, dest, TRANSFER_DSGU_DATA_TAG, comm);
}

/**
* Recvs the raw dsg data from the source process in communicator comm.
*/
template <typename FG_ELEMENT>
static void recvDsgData(DistributedSparseGridUniform<FG_ELEMENT> * dsgu,
                          RankType source, CommunicatorType comm) {
  assert(dsgu->getRawDataSize() < INT_MAX && "Dsg is too large and can not be "
                                            "transferred in a single MPI Call (not "
                                            "supported yet) try a more refined"
                                            "decomposition");

  FG_ELEMENT* data = dsgu->getRawData();
  int dataSize  = static_cast<int>(dsgu->getRawDataSize());
  MPI_Datatype dataType = getMPIDatatype(abstraction::getabstractionDataType<FG_ELEMENT>());

  MPI_Recv(data, dataSize, dataType, source, TRANSFER_DSGU_DATA_TAG, comm, MPI_STATUS_IGNORE);
}

/**
 * Asynchronous Bcast of the raw dsg data in the communicator comm.
 */
template <typename FG_ELEMENT>
static MPI_Request asyncBcastDsgData(DistributedSparseGridUniform<FG_ELEMENT> * dsgu,
                              RankType root, CommunicatorType comm) {
  assert(dsgu->getRawDataSize() < INT_MAX && "Dsg is too large and can not be "
                                            "transferred in a single MPI Call (not "
                                            "supported yet) try a more refined"
                                            "decomposition");

  FG_ELEMENT* data = dsgu->getRawData();
  int dataSize  = static_cast<int>(dsgu->getRawDataSize());
  MPI_Datatype dataType = getMPIDatatype(abstraction::getabstractionDataType<FG_ELEMENT>());
  MPI_Request request;

  MPI_Ibcast(data, dataSize, dataType, root, comm, &request);
  return request;
}

/**
 * Bcast of the raw dsg data in the communicator comm.
 */
template <typename FG_ELEMENT>
static void bcastDsgData(DistributedSparseGridUniform<FG_ELEMENT> * dsgu,
                              RankType root, CommunicatorType comm) {
  assert(dsgu->getRawDataSize() < INT_MAX && "Dsg is too large and can not be "
                                            "transferred in a single MPI Call (not "
                                            "supported yet) try a more refined"
                                            "decomposition");

  FG_ELEMENT* data = dsgu->getRawData();
  int dataSize  = static_cast<int>(dsgu->getRawDataSize());
  MPI_Datatype dataType = getMPIDatatype(abstraction::getabstractionDataType<FG_ELEMENT>());

  MPI_Bcast(data, dataSize, dataType, root, comm);
}

/** Performs an in place allreduce on the dsgu data with all procs in
 * communicator comm.
 * This corresponds to a sparse grid reduce, cf. Heene
 */
template <typename FG_ELEMENT>
static void reduceDsgData(DistributedSparseGridUniform<FG_ELEMENT> * dsgu,
                               CommunicatorType comm) {
  assert(dsgu->getRawDataSize() < INT_MAX && "Dsg is too large and can not be "
                                            "transferred in a single MPI Call (not "
                                            "supported yet) try a more refined"
                                            "decomposition");

  // prepare for MPI call in globalReduceComm
  MPI_Datatype dtype = getMPIDatatype(
                        abstraction::getabstractionDataType<size_t>());
  const std::vector<size_t>& dsguData = dsgu->getRawData();

  // perform allreduce
  assert(dsguData.size() < static_cast<size_t>(std::numeric_limits<int>::max()));
  MPI_Allreduce(MPI_IN_PLACE, dsguData.data(), static_cast<int>(dsguData.size()), dtype, MPI_MAX, comm);
}

/**
* Sends all subspace data sizes to the receiver in communicator comm.
*/
template <typename FG_ELEMENT>
static void sendSubspaceDataSizes(DistributedSparseGridUniform<FG_ELEMENT> * dsgu,
                          const std::vector<LevelVector>& subspaces,
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


/** Performs a max allreduce in comm with subspace sizes of each dsg
 *
 * After calling, all workers which share the same spatial decomposition will
 * have the same subspace sizes and therefor in the end have equally sized dsgs.
 */
template <typename FG_ELEMENT>
static void reduceSubspaceSizes(DistributedSparseGridUniform<FG_ELEMENT> * dsgu,
                               CommunicatorType comm) {
  assert(dsgu->getNumSubspaces() > 0);

  // prepare for MPI call in globalReduceComm
  MPI_Datatype dtype = getMPIDatatype(
                        abstraction::getabstractionDataType<size_t>());

  const std::vector<size_t>& subspacesDataSizes = dsgu->getSubspaceDataSizes();
  std::vector<size_t> buf(subspacesDataSizes.size());

  // perform allreduce
  assert(buf.size() < static_cast<size_t>(std::numeric_limits<int>::max()));
  MPI_Allreduce(subspacesDataSizes.data(), buf.data(), static_cast<int>(buf.size()), dtype, MPI_MAX, comm);

  // set updated sizes in dsg
  for (size_t i = 0; i < subspacesDataSizes.size(); i++) {
    dsgu->setDataSize(i, buf[i]);
  }
}

} /* namespace combigrid */

#endif /* SRC_SGPP_COMBIGRID_SPARSEGRID_DISTRIBUTEDSPARSEGRID_HPP_ */