#pragma once

#include <assert.h>
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <algorithm>
#include <numeric>
#include <vector>

#include "discotec/fullgrid/Tensor.hpp"
#include "discotec/utils/Types.hpp"

namespace combigrid {

template <DimType DIM>
class MPICartesianUtils {
 public:
  MPICartesianUtils() = default;

  MPICartesianUtils(CommunicatorType comm) : comm_(comm) {
    // check if communicator is cartesian
    int status;
    MPI_Topo_test(comm, &status);
    if (status == MPI_CART) {
      int ndims = 0;
      MPI_Cartdim_get(comm, &ndims);
      assert(ndims == static_cast<int>(DIM));
      MPI_Cart_get(comm, DIM, cartdims_.data(), periods_.data(), localCoords_.data());
      // fill the partitionCoords_
      {
        IndexArray<DIM> extents;
        for (DimType j = 0; j < DIM; ++j) {
          extents[j] = cartdims_[j];
        }
        partitionCoordsIndexer_ = TensorIndexer<DIM>(extents);
        partitionCoords_.resize(partitionCoordsIndexer_.size());
        // fill partition coords vector, only once
        for (int i = 0; i < this->getCommunicatorSize(); ++i) {
          std::array<int, DIM> tmp{};
          MPI_Cart_coords(comm, i, static_cast<int>(DIM), tmp.data());
          IndexArray<DIM> tmpIndex;
          for (DimType j = 0; j < DIM; ++j) {
            tmpIndex[j] = tmp[j];
          }
          auto sequentialIndex = partitionCoordsIndexer_.sequentialIndex(tmpIndex);
          partitionCoords_[sequentialIndex] = i;
        }
      }
    } else {
      comm_ = MPI_COMM_NULL;
      partitionCoords_.clear();
      throw std::runtime_error("MPICartesianUtils: communicator is not cartesian");
    }
    MPI_Comm_rank(comm, &rank_);

    int size = 0;
    MPI_Comm_size(comm, &size);
    if (size != this->getCommunicatorSize()) {
      throw std::runtime_error(
          "MPICartesianUtils: communicator size does not match cartesian "
          "dimensions");
    }
  }

  explicit MPICartesianUtils(const MPICartesianUtils& other) = delete;
  MPICartesianUtils& operator=(const MPICartesianUtils&) = delete;
  explicit MPICartesianUtils(MPICartesianUtils&& other) = default;
  MPICartesianUtils& operator=(MPICartesianUtils&& other) = default;
  virtual ~MPICartesianUtils() = default;

  CommunicatorType getComm() const { return comm_; }

  /**
   * @brief Get the cartesian coordinates of a rank in the local cartesian communicator
   *
   * @param r         local rank
   * @return          coordinates in the cartesian grid of processes in local comm
   */
  inline IndexArray<DIM> getPartitionCoordsOfRank(RankType r) const {
    assert(r >= 0 && r < getCommunicatorSize());
    assert(!partitionCoords_.empty());
    // find rank r in partitionCoords_
    auto findIt = std::find(partitionCoords_.begin(), partitionCoords_.end(), r);
    IndexType rIndex = static_cast<IndexType>(std::distance(partitionCoords_.begin(), findIt));
    return partitionCoordsIndexer_.getArrayIndex(rIndex);
  }

  inline const std::array<int, DIM>& getPartitionCoordsOfLocalRank() const { return localCoords_; }

  inline bool isOnLowerBoundaryInDimension(DimType d) const { return localCoords_[d] == 0; }

  inline bool isOnUpperBoundaryInDimension(DimType d) const {
    return localCoords_[d] + 1 == cartdims_[d];
  }

  inline RankType getRankFromPartitionCoords(const std::array<int, DIM>& partitionCoordsInt) const {
    for (DimType d = 0; d < DIM; ++d) assert(partitionCoordsInt[d] < cartdims_[d]);

    assert(!partitionCoords_.empty());
    IndexArray<DIM> idx;
    for (DimType d = 0; d < DIM; ++d) {
      idx[d] = partitionCoordsInt[d];
    }
    return partitionCoords_[partitionCoordsIndexer_.sequentialIndex(idx)];
  }

  // backward-compatible overload accepting vector
  inline RankType getRankFromPartitionCoords(const std::vector<int>& partitionCoordsInt) const {
    assert(partitionCoordsInt.size() == DIM);
    std::array<int, DIM> arr;
    std::copy_n(partitionCoordsInt.begin(), DIM, arr.begin());
    return getRankFromPartitionCoords(arr);
  }

  RankType getNeighbor1dFromPartitionIndex(DimType dim, int idx1d) const {
    assert(idx1d >= 0);
    assert(idx1d < cartdims_[dim]);

    std::array<int, DIM> neighborPartitionCoords = localCoords_;
    neighborPartitionCoords[dim] = idx1d;
    return this->getRankFromPartitionCoords(neighborPartitionCoords);
  }

  /**
   * @brief get a vector containing the ranks of all my cartesian neighboring
   *        ranks in dimension dim (not only the direct neighbors, all of them)
   */
  inline std::vector<RankType> getAllMyPoleNeighborRanks(DimType dim) const {
    auto ranks = std::vector<RankType>();
    ranks.reserve(cartdims_[dim] - 1);
    const auto& myPartitionCoords = this->getPartitionCoordsOfLocalRank();
    for (int i = 0; i < myPartitionCoords[dim]; ++i) {
      std::array<int, DIM> neighborPartitionCoords = myPartitionCoords;
      neighborPartitionCoords[dim] = i;
      ranks.push_back(getRankFromPartitionCoords(neighborPartitionCoords));
    }
    for (int i = myPartitionCoords[dim] + 1; i < cartdims_[dim]; ++i) {
      std::array<int, DIM> neighborPartitionCoords = myPartitionCoords;
      neighborPartitionCoords[dim] = i;
      ranks.push_back(getRankFromPartitionCoords(neighborPartitionCoords));
    }
    return ranks;
  }

  const std::array<int, DIM>& getCartesianDimensions() const { return cartdims_; }

  // backward-compatible: return as vector
  std::vector<int> getCartesianDimensionsVector() const {
    return std::vector<int>(cartdims_.begin(), cartdims_.end());
  }

  inline int getCommunicatorSize() const { return static_cast<int>(partitionCoords_.size()); }

  inline RankType getCommunicatorRank() const { return rank_; }

 private:
  CommunicatorType comm_ = MPI_COMM_NULL;

  RankType rank_ = MPI_UNDEFINED;

  /** the cartesian dimensions (number of procs per dimension) */
  std::array<int, DIM> cartdims_{};

  /** the periodicity of the cartesian communicator (periodic in each dimension) */
  std::array<int, DIM> periods_{};

  /** the cartesian coordinates of the calling process (where am I in each dimension) */
  std::array<int, DIM> localCoords_{};

  /** the coordinates of each rank on the cartesian communicator*/
  std::vector<int> partitionCoords_;

  /** multi-dim indexing for partitionCoords_ */
  TensorIndexer<DIM> partitionCoordsIndexer_;
};

}  // namespace combigrid
