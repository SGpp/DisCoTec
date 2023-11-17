#pragma once

#include <assert.h>
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <algorithm>
#include <numeric>
#include <vector>

#include "fullgrid/Tensor.hpp"
#include "utils/Types.hpp"

namespace combigrid {

class MPICartesianUtils {
 public:
  MPICartesianUtils() = default;

  MPICartesianUtils(CommunicatorType comm) : comm_(comm) {
    // check if communicator is cartesian
    int status;
    MPI_Topo_test(comm, &status);
    std::vector<int> cartdims;
    if (status == MPI_CART) {
      int ndims = 0;
      MPI_Cartdim_get(comm, &ndims);
      dim_ = static_cast<DimType>(ndims);
      cartdims.resize(dim_);
      periods_.resize(dim_);
      localCoords_.resize(dim_);
      MPI_Cart_get(comm, dim_, &cartdims[0], &periods_[0], &localCoords_[0]);
      // fill the partitionCoords_
      {
        partitionCoords_.resize(std::accumulate(cartdims.begin(), cartdims.end(), 1,
                                                std::multiplies<int>()));
        partitionCoordsTensor_ = Tensor(partitionCoords_.data(), std::move(cartdims));
        // fill partition coords vector, only once
        for (int i = 0; i < this->getCommunicatorSize(); ++i) {
          std::vector<int> tmp(dim_);
          MPI_Cart_coords(comm, i, static_cast<int>(dim_), &tmp[0]);
          auto sequentialIndex = partitionCoordsTensor_.sequentialIndex(tmp);
          partitionCoords_[sequentialIndex] = i;
        }
      }
    } else {
      comm_ = MPI_COMM_NULL;
      periods_.clear();
      localCoords_.clear();
      partitionCoords_.clear();
      partitionCoordsTensor_ = Tensor(partitionCoords_.data(), std::move(cartdims));
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

  explicit MPICartesianUtils(const MPICartesianUtils& other) = default;
  MPICartesianUtils& operator=(const MPICartesianUtils&) = default;
  explicit MPICartesianUtils(MPICartesianUtils&& other) = default;
  MPICartesianUtils& operator=(MPICartesianUtils&& other) = default;
  virtual ~MPICartesianUtils() = default;

  CommunicatorType getComm() const { return comm_; }

  /**
   * @brief Get the cartesian coordinates of a rank in the local cartesian communicator
   *
   * @param r         local rank
   * @param coords    out: coordinates in the cartesian grid of processes in local comm
   */
  inline const IndexVector& getPartitionCoordsOfRank(RankType r) const {
    assert(r >= 0 && r < getCommunicatorSize());
    assert(!partitionCoords_.empty());
    // find rank r in partitionCoords_
    auto findIt = std::find(partitionCoords_.begin(), partitionCoords_.end(), r);
    IndexType rIndex = static_cast<IndexType>(std::distance(partitionCoords_.begin(), findIt));
    return partitionCoordsTensor_.getVectorIndex(rIndex);
  }

  inline const IndexVector& getPartitionCoordsOfLocalRank() const {
    assert(!localCoords_.empty());
    return localCoords_;
  }

  inline bool isOnLowerBoundaryInDimension(DimType d) const { return localCoords_[d] == 0; }

  inline bool isOnUpperBoundaryInDimension(DimType d) const {
    return localCoords_[d] + 1 == this->getCartesianDimensions()[d];
  }

  inline RankType getRankFromPartitionCoords(const std::vector<int>& partitionCoordsInt) const {
    // check wheter the partition coords are valid
    assert(partitionCoordsInt.size() == dim_);

    for (DimType d = 0; d < dim_; ++d)
      assert(partitionCoordsInt[d] < this->getCartesianDimensions()[d]);

    assert(!partitionCoords_.empty());
    return partitionCoords_[partitionCoordsTensor_.sequentialIndex(partitionCoordsInt)];
  }

  RankType getNeighbor1dFromPartitionIndex(DimType dim, int idx1d) const {
    assert(idx1d >= 0);
    assert(idx1d < this->getCartesianDimensions()[dim]);

    static thread_local IndexVector neighborPartitionCoords;
    neighborPartitionCoords = this->localCoords_;
    neighborPartitionCoords[dim] = idx1d;
    RankType r = this->getRankFromPartitionCoords(neighborPartitionCoords);
    return r;
  }

  /**
   * @brief get a vector containing the ranks of all my cartesian neighboring
   *        ranks in dimension dim (not only the direct neighbors, all of them)
   */
  inline std::vector<RankType> getAllMyPoleNeighborRanks(DimType dim) const {
    auto ranks = std::vector<RankType>();
    ranks.reserve(this->getCartesianDimensions()[dim] - 1);
    const auto& myPartitionCoords = this->getPartitionCoordsOfLocalRank();
    for (int i = 0; i < myPartitionCoords[dim]; ++i) {
      auto neighborPartitionCoords = myPartitionCoords;
      neighborPartitionCoords[dim] = i;
      ranks.push_back(getRankFromPartitionCoords(neighborPartitionCoords));
    }
    for (int i = myPartitionCoords[dim] + 1; i < this->getCartesianDimensions()[dim]; ++i) {
      auto neighborPartitionCoords = myPartitionCoords;
      neighborPartitionCoords[dim] = i;
      ranks.push_back(getRankFromPartitionCoords(neighborPartitionCoords));
    }
    return ranks;
  }

  const std::vector<int>& getCartesianDimensions() const {
    return partitionCoordsTensor_.getExtentsVector();
  }

  inline int getCommunicatorSize() const {
    return static_cast<int>(partitionCoords_.size());
  }

  inline RankType getCommunicatorRank() const { return rank_; }

 private:
  CommunicatorType comm_ = MPI_COMM_NULL;

  RankType rank_ = MPI_UNDEFINED;

  DimType dim_;

  /** the periodicity of the cartesian communicator (periodic in each dimension) */
  std::vector<int> periods_;

  /** the cartesian coordinates of the calling process (where am I in each dimension) */
  std::vector<int> localCoords_;

  /** the coordinates of each rank on the cartesian communicator*/
  std::vector<int> partitionCoords_;

  /** multi-dim indexing for partitionCoords_ */
  Tensor<int> partitionCoordsTensor_;
};

}  // namespace combigrid