#pragma once

#include <assert.h>
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <algorithm>
#include <numeric>
#include <vector>

#include "utils/Types.hpp"

namespace combigrid {

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
      dim_ = static_cast<DimType>(ndims);
      cartdims_.resize(dim_);
      periods_.resize(dim_);
      localCoords_.resize(dim_);
      MPI_Cart_get(comm, dim_, &cartdims_[0], &periods_[0], &localCoords_[0]);
      // fill the partitionCoords_
      {
        partitionCoords_.resize(this->getCommunicatorSize());
        // fill partition coords vector, only once
        for (int i = 0; i < this->getCommunicatorSize(); ++i) {
          std::vector<int> tmp(dim_);
          MPI_Cart_coords(comm, i, static_cast<int>(dim_), &tmp[0]);

          // important: reverse ordering of partition coords!
          if (reverseOrderingDFGPartitions) {
            partitionCoords_[i].assign(tmp.rbegin(), tmp.rend());
          } else {
            partitionCoords_[i].assign(tmp.begin(), tmp.end());
          }
        }
      }
    } else {
      cartdims_.clear();
      periods_.clear();
      localCoords_.clear();
      partitionCoords_.clear();
    }
  }

  explicit MPICartesianUtils(const MPICartesianUtils& other) = default;
  MPICartesianUtils& operator=(const MPICartesianUtils&) = default;
  explicit MPICartesianUtils(MPICartesianUtils&& other) = default;
  MPICartesianUtils& operator=(MPICartesianUtils&& other) = default;
  virtual ~MPICartesianUtils() = default;

  CommunicatorType getComm() const {
    return comm_;
  }

  /**
   * @brief Get the cartesian coordinates of a rank in the local cartesian communicator
   *
   * @param r         local rank
   * @param coords    out: coordinates in the cartesian grid of processes in local comm
   */
  inline void getPartitionCoordsOfRank(RankType r, std::vector<int>& coords) const {
    assert(r >= 0 && r < getCommunicatorSize());
    assert(!partitionCoords_.empty());
    coords.assign(partitionCoords_[r].begin(), partitionCoords_[r].end());
  }

  inline void getPartitionCoordsOfLocalRank(std::vector<int>& coords) const {
    assert(!localCoords_.empty());
    coords = localCoords_;
  }

  inline bool isOnLowerBoundaryInDimension(DimType d) const { return localCoords_[d] == 0; }

  inline bool isOnUpperBoundaryInDimension(DimType d) const {
    return localCoords_[d] + 1 == cartdims_[d];
  }

  inline RankType getRankFromPartitionCoords(const std::vector<int>& partitionCoordsInt) const {
    // check wheter the partition coords are valid
    assert(partitionCoordsInt.size() == dim_);

    for (DimType d = 0; d < dim_; ++d) assert(partitionCoordsInt[d] < cartdims_[d]);

    assert(!partitionCoords_.empty());
    auto dim = dim_;
    auto it = std::find_if(partitionCoords_.begin(), partitionCoords_.end(),
                           [&partitionCoordsInt, &dim](const std::vector<int>& pcoord) {
                             return (pcoord == partitionCoordsInt);
                           });
    // check if found
    assert(it != partitionCoords_.end());

    RankType rank = static_cast<RankType>(it - partitionCoords_.begin());

    return rank;
  }

  RankType getNeighbor1dFromPartitionIndex(DimType dim, int idx1d) const {
    assert(idx1d >= 0);
    assert(idx1d < this->getCartesianDimensions()[dim]);

    auto neighborPartitionCoords = this->localCoords_;
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
    ranks.reserve(cartdims_[dim] - 1);
    auto myPartitionCoords = std::vector<int>();
    this->getPartitionCoordsOfLocalRank(myPartitionCoords);
    for (int i = 0; i < myPartitionCoords[dim]; ++i) {
      auto neighborPartitionCoords = myPartitionCoords;
      neighborPartitionCoords[dim] = i;
      ranks.push_back(getRankFromPartitionCoords(neighborPartitionCoords));
    }
    for (int i = myPartitionCoords[dim] + 1; i < cartdims_[dim]; ++i) {
      auto neighborPartitionCoords = myPartitionCoords;
      neighborPartitionCoords[dim] = i;
      ranks.push_back(getRankFromPartitionCoords(neighborPartitionCoords));
    }
    return ranks;
  }

  const std::vector<int>& getCartesianDimensions() const {
    return cartdims_;
  }

  inline int getCommunicatorSize() const {
    if (cartdims_.empty()) return 0;
    return std::accumulate(cartdims_.begin(), cartdims_.end(), 1, std::multiplies<int>());
  }

 private:
  CommunicatorType comm_ = MPI_COMM_NULL;

  DimType dim_;

  /** the extents of the cartesian communicator (how many ranks in each dimension) */
  std::vector<int> cartdims_;

  /** the periodicity of the cartesian communicator (periodic in each dimension) */
  std::vector<int> periods_;

  /** the cartesian coordinates of the calling process (where am I in each dimension) */
  std::vector<int> localCoords_;

  /** the coordinates of each rank on the cartesian communicator*/
  std::vector<std::vector<int>> partitionCoords_;
};

}  // namespace combigrid