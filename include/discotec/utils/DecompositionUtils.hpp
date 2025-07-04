#pragma once
#include <vector>

#include "Types.hpp"
namespace combigrid {

/* a regular (equidistant) domain decompositioning for an even number of processes
 * leads to grid points on the (geometrical) process boundaries.
 * with the forwardDecomposition flag it can be decided if the grid points on
 * the process boundaries belong to the process on the right-hand side (true)
 * of the process boundary, or to the one on the left-hand side (false).
 */
static inline std::vector<IndexVector> getDefaultDecomposition(
    const IndexVector& globalNumPointsPerDimension,
    const std::vector<int>& cartesianProcsPerDimension, bool forwardDecomposition) {
  auto dim = static_cast<DimType>(globalNumPointsPerDimension.size());
  assert(cartesianProcsPerDimension.size() == dim);

  // create decomposition vectors
  std::vector<IndexVector> decomposition(dim);
  for (DimType i = 0; i < dim; ++i) {
    if (cartesianProcsPerDimension[i] > globalNumPointsPerDimension[i]) {
      throw std::invalid_argument(
          "change to coarser parallelization! currently not all processes can have points.");
    }
    IndexVector& llbnd = decomposition[i];
    llbnd.resize(cartesianProcsPerDimension[i]);

    for (int j = 0; j < cartesianProcsPerDimension[i]; ++j) {
      double tmp = static_cast<double>(globalNumPointsPerDimension[i]) * static_cast<double>(j) /
                   static_cast<double>(cartesianProcsPerDimension[i]);

      if (forwardDecomposition)
        llbnd[j] = static_cast<IndexType>(std::ceil(tmp));
      else
        llbnd[j] = static_cast<IndexType>(std::floor(tmp));
    }
  }
  return decomposition;
}

inline static std::vector<IndexVector> getStandardDecomposition(LevelVector lref,
                                                                std::vector<int>& procsRef) {
  assert(lref.size() == procsRef.size());
  std::vector<IndexVector> decomposition;
  for (DimType d = 0; d < static_cast<DimType>(lref.size()); ++d) {
    IndexVector di;
    if (procsRef[d] == 1) {
      di = {0};
    } else if (procsRef[d] == 2) {
      di = {0, powerOfTwo[lref[d]] / procsRef[d] + 1};
    } else if (procsRef[d] == 3) {
      di = {0, powerOfTwo[lref[d]] / procsRef[d] + 1, 2 * powerOfTwo[lref[d]] / procsRef[d] + 1};
    } else if (procsRef[d] == 4) {
      di = {0, powerOfTwo[lref[d]] / procsRef[d] + 1, 2 * powerOfTwo[lref[d]] / procsRef[d] + 1,
            3 * powerOfTwo[lref[d]] / procsRef[d] + 1};
    } else {
      throw std::runtime_error("please implement a test decomposition matching procs and lref");
    }
    decomposition.push_back(di);
  }
  return decomposition;
}

static inline std::vector<IndexVector> downsampleDecomposition(
    const std::vector<IndexVector>& decomposition, const LevelVector& referenceLevel,
    const LevelVector& newLevel, const std::vector<BoundaryType>& boundary) {
  auto newDecomposition = decomposition;
  if (decomposition.size() > 0) {
    for (DimType d = 0; d < static_cast<DimType>(referenceLevel.size()); ++d) {
      // for now, assume that we never want to interpolate on a level finer than referenceLevel
      assert(referenceLevel[d] >= newLevel[d]);
      auto levelDiff = referenceLevel[d] - newLevel[d];
      auto stepFactor = oneOverPowOfTwo[levelDiff];
      if (boundary[d] > 0) {
        // all levels contain the boundary points -> point 0 is the same
        for (auto& dec : newDecomposition[d]) {
          dec = static_cast<IndexType>(std::ceil(static_cast<double>(dec) * stepFactor));
        }
      } else {
        // all levels do not contain the boundary points -> mid point is the same
        auto leftProtrusion = powerOfTwo[levelDiff] - 1;
        for (auto& dec : newDecomposition[d]) {
          // same as before, but subtract the "left" protrusion on the finer level
          dec = static_cast<IndexType>(
              std::max(0., std::ceil(static_cast<double>(dec - leftProtrusion) * stepFactor)));
        }
      }
    }
  }
  return newDecomposition;
}

}  // namespace combigrid