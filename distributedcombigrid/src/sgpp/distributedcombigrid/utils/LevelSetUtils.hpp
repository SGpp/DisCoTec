#pragma once

#include <vector>

#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/utils/PowerOfTwo.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"

namespace combigrid {

inline IndexType getNumDofNodal(const LevelType& l_i, const BoundaryType& boundary) {
  assert(l_i > -1);
  if (l_i > 0) {
    return powerOfTwoByBitshift(l_i) + boundary - 1;
  } else {
    if (boundary == 0) {
      throw std::runtime_error("without boundary, there are no points on level 0");
    }
    return boundary;
  }
}

inline IndexType getNumDofNodal(const LevelVector& l, const std::vector<BoundaryType>& boundary) {
  assert(l.size() == boundary.size());
  IndexType numDofPerGrid = 1;
  auto b = boundary.cbegin();
  for (const auto& l_i : l) {
    numDofPerGrid *= getNumDofNodal(l_i, *b);
    ++b;
  }
  return numDofPerGrid;
}

inline IndexType getNumDofNodal(const std::vector<LevelVector>& levelVectors,
                                const std::vector<BoundaryType>& boundary) {
  IndexType numDof = 0;
  for (const auto& l : levelVectors) {
    numDof += getNumDofNodal(l, boundary);
  }
  return numDof;
}

inline IndexType getNumDofHierarchical(const LevelType& l_i, const BoundaryType& boundary) {
  assert(l_i > 0);
  if (l_i == 1) {
    return 1 + boundary;
  } else {
    return powerOfTwoByBitshift(static_cast<LevelType>(l_i - 1));
  }
}

inline IndexType getNumDofHierarchical(const LevelVector& l, const std::vector<BoundaryType>& boundary) {
  assert(l.size() == boundary.size());
  auto b = boundary.cbegin();
  IndexType numDofPerSubspace = 1;
  for (const auto& l_i : l) {
    numDofPerSubspace *= getNumDofHierarchical(l_i, *b);
    ++b;
  }
  return numDofPerSubspace;
}

inline IndexType getNumDofHierarchical(const std::vector<LevelVector>& levelVectors,
                                       const std::vector<BoundaryType>& boundary) {
  IndexType numDof = 0;
  for (const auto& l : levelVectors) {
    numDof += getNumDofHierarchical(l, boundary);
  }
  return numDof;
}

template <typename T>
inline bool isPowerOfTwo(const T& x) {
  return (x != 0) && ((x & (x - 1)) == 0);
}

// get downward closed set of a single LevelVector
std::vector<LevelVector> getDownSet(combigrid::LevelVector const& l);

struct AllKOutOfDDimensions {
  /**
   * @brief Get all combinations of k out of d dimensions (from 0 to d-1)
   */
  static const std::vector<std::vector<DimType>>& get(DimType k, DimType d);

  static std::map<std::pair<DimType, DimType>, std::vector<std::vector<DimType>>> cache_;
};

/**
 * @brief recursively generate a downward-closed set of hierarchical level vectors
 *
 * @param dim : the currently recursively iterated dimension
 * @param n : the "regular level", here the minimum of the difference of lmax and lmin
 * @param l : the currently populated level vector (entries filled only from dim to d-1)
 * @param lmax
 * @param lmin
 * @param created : output vector of subspaces
 *
 *  start recursion by setting dim=d=dimensionality of the vector space.
    For instance, for correct subspace restriction in d > 2, we need 3 criteria:
    the diagonal hyperplane that restricts to the simplex,
    the maximum of lmax in every dimension, and
    the mixed dimension sum restrictions
 */
void createTruncatedHierarchicalLevelsRec(DimType dim, size_t n, LevelVector& l,
                                          const LevelVector& lmax, const LevelVector& lmin,
                                          std::vector<LevelVector>& created);

void createTruncatedHierarchicalLevels(const LevelVector& lmax, const LevelVector& lmin,
                                       std::vector<LevelVector>& created);

}  // namespace combigrid
