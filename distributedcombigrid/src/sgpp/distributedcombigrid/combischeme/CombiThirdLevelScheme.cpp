#include "sgpp/distributedcombigrid/combischeme/CombiThirdLevelScheme.hpp"

namespace combigrid {
/**
* Computes the distribution of a classical scheme to the systems of the
* third level combination.
*/
void CombiThirdLevelScheme::createThirdLevelScheme(std::vector<LevelVector>& levels,
                                                   std::vector<real>& coeffs,
                                                   std::vector<LevelVector>& commonSubspaces,
                                                   const std::vector<bool>& boundary,
                                                   unsigned int systemNumber,
                                                   unsigned int numSystems) {
  assert(!levels.empty() && !coeffs.empty());

  // decompose scheme
  std::vector<std::vector<LevelVector>> decomposedScheme;
  std::vector<std::vector<combigrid::real>> decomposedCoeffs;
  decomposeScheme(levels, coeffs, decomposedScheme, decomposedCoeffs, numSystems);

  // assign part to system
  levels  = decomposedScheme[systemNumber];
  coeffs = decomposedCoeffs[systemNumber];

  commonSubspaces = computeCommonSubspaces(decomposedScheme, boundary);
}


/**
* Computes an optimal disjunct decomposition of the given combination scheme,
* such that each part can be assigned to a system in the third level reduce and
* the amount shared grid points is minimal.
*
* We assume only 2 systems participating.
* For example purpose we just split the scheme in half and later assign each half
* to a system.
* TODO Implement for arbitrary number of systems
*/
void CombiThirdLevelScheme::decomposeScheme(std::vector<LevelVector>& fullScheme,
                                            std::vector<real> fullSchemeCoeffs,
                                            std::vector<std::vector<LevelVector>>& decomposedScheme,
                                            std::vector<std::vector<real>>& decomposedCoeffs,
                                            size_t numSystems) {
  // TODO
  auto mid = fullScheme.begin() + fullScheme.size()/2;
  auto midC = fullSchemeCoeffs.begin() + fullSchemeCoeffs.size()/2;
  std::vector<LevelVector> lowerHalf(fullScheme.begin(), mid);
  std::vector<real> lowerCoeffs(fullSchemeCoeffs.begin(), midC);
  std::vector<LevelVector> upperHalf(mid, fullScheme.end());
  std::vector<real> upperCoeffs(midC, fullSchemeCoeffs.end());

  assert( !lowerHalf.empty() && !upperHalf.empty() );
  //

  decomposedScheme.reserve(numSystems);
  decomposedCoeffs.reserve(numSystems);
  decomposedScheme.push_back(std::move(lowerHalf));
  decomposedScheme.push_back(std::move(upperHalf));
  decomposedCoeffs.push_back(std::move(lowerCoeffs));
  decomposedCoeffs.push_back(std::move(upperCoeffs));
}


/**
* Computes the common subspaces for a given decomposed scheme.
*/
std::vector<LevelVector> CombiThirdLevelScheme::computeCommonSubspaces(
                  const std::vector<std::vector<LevelVector>>& decomposedScheme,
                  const std::vector<bool>& boundary){
  // therefore we compute the component wise maximum level which is contained
  // in all sets
  long longMax = 2147483647;
  size_t dim = decomposedScheme[0][0].size();
  LevelVector maxLevel(dim, longMax);
  for (size_t d = 0; d < dim; d++) {
    for (const auto& levels : decomposedScheme) {
      long dimMax = 0;
      for (const auto& level : levels) {
        // find component wise maximum
        if (level[d] > dimMax)
          dimMax = level[d];
      }
      // minimum of all largest ss is contained in all sets
      if (dimMax < maxLevel[d])
        maxLevel[d] = dimMax;
    }
  }

  // by creating a dummy sg with this level we can extract the subspaces
  std::vector<LevelVector> commonSubspaces;
  SGrid<real> sg(dim, maxLevel, maxLevel, boundary);

  for (size_t ssID = 0; ssID < sg.getSize(); ++ssID) {
    const LevelVector& ss = sg.getLevelVector(ssID);
    commonSubspaces.push_back(ss);
    std::cout << toString(ss) << std::endl;
  }

  return commonSubspaces;
}

}
