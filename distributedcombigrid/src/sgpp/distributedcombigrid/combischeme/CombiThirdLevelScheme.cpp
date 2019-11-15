#include "sgpp/distributedcombigrid/combischeme/CombiThirdLevelScheme.hpp"

namespace combigrid {
/**Computes the distribution of a classical scheme to the systems of the
 * third level combination.
 */
void CombiThirdLevelScheme::createThirdLevelScheme(const std::vector<LevelVector>& levels,
                                                   const std::vector<real>& coeffs,
                                                   const std::vector<bool>& boundary,
                                                   unsigned int systemNumber,
                                                   unsigned int numSystems,
                                                   std::vector<LevelVector>& newLevels,
                                                   std::vector<real>& newCoeffs) {
  assert(!levels.empty() && !coeffs.empty());

  if (numSystems < 2) {
    newLevels = levels;
    newCoeffs = coeffs;
    return;
  }

  assert(numSystems <= 2 && "Implemented for maximal 2 systems");

  // decompose scheme
  std::vector<std::vector<LevelVector>> decomposedScheme;
  std::vector<std::vector<combigrid::real>> decomposedCoeffs;
  decomposeScheme(levels, coeffs, decomposedScheme, decomposedCoeffs, numSystems);

  // assign part to system
  newLevels  = decomposedScheme[systemNumber];
  newCoeffs = decomposedCoeffs[systemNumber];

}


/**Computes an optimal disjunct decomposition of the given combination scheme.
 * Each part can be assigned to a system in the third level reduce and
 * the amount of shared grid points is minimal.
 *
 * We assume only 2 systems participating.
 * For example purpose we just split the scheme in half and later assign each half
 * to a system.
 * TODO Implement for arbitrary number of systems
 */
void CombiThirdLevelScheme::decomposeScheme(const std::vector<LevelVector>& fullScheme,
                                            const std::vector<real> fullSchemeCoeffs,
                                            std::vector<std::vector<LevelVector>>& decomposedScheme,
                                            std::vector<std::vector<real>>& decomposedCoeffs,
                                            size_t numSystems) {
  auto mid = fullScheme.begin() + fullScheme.size()/2;
  auto midC = fullSchemeCoeffs.begin() + fullSchemeCoeffs.size()/2;
  std::vector<LevelVector> lowerHalf(fullScheme.begin(), mid);
  std::vector<real> lowerCoeffs(fullSchemeCoeffs.begin(), midC);
  std::vector<LevelVector> upperHalf(mid, fullScheme.end());
  std::vector<real> upperCoeffs(midC, fullSchemeCoeffs.end());

  assert( !lowerHalf.empty() && !upperHalf.empty() );
  assert(lowerHalf.size() + upperHalf.size() == fullScheme.size());
  assert(lowerCoeffs.size() + upperCoeffs.size() == fullSchemeCoeffs.size());

  decomposedScheme.reserve(numSystems);
  decomposedCoeffs.reserve(numSystems);
  decomposedScheme.push_back(std::move(lowerHalf));
  decomposedScheme.push_back(std::move(upperHalf));
  decomposedCoeffs.push_back(std::move(lowerCoeffs));
  decomposedCoeffs.push_back(std::move(upperCoeffs));
}
}
