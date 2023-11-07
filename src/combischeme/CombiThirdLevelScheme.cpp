#include "combischeme/CombiThirdLevelScheme.hpp"

#include <numeric>
#include <vector>

namespace combigrid {
void CombiThirdLevelScheme::createThirdLevelScheme(
    const std::vector<LevelVector>& levels, const std::vector<real>& coeffs,
    unsigned int systemNumber, unsigned int numSystems, std::vector<LevelVector>& newLevels,
    std::vector<real>& newCoeffs, const std::vector<real>& fractionsOfScheme) {
  assert(!levels.empty() && !coeffs.empty());
  assert(levels.size() == coeffs.size());

  if (numSystems < 2) {
    newLevels = levels;
    newCoeffs = coeffs;
    return;
  }

  // decompose scheme
  std::vector<std::vector<LevelVector>> decomposedScheme;
  std::vector<std::vector<combigrid::real>> decomposedCoeffs;
  decomposeScheme(levels, coeffs, decomposedScheme, decomposedCoeffs, numSystems,
                  fractionsOfScheme);

  // assign part to system
  newLevels = decomposedScheme[systemNumber];
  newCoeffs = decomposedCoeffs[systemNumber];
}

void CombiThirdLevelScheme::decomposeScheme(const std::vector<LevelVector>& fullScheme,
                                            const std::vector<real>& fullSchemeCoeffs,
                                            std::vector<std::vector<LevelVector>>& decomposedScheme,
                                            std::vector<std::vector<real>>& decomposedCoeffs,
                                            size_t numSystems,
                                            const std::vector<real>& fractionsOfScheme) {
  [[maybe_unused]] auto fracSum =
      std::accumulate(fractionsOfScheme.begin(), fractionsOfScheme.end(), 0., std::plus<real>());
  assert(std::abs(fracSum - 1.) < 1e-3);
  assert(fractionsOfScheme.size() == numSystems);
  decomposedScheme.reserve(numSystems);
  decomposedCoeffs.reserve(numSystems);

  real scannedFrac = 0.;
  auto beginNextL = fullScheme.begin();
  auto beginNextC = fullSchemeCoeffs.begin();
  for (auto frac : fractionsOfScheme) {
    scannedFrac += frac;
    auto currentSystemUpToIndex = std::round(static_cast<real>(fullScheme.size()) * scannedFrac);
    auto endIntervalL = fullScheme.begin() + currentSystemUpToIndex;
    auto endIntervalC = fullSchemeCoeffs.begin() + currentSystemUpToIndex;
    decomposedScheme.emplace_back(beginNextL, endIntervalL);
    decomposedCoeffs.emplace_back(beginNextC, endIntervalC);
    beginNextL = endIntervalL;
    beginNextC = endIntervalC;
  }

  // assert that all are assigned
  assert(std::accumulate(decomposedScheme.begin(), decomposedScheme.end(), static_cast<size_t>(0),
                         [](size_t a, const std::vector<LevelVector>& l) {
                           return a + l.size();
                         }) == fullScheme.size());
  assert(std::accumulate(decomposedCoeffs.begin(), decomposedCoeffs.end(), 0,
                         [](size_t a, const std::vector<real>& c) { return a + c.size(); }) ==
         fullSchemeCoeffs.size());
}
}  // namespace combigrid
