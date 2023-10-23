#include "sparsegrid/SGrid.hpp"

namespace combigrid {

class CombiThirdLevelScheme {
 public:
  /** Computes the distribution of a classical scheme to the systems of the
   * third level combination.
   */
  static void createThirdLevelScheme(const std::vector<LevelVector>& levels,
                                     const std::vector<real>& coeffs, unsigned int systemNumber,
                                     unsigned int numSystems, std::vector<LevelVector>& newLevels,
                                     std::vector<real>& newCoeffs,
                                     const std::vector<real>& fractionsOfScheme = {0.5, 0.5});

 private:
  /** Computes a (non-optimal) disjunct decomposition of the given combination scheme.
   * Each part can be assigned to a system in the third level reduce.
   */
  static void decomposeScheme(const std::vector<LevelVector>& fullScheme,
                              const std::vector<real> fullSchemeCoeffs,
                              std::vector<std::vector<LevelVector>>& decomposedScheme,
                              std::vector<std::vector<real>>& decomposedCoeffs,
                              size_t numSystems = 2,
                              const std::vector<real>& fractionsOfScheme = {0.5, 0.5});
};

}  // namespace combigrid
