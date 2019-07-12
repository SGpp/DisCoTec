#include "sgpp/distributedcombigrid/sparsegrid/SGrid.hpp"
#include <algorithm>

namespace combigrid {

class CombiThirdLevelScheme {

  public:
    /*
     * Creates system specific scheme for third level combination based on the
     * given full scheme (levels, coeffs). Additionally returns the list of
     * subspaces which all participating systems have in common.
     */
    static void createThirdLevelScheme(std::vector<LevelVector>& levels,
                                       std::vector<real>& coeffs,
                                       std::vector<LevelVector>& commonSubspaces,
                                       const std::vector<bool>& boundary,
                                       unsigned int systemNumber,
                                       unsigned int numSystems);

  private:
    static void decomposeScheme(std::vector<LevelVector>& fullScheme,
                                std::vector<real> fullSchemeCoeffs,
                                std::vector<std::vector<LevelVector>>& decomposedScheme,
                                std::vector<std::vector<real>>& decomposedCoeffs,
                                size_t numSystems = 2);

    static std::vector<LevelVector> computeCommonSubspaces(
                    const std::vector<std::vector<LevelVector>>& splittedScheme,
                    const std::vector<bool>& boundary);
};

}
