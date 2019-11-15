#include "sgpp/distributedcombigrid/sparsegrid/SGrid.hpp"
#include <algorithm>

namespace combigrid {

class CombiThirdLevelScheme {

  public:
    /**Creates system specific scheme for third level combination based on the
     * given full scheme (levels, coeffs). Additionally returns the list of
     * subspaces which all participating systems have in common.
     *
     * I packed that into a separate class because scheme generation should be
     * supported for other schemes than minmax e.g. schemes read directly from
     * file.
     */
    static void createThirdLevelScheme(const std::vector<LevelVector>& levels,
                                       const std::vector<real>& coeffs,
                                       const std::vector<bool>& boundary,
                                       unsigned int systemNumber,
                                       unsigned int numSystems,
                                       std::vector<LevelVector>& newLevels,
                                       std::vector<real>& newCoeffs);

  private:
    static void decomposeScheme(const std::vector<LevelVector>& fullScheme,
                                const std::vector<real> fullSchemeCoeffs,
                                std::vector<std::vector<LevelVector>>& decomposedScheme,
                                std::vector<std::vector<real>>& decomposedCoeffs,
                                size_t numSystems = 2);
};

}
