#ifndef LEVELVECTOR_HPP_
#define LEVELVECTOR_HPP_

// #include <sstream>
// #include <string>
#include "sgpp/distributedcombigrid/utils/IndexVector.hpp"

namespace combigrid {

typedef IndexVector LevelVector;
std::string toString(combigrid::LevelVector const& l);

// get downward closed set of a single LevelVector
std::vector<LevelVector> getDownSet(combigrid::LevelVector const& l);

/**
 * @brief recursively generate a downward-closed set of hierarchical level vectors
 *
 * @param dim : the currently recursively iterated dimension
 * @param n : the "regular level", here the minimum of the difference of lmax and lmin
 * @param d : the dimensionality
 * @param l : the currently populated level vector (entries filled only from dim to d-1)
 * @param lmax
 * @param lmin
 *
 *  start recursion by setting dim=d=dimensionality of the vector space
    for correct subspace restriction in d > 2, we need 3 criteria:
    the diagonal hyperplane that restricts to the simplex,
    the maximum of lmax in every dimension, and
    the mixed dimension sum restrictions
 */
void createTruncatedHierarchicalLevelsRec(size_t dim, size_t n, LevelVector& l,
                                          const LevelVector& lmax, const LevelVector& lmin,
                                          std::vector<LevelVector>& created);

void createTruncatedHierarchicalLevels(const LevelVector& lmax, const LevelVector& lmin,
                                       std::vector<LevelVector>& created);

}  // namespace combigrid

#endif /* LEVELVECTOR_HPP_ */
