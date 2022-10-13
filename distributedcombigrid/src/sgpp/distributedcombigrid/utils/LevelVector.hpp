#ifndef LEVELVECTOR_HPP_
#define LEVELVECTOR_HPP_

#include <numeric>
#include <string>
#include <vector>

#include "sgpp/distributedcombigrid/utils/IndexVector.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"

namespace combigrid {

typedef std::vector<LevelType> LevelVector;
std::string toString(combigrid::LevelVector const& l);

inline LevelType levelSum(const LevelVector& l) {
  return std::accumulate(l.begin(), l.end(), static_cast<LevelType>(0));
}

}  // namespace combigrid

#endif /* LEVELVECTOR_HPP_ */
