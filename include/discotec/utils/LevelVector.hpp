#ifndef LEVELVECTOR_HPP_
#define LEVELVECTOR_HPP_

#include <numeric>
#include <string>
#include <vector>

#include "IndexVector.hpp"
#include "Types.hpp"

namespace combigrid {

typedef std::vector<LevelType> LevelVector;

template <DimType NumDimensions>
using LevelArray = std::array<LevelType, NumDimensions>;

std::string toString(combigrid::LevelVector const& l);

inline LevelType levelSum(const LevelVector& l) {
  return std::accumulate(l.begin(), l.end(), static_cast<LevelType>(0));
}

}  // namespace combigrid

#endif /* LEVELVECTOR_HPP_ */
