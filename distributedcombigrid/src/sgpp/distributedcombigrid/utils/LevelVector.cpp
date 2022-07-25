
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"

namespace combigrid {
std::string toString(combigrid::LevelVector const& l) {
  std::stringstream ss;
  for (size_t i = 0; i < l.size(); ++i) {
    if (i != 0) ss << ",";
    ss << l[i];
  }
  return ss.str();
}

void getDownSetRecursively(combigrid::LevelVector const& l, combigrid::LevelVector fixedDimensions,
                           std::vector<LevelVector>& downSet) {
  if (fixedDimensions.size() == l.size()) {
    downSet.push_back(fixedDimensions);
  } else {
    // for the whole range of values in dimension d
    size_t d = fixedDimensions.size();
    for (LevelType i = 1; i <= l[d]; ++i) {  // levels start at 1 here, in compliance with our
                                         // definition of DistributedSparseGrid
      auto moreDimensionsFixed = fixedDimensions;
      moreDimensionsFixed.push_back(i);
      getDownSetRecursively(l, moreDimensionsFixed, downSet);
    }
  }
}

std::vector<LevelVector> getDownSet(combigrid::LevelVector const& l) {
  std::vector<LevelVector> downSet;
  getDownSetRecursively(l, LevelVector(), downSet);
  return downSet;
}
}  // namespace combigrid
