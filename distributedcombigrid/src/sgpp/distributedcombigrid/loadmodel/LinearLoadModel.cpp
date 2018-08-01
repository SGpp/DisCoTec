#include "sgpp/distributedcombigrid/loadmodel/LinearLoadModel.hpp"

namespace combigrid {

real LinearLoadModel::eval(const LevelVector& l) const {
  real ret(1.0);

  for (size_t i = 0; i < l.size(); ++i) {
    ret *= std::pow(real(2.0), static_cast<real>(l[i]));
  }

  return ret;
}
}
