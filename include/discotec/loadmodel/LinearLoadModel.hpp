#ifndef LINEARLOADMODEL_HPP_
#define LINEARLOADMODEL_HPP_

#include "LoadModel.hpp"
#include "../utils/LevelVector.hpp"
#include "../utils/Types.hpp"

namespace combigrid {

class LinearLoadModel : public LoadModel {
 public:
  LinearLoadModel() = default;

  ~LinearLoadModel() = default;

  virtual real eval(const LevelVector& l);
};

} /* namespace combigrid */
#endif /* LINEARLOADMODEL_HPP_ */
