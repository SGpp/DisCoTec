#ifndef LINEARLOADMODEL_HPP_
#define LINEARLOADMODEL_HPP_

#include "discotec/loadmodel/LoadModel.hpp"
#include "discotec/utils/LevelVector.hpp"
#include "discotec/utils/Types.hpp"

namespace combigrid {

class LinearLoadModel : public LoadModel {
 public:
  LinearLoadModel() = default;

  ~LinearLoadModel() = default;

  virtual real eval(const LevelVector& l);
};

} /* namespace combigrid */
#endif /* LINEARLOADMODEL_HPP_ */
