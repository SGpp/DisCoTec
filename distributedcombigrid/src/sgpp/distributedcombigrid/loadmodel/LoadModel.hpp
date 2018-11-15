/*
 * LoadModel.hpp
 *
 *  Created on: Oct 9, 2013
 *      Author: heenemo
 */

#ifndef LOADMODEL_HPP_
#define LOADMODEL_HPP_

#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"

namespace combigrid {

class LoadModel {
 public:
  virtual ~LoadModel() = default;
  virtual real eval(const LevelVector& l) = 0;

};

} /* namespace combigrid */
#endif /* LOADMODEL_HPP_ */
