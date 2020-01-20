#ifndef LOADMODEL_HPP_
#define LOADMODEL_HPP_

#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/utils/Config.hpp"

namespace combigrid {

/**
 * The LoadModel is an interface to be implemented to provide the expected load
 * of specific tasks. It is to be used for scheduling the Tasks.
 */
class LoadModel {
 public:
  virtual ~LoadModel() = default;

  /**
   * Calculates the expected load of a task.
   *
   * The calculated load is in no specific unit. However, the 
   * "greater than"-relation for two tasks has to be correct.
   *
   * @param l The level vector corresponding to a task.
   * @returns Expected load of the given task.
   */
  virtual real eval(const LevelVector& l) = 0;

};

} /* namespace combigrid */
#endif /* LOADMODEL_HPP_ */
