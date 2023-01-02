#ifndef LOADMODEL_HPP_
#define LOADMODEL_HPP_

#include "utils/LevelVector.hpp"
#include "utils/Config.hpp"

namespace combigrid {

/**
 * The LoadModel is an interface to be implemented to provide the expected load
 * of specific tasks. It is to be used for scheduling the Tasks.
 */
class LoadModel {
 public:
  virtual ~LoadModel() = default;

  /**
   * Calculates the relative expected load of a given task compared to other 
   * tasks.
   *
   * The calculated load is in no specific unit. However, the 
   * "greater than"-relation for two tasks has to be correct.
   *
   * @param lvlVec The level vector corresponding to the task.
   * @returns Relative value (in comparison of different tasks) of the expected
   *          load for the given task.
   */
  virtual real eval(const LevelVector& lvlVec) = 0;
};

} /* namespace combigrid */
#endif /* LOADMODEL_HPP_ */
