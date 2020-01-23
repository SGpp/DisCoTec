#ifndef TASKRESCHEDULER_HPP_
#define TASKRESCHEDULER_HPP_

#include <vector>
#include <map>

#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"

namespace combigrid {

/**
 * Interface of a task re-scheduler.
 */
class TaskRescheduler {
 public:
  /**
   * Calculates changes in task distribution with the goal of a more balanced
   * runtime of process groups.
   *
   * @param levelVectorToProcessGroupIndex The current task distribution.
   * @param levelVectorToTaskDuration The last measured durations of the tasks 
   *                                  run function.
   * @returns Vector of pair with level vector and process group. To realize 
   *          the calculated optimized task distribution; the task (defined by 
   *          its level vector) has to be moved to the process group (defined 
   *          by its process group id).
   */
  virtual std::vector<std::pair<LevelVector, int>> eval(
      const std::map<LevelVector, int>& levelVectorToProcessGroupIndex,
      const std::map<LevelVector, unsigned long>& levelVectorToTaskDuration) = 0;
};

} /* namespace combigrid */

#endif /* TASKRESCHEDULER_HPP_ */
