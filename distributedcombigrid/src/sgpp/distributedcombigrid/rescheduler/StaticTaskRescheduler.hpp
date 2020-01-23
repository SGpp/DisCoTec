#ifndef STATICTASKRESCHEDULER_CPP_
#define STATICTASKRESCHEDULER_CPP_

#include "sgpp/distributedcombigrid/rescheduler/TaskRescheduler.hpp"

namespace combigrid {

/**
 * Static task rescheduler implementing the task rescheduler interface.
 *
 * The static task rescheduler NEVER instructs a change in the task distribution.
 */
class StaticTaskRescheduler : public TaskRescheduler {
 public:
  /**
   * Does not calculate a change in task distribution. It never instructs a 
   * change in the task distribution.
   *
   * @param levelVectorToProcessGroupIndex The current task distribution.
   * @param levelVectorToTaskDuration The last measured durations of the tasks 
   *                                  run function.
   * @returns An empty vector.
   */
  virtual std::vector<std::pair<LevelVector, int>> eval(
      const std::map<LevelVector, int>& levelVectorToProcessGroupIndex,
      const std::map<LevelVector, unsigned long>& levelVectorToTaskDuration) {
    return {};
  }
};

} /* namespace combigrid */

#endif /* STATICTASKRESCHEDULER_CPP_ */
