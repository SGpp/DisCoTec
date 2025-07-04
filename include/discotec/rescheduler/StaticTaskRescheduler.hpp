#ifndef STATICTASKRESCHEDULER_CPP_
#define STATICTASKRESCHEDULER_CPP_

#include "TaskRescheduler.hpp"

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
   * @param levelVectorToProcessGroupIndex The current task distribution -- 
   *                                       ignored.
   * @param levelVectorToTaskDuration The last measured durations of the tasks 
   *                                  run function -- ignored.
   * @param loadModel The load model to use for the prognosis of future task 
   *                  loads -- ignored.
   * @returns An empty vector.
   */
  std::vector<std::pair<LevelVector, int>> eval(
      const std::map<LevelVector, int>& levelVectorToProcessGroupIndex,
      const std::map<LevelVector, unsigned long>& levelVectorToTaskDuration,
      LoadModel *loadModel) override {
    return {};
  }
};

} /* namespace combigrid */

#endif /* STATICTASKRESCHEDULER_CPP_ */
