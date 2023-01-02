#ifndef REBALANCINGTASKRESCHEDULER_HPP_
#define REBALANCINGTASKRESCHEDULER_HPP_

#include <memory>

#include "rescheduler/TaskRescheduler.hpp"
#include "loadmodel/MicrosecondsLearningLoadModel.hpp"

namespace combigrid {

/**
 * A task rescheduler improving the load balance by using a load model.
 * The resulting load balance is not necessarily optimal.
 */
class RebalancingTaskRescheduler : public TaskRescheduler {
 public:
  /**
   * Constructor for a rebalancing task rescheduler with specific scheduling 
   * parameters.
   *
   * @param max_iterations The maximum of allowed iterations. In every 
   *                       iteration a maximum of one task is going to be able 
   *                       to be rescheduled.
   * @param min_inbalance The minimum of load inbalance necessary to start the 
   *                      search for a task to reschedule.
   */
  RebalancingTaskRescheduler(unsigned int max_iterations,
                             double min_inbalance) : 
    max_iterations_{max_iterations}, 
    min_inbalance_{min_inbalance} 
  {};

  /**
   * Calculates changes in task distribution with the goal of a more balanced
   * runtime of process groups.
   *
   * The resulting load balance is not necessarily going to be the optimal 
   * distribution.
   *
   * @param levelVectorToProcessGroupIndex The current task distribution.
   * @param levelVectorToTaskDuration The last measured durations of the tasks 
   *                                  run function.
   * @param loadModel The load model to use for the prognosis of future task 
   *                  loads. Needs to be a instance of SpecificUOTLearningLoadModel 
   *                  otherwise the scheduling is aborted.
   * @returns Vector of pairs of level vector and process group. 
   *          The semantic meaning is: The task (defined by its level vector) 
   *          should be moved to the process group (defined by its process 
   *          group id).
   */
  std::vector<std::pair<LevelVector, int>> eval(
      const std::map<LevelVector, int>& levelVectorToProcessGroupIndex,
      const std::map<LevelVector, unsigned long>& levelVectorToTaskDuration,
      LoadModel *loadModel) override;

 private:
  unsigned int max_iterations_;
  double min_inbalance_;
};

} /* namespace combigrid */

#endif /* REBALANCINGTASKRESCHEDULER_HPP_ */



