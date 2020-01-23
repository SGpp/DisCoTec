#ifndef REBALANCINGTASKRESCHEDULER_HPP_
#define REBALANCINGTASKRESCHEDULER_HPP_

#include <memory>

#include "sgpp/distributedcombigrid/rescheduler/TaskRescheduler.hpp"
#include "sgpp/distributedcombigrid/loadmodel/SpecificUOTLearningLoadModel.hpp"

namespace combigrid {

/**
 * A task rescheduler improving the load balance by using a load model.
 * The resulting load balance is not necessarily optimal.
 */
class RebalancingTaskRescheduler : public TaskRescheduler {
 public:
  /**
   * Constructor for a rebalancing task rescheduler with the specific UOT load 
   * model to use.
   *
   * @param sllm The load model to use for prognosis of future task loads.
   * @param max_iterations The maximum of allowed iterations. In every 
   *                       iteration a maximum of one task is going to be able 
   *                       to be rescheduled.
   * @param min_inbalance The minimum of load inbalance necessary to start the 
   *                      search for a task to reschedule.
   */
  RebalancingTaskRescheduler(std::shared_ptr<SpecificUOTLearningLoadModel> sllm,
                             unsigned int max_iterations,
                             double min_inbalance) : 
    sllm_{sllm}, 
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
   * @returns Vector of pair with level vector and process group. To realize 
   *          the calculated optimized task distribution; the task (defined by 
   *          its level vector) has to be moved to the process group (defined 
   *          by its process group id).
   */
  virtual std::vector<std::pair<LevelVector, int>> eval(
      const std::map<LevelVector, int>& levelVectorToProcessGroupIndex,
      const std::map<LevelVector, unsigned long>& levelVectorToTaskDuration);

 private:
  std::shared_ptr<SpecificUOTLearningLoadModel> sllm_;
  unsigned int max_iterations_;
  double min_inbalance_;
};

} /* namespace combigrid */

#endif /* REBALANCINGTASKRESCHEDULER_HPP_ */



