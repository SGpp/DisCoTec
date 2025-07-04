#ifndef LEARNINGLOADMODEL_HPP_
#define LEARNINGLOADMODEL_HPP_

#include "LoadModel.hpp"
#include "../utils/LevelVector.hpp"
#include "../utils/Config.hpp"

namespace combigrid {

/**
 * POD containing fields of relevant values for the calculation of expected load.
 * This POD only contains fields that change over the runtime of a simulation.
 */
struct DurationInformation {
  size_t task_id; /**< included for easy identification of the corresponding task */
  unsigned long duration;
  real simtime_now;
  real real_dt;
  int pgroup_id;
  unsigned int nProcesses;

  /**
   * Function to serialize the struct. See "Boost Serialization".
   */
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar & task_id;
    ar & duration;
    ar & simtime_now;
    ar & real_dt;
    ar & pgroup_id;
    ar & nProcesses;
  }
};

/**
 * The LearningLoadModel extends the interface of a LoadModel. 
 * A LearningLoadModel is able to receive information about tasks that can be 
 * used to calculate a more accurate expectation of load.
 */
class LearningLoadModel : public LoadModel {
 public:

  /**
   * Adds duration information about a task.
   *
   * @param info The duration information to add.
   * @param lvlVec The level vector of the task to add the information to.
   */
  virtual void addDurationInformation(const DurationInformation& info, 
                                      const LevelVector& lvlVec) = 0;
};

} /* namespace combigrid */

#endif /* LEARNINGLOADMODEL_HPP_ */
