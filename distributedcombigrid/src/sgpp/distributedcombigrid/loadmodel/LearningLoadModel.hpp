#ifndef LEARNINGLOADMODEL_HPP_
#define LEARNINGLOADMODEL_HPP_

#include "sgpp/distributedcombigrid/loadmodel/LoadModel.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/utils/Config.hpp"

namespace combigrid {

/**
 * POD containing fields of relevant values for the calculation of expected load.
 * This POD only contains fields that change over the runtime of the simulation.
 */
struct DurationInformation {
  long unsigned int duration;
  real simtime_now;
  real real_dt;
  int pgroup_id;
  uint nProcesses;

  /**
   * Function to serialize the struct into the given archive. It is used to 
   * send this POD via an MPI message.
   */
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
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
   * @param info the duration information to add
   * @param lvlVec the level vector of the task to add the information to
   */
  virtual void addDurationInformation(DurationInformation info, 
                                      LevelVector lvlVec) = 0;
};

} /* namespace combigrid */

#endif /* LEARNINGLOADMODEL_HPP_ */
