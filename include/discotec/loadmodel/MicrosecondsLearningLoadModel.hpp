#ifndef MICROSECONDSLEARNINGLOADMODEL_HPP_
#define MICROSECONDSLEARNINGLOADMODEL_HPP_

#include "../../../../../.cache/JetBrains/CLion2024.1/.remote/ipvs-epyc2_22/347dacd9-0950-4b02-b7f7-b2d9ba9d9d1d/usr/include/c++/9/chrono"

#include "LearningLoadModel.hpp"
#include "../utils/LevelVector.hpp"

namespace combigrid {

/**
 * The SpecificUOTLearningModel (UOT = unit of time) extends the interface of a 
 * LearningLoadModel.
 * A SpecificUOTLearningLoadModel needs to be able to calculate the expected 
 * load in the specific time unit of microseconds.
 */
class MicrocsecondsLearningLoadModel : public LearningLoadModel {
 public:

  /**
   * Calculates the expected load of a given task in microseconds.
   *
   * @param lvlVec The level vector corresponding to the task.
   * @returns Expected load of the given task. For consistency; if no accurate 
   *          calculation is possible return 0 seconds.
   */
  virtual std::chrono::microseconds evalSpecificUOT(const LevelVector& lvlVec) = 0;
};

} /* namespace combigrid */

#endif /* MICROSECONDSLEARNINGLOADMODEL_HPP_ */

