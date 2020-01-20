#ifndef SPECIFICUOTLEARNINGLOADMODEL_HPP_
#define SPECIFICUOTLEARNINGLOADMODEL_HPP_

#include <chrono>

#include "sgpp/distributedcombigrid/loadmodel/LearningLoadModel.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"

namespace combigrid {

/**
 * The SpecificUOTLearningModel (UOT = unit of time) extends the interface of a 
 * LearningLoadModel.
 * A SpecificUOTLearningLoadModel needs to be able to calculate the expected 
 * load in the specific time unit of microseconds.
 */
class SpecificUOTLearningLoadModel : public LearningLoadModel {
 public:
  virtual ~SpecificUOTLearningLoadModel() = default;

   /**
    * Calculates the expected load in microseconds of a given task.
    *
    * @param lvlVec The level vector corresponding to the task.
    * @returns Expected load of the given task.
    */
   virtual std::chrono::microseconds evalSpecificUOT(const LevelVector& lvlVec) = 0;

};

} /* namespace combigrid */

#endif /* SPECIFICUOTLEARNINGLOADMODEL_HPP_ */

