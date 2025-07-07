#ifndef AVERAGEOFLASTNLOADMODEL_HPP_
#define AVERAGEOFLASTNLOADMODEL_HPP_

#include <map>
#include <deque>

#include "MicrosecondsLearningLoadModel.hpp"
#include "LinearLoadModel.hpp"

namespace combigrid {

/**
 * Load model implementing a specific UOT learning load model.
 * The expected task load is the average load of the last n run calculations.
 */
class AverageOfLastNLoadModel : public MicrocsecondsLearningLoadModel {
 public:
  /**
   * The constructor for this load model.
   *
   * @param n Count of past run times to be used for the calculation of the
   *          load prognosis.
   * @param tasks The level vectors of all tasks in this simulation.
   * @param loadModelIfNoHistory The load model used for the calculation of
   *                             LoadModel::eval when no past duration
   *                             information is available.
   */
  AverageOfLastNLoadModel(unsigned int n, const std::vector<LevelVector>& tasks,
                          std::unique_ptr<LoadModel> loadModelIfNoHistory);

  /**
   * Adds duration information about a task.
   *
   * @param info The duration information to add.
   * @param lvlVec The level vector of the task to add the information to.
   */
  void addDurationInformation(const DurationInformation& info,
                              const LevelVector& lvlVec) override;

  /**
   * Calculates the expected load of a given task in microseconds by using the
   * average of the last n received duration informations.
   *
   * @param lvlVec The level vector corresponding to the task.
   * @returns Expected load of the given task. If no accurate calculation is
   *          possible zero seconds are returned.
   */
  std::chrono::microseconds evalSpecificUOT(const LevelVector& lvlVec) override;

  /**
   * Calculates the relative expected load of a given task compared to other
   * tasks.
   *
   * @param lvlVec The level vector corresponding to the task.
   * @returns Relative value of expected load for the given task.
   */
  real eval(const LevelVector& lvlVec) override;

 private:
  // TODO Optimization possibility: Instead of a deque use a circular buffer.
  std::map<LevelVector, std::deque<unsigned long>> levelVectorToLastNDurations_;
  unsigned int lastN_;
  std::unique_ptr<LoadModel> loadModelIfNoHistory_;
};

} /* namespace combigrid */

#endif /* AVERAGEOFLASTNLOADMODEL_HPP_ */
