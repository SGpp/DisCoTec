#include "sgpp/distributedcombigrid/loadmodel/AverageOfLastNLoadModel.hpp"

#include <numeric>

namespace combigrid {

template<typename LM>
AverageOfLastNLoadModel<LM>::AverageOfLastNLoadModel(
    unsigned int n, std::vector<LevelVector> tasks, LM loadModelIfNoHistory)
    : lastN_{n}, loadModelIfNoHistory_{loadModelIfNoHistory} {

  for (const auto& l : tasks) {
    this->levelVectorToLastNDurations_.insert({l, {}});
  }
}

template<typename LM>
void AverageOfLastNLoadModel<LM>::addDurationInformation(
    DurationInformation info, LevelVector lvlVec) {

  auto& durations = this->levelVectorToLastNDurations_[lvlVec];
  if (durations.size() >= this->n_) {
    durations.pop_front();
  }
  durations.push_back(info.duration);
}

template<typename LM>
std::chrono::microseconds AverageOfLastNLoadModel<LM>::evalSpecificUOT(
    const LevelVector& lvlVec) {

  const auto& durations = this->levelVectorToLastNDurations_[lvlVec];
  const auto average = (static_cast<real>(std::accumulate(durations.begin(), 
                                                          durations.end(), 0)) 
                        / static_cast<real>(durations.size()));
  return std::chrono::microseconds{static_cast<long>(average)};
}

template<typename LM>
real AverageOfLastNLoadModel<LM>::eval(const LevelVector& lvlVec) {
  if (this->levelVectorToLastNDurations_.at(lvlVec).size() == 0) {
    return this->loadModelIfNoHistory_.eval(lvlVec);
  } else {
    return this->evalSpecificUOT(lvlVec);
  }
}

} /* namespace combigrid */
