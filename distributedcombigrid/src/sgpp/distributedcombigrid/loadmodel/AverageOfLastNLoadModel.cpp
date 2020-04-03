#include "sgpp/distributedcombigrid/loadmodel/AverageOfLastNLoadModel.hpp"

#include <numeric>

namespace combigrid {

AverageOfLastNLoadModel::AverageOfLastNLoadModel(
    unsigned int n, const std::vector<LevelVector>& tasks, 
    std::unique_ptr<LoadModel> loadModelIfNoHistory)
    : lastN_{n}, loadModelIfNoHistory_{std::move(loadModelIfNoHistory)} {

  for (const auto& l : tasks) {
    this->levelVectorToLastNDurations_.insert({l, {}});
  }
}

void AverageOfLastNLoadModel::addDurationInformation(
    const DurationInformation& info, const LevelVector& lvlVec) {

  auto& durations = this->levelVectorToLastNDurations_[lvlVec];
  if (durations.size() >= this->lastN_) {
    durations.pop_front();
  }
  durations.push_back(info.duration);
}

std::chrono::microseconds AverageOfLastNLoadModel::evalSpecificUOT(
    const LevelVector& lvlVec) {

  const auto& durations = this->levelVectorToLastNDurations_[lvlVec];
  const auto average = (static_cast<real>(std::accumulate(durations.begin(), 
                                                          durations.end(), 0)) 
                        / static_cast<real>(durations.size()));
  return std::chrono::microseconds{static_cast<long>(average)};
}

real AverageOfLastNLoadModel::eval(const LevelVector& lvlVec) {
  if (this->levelVectorToLastNDurations_.at(lvlVec).size() == 0) {
    return this->loadModelIfNoHistory_->eval(lvlVec);
  } else {
    return this->evalSpecificUOT(lvlVec).count();
  }
}

} /* namespace combigrid */
