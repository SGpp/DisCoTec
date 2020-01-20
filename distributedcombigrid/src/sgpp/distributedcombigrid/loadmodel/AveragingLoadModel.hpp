#ifndef AVERAGINGLOADMODEL_HPP_
#define AVERAGINGLOADMODEL_HPP_

#include <deque>
#include <limits>
#include <string>

#include "sgpp/distributedcombigrid/loadmodel/LearningLoadModel.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LinearLoadModel.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/utils/Stats.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/mpi/MPIUtils.hpp"

namespace combigrid {

static std::string getFilename(const LevelVector& levelVector){
  return "./l_" + toString(levelVector) + ".durations";//TODO which directory
}

unsigned long int getAverageOfFirstColumn(std::istream& fs);

/**
 * AveragingLoadModel is a LearningLoadModel that is saving received duration 
 * information to files. For its expected load calculation it is using the 
 * average of previously received duration information (from files). The 
 * expected load calculation of the LinearLoadModel is used if no previous 
 * duration information is available.
 */
class AveragingLoadModel : public LearningLoadModel {
 public:
  AveragingLoadModel(std::vector<LevelVector> levelVectors) :
  writeEveryNCombis_(50)
  {
    durationsOfLevels_ =
        std::unique_ptr<std::map<LevelVector, std::deque<DurationInformation>>>(
            new std::map<LevelVector, std::deque<DurationInformation>>);
    for (const auto& lv : levelVectors) {
      addLevelVectorToLoadModel(lv);
    }
  }

  inline real eval(const LevelVector& l);

  ~AveragingLoadModel() {
    // upon destruction, write data gathered to file
    for (const auto& item : *durationsOfLevels_) {
      LevelVector lv = item.first;
      writeDurationsToFile(lv);
    }
  }

  void addDurationInformation(DurationInformation dI, LevelVector l) {
    durationsOfLevels_->at(l).push_back(dI);
    // write to file intermediately for long simulations
    // if (durationsOfLevels_->at(l).size() > writeEveryNCombis_*2) {
    if (durationsOfLevels_->at(l).size() > writeEveryNCombis_) {
      writeDurationsToFile(l, durationsOfLevels_->at(l), writeEveryNCombis_);
    }
    assert(durationsOfLevels_->at(l).size() < writeEveryNCombis_ + 1);
  }

 private:
  size_t writeEveryNCombis_;

  void addLevelVectorToLoadModel(const LevelVector& lv);

  // durations that are written to file will be removed from the referenced deque
  void writeDurationsToFile(const LevelVector& lv,
                            std::deque<DurationInformation>& durations,
                            size_t num_entries = std::numeric_limits<size_t>::max()) {
    std::fstream fs;
    fs.open(getFilename(lv), std::fstream::out | std::fstream::app);
    // std::cout << "writing data points for " << toString(lv) << " to load model" << std::endl;
    // assert(fs.good());
    for (size_t i = 0; i < num_entries && !durations.empty(); ++i) {
      DurationInformation duration = durations.front();
      fs << std::to_string(duration.duration) << " " 
         << std::to_string(duration.real_dt) << " " 
         // << std::to_string(duration.task_id)  << " "
         << std::to_string(duration.simtime_now)  << " " 
         << std::to_string(duration.pgroup_id) << " "
         << std::endl;
      durations.pop_front();
    }
    fs.close();
  }

  // write all durations for lv to file
  void writeDurationsToFile(const LevelVector& lv) {
    writeDurationsToFile(lv, durationsOfLevels_->at(lv));
  }

  std::unique_ptr<std::map<LevelVector, std::deque<DurationInformation>>>
      durationsOfLevels_;
  std::map<LevelVector, long unsigned int> last_durations_avg_;
};

inline real AveragingLoadModel::eval(const LevelVector& l) {
  // std::cout << "eval loadmodel at " << toString(l) << std::endl;
  real ret(0.0);
  if (durationsOfLevels_->at(l).empty()) {
    // if no data yet, use linear load model
    if (!last_durations_avg_.count(l)) {
      // std::cout << "using linear load model for initial guess" << std::endl;
      LinearLoadModel llm = LinearLoadModel();
      ret = llm.eval(l);
    } else {
      // use data from last time
      // std::cout << "using recoded average of " << last_durations_avg_[l] << " for " << toString(l)
                // << std::endl;
      ret = last_durations_avg_[l];
    }
  } else {  // use simple averaging for now //TODO do fancier things, cache results
    for (auto duration : durationsOfLevels_->at(l)) {
      ret += duration.duration;
    }
    ret /= static_cast<double>(durationsOfLevels_->at(l).size());
  }

  return ret;
}

} /* namespace combigrid */

#endif /* AVERAGINGLOADMODEL_HPP_ */
