/*
 * LearningLoadModel.hpp
 *
 *      Author: pollinta
 */

#ifndef LEARNINGLOADMODEL_HPP_
#define LEARNINGLOADMODEL_HPP_

#include <deque>
#include <limits>
#include <string>

#include "sgpp/distributedcombigrid/loadmodel/LinearLoadModel.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LoadModel.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/utils/Stats.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/mpi/MPIUtils.hpp"

namespace combigrid {
// TODO include "metadata" in model: nrg.dat, parameters etc. ?
struct durationInformation {
  int task_id;
  long unsigned int duration;
  real real_dt;
  int pgroup_id;
  uint nProcesses;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& task_id;
    ar& duration;
    ar& real_dt;
    ar& pgroup_id;
    ar& nProcesses;
  }
};


static std::string getFilename(const LevelVector& levelVector){
  return "./l_" + toString(levelVector) + ".durations";//TODO which directory
}

unsigned long int getAverageOfFirstColumn(std::istream& fs);

class LearningLoadModel : public LoadModel {
 public:
  LearningLoadModel(std::vector<LevelVector> levelVectors) :
  writeEveryNCombis_(20)
  {
    durationsOfLevels_ =
        std::unique_ptr<std::map<LevelVector, std::deque<durationInformation>>>(
            new std::map<LevelVector, std::deque<durationInformation>>);
    for (const auto& lv : levelVectors) {
      addLevelVectorToLoadModel(lv);
    }
  }

  inline real eval(const LevelVector& l);

  ~LearningLoadModel() {
    // upon destruction, write data gathered to file
    for (const auto& item : *durationsOfLevels_) {
      LevelVector lv = item.first;
      writeDurationsToFile(lv);
    }
  }

  void addDataPoint(durationInformation dI, LevelVector l) {
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
                            std::deque<durationInformation>& durations,
                            size_t num_entries = std::numeric_limits<size_t>::max()) {
    std::fstream fs;
    fs.open(getFilename(lv), std::fstream::out | std::fstream::app);
    // std::cout << "writing data points for " << toString(lv) << " to load model" << std::endl;
    // assert(fs.good());
    for (size_t i = 0; i < num_entries && !durations.empty(); ++i) {
      durationInformation duration = durations.front();
      fs << std::to_string(duration.duration) << " " << std::to_string(duration.real_dt) << " " 
        //  << std::to_string(duration.task_id)  << " " 
        << std::to_string(duration.pgroup_id) << " " << std::to_string(duration.nProcesses) << " "
         << std::endl;
      durations.pop_front();
    }
    fs.close();
  }

  // write all durations for lv to file
  void writeDurationsToFile(const LevelVector& lv) {
    writeDurationsToFile(lv, durationsOfLevels_->at(lv));
  }

  std::unique_ptr<std::map<LevelVector, std::deque<durationInformation>>>
      durationsOfLevels_;
  std::map<LevelVector, long unsigned int> last_durations_avg_;
};

inline real LearningLoadModel::eval(const LevelVector& l) {
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

#endif /* LEARNINGLOADMODEL_HPP_ */
