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

namespace combigrid {
// TODO include "metadata" in model: nrg.dat, parameters etc. ?
struct durationInformation {
  // MPI_INT task_id;
  int task_id;
  // MPI_UNSIGNED_LONG duration;
  long unsigned int duration;
  // MPI_UNSIGNED nProcesses;
  uint nProcesses;
  // long int order;
};

// this data type needs to be communicated to MPI in every process using it
MPI_Datatype createMPIDurationType();

class DurationType {
 public:
  DurationType() { durationType_ = createMPIDurationType(); }
  ~DurationType() { MPI_Type_free(&durationType_); }
  MPI_Datatype get() const { return durationType_; }

 private:
  MPI_Datatype durationType_;
};

std::string getFilename(const LevelVector& levelVector);

unsigned long int getAverageOfFirstColumn(std::istream& fs);

class LearningLoadModel : public LoadModel {
 public:
  LearningLoadModel(std::vector<LevelVector> levelVectors) {
    durationsOfLevels_ =
        std::unique_ptr<std::map<LevelVector, std::deque<std::pair<long unsigned int, uint>>>>(
            new std::map<LevelVector, std::deque<std::pair<long unsigned int, uint>>>);
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
    // std::cout << "adding data point for " << toString(l) << " to load model" << std::endl;
    addDataPoint(dI.duration, dI.nProcesses, l);
  }

  void addDataPoint(long unsigned int duration, uint nProcesses, LevelVector l) {
    std::pair<long unsigned int, uint> value = std::make_pair(duration, nProcesses);
    durationsOfLevels_->at(l).push_back(value);
    // write to file intermediately for long simulations
    int maxsize = 500;
    if (durationsOfLevels_->at(l).size() > maxsize) {
      writeDurationsToFile(l, durationsOfLevels_->at(l), 200);
    }
    assert(durationsOfLevels_->at(l).size() < maxsize + 1);
  }

 private:
  void addLevelVectorToLoadModel(const LevelVector& lv) {
    durationsOfLevels_->insert(make_pair(lv, std::deque<std::pair<long unsigned int, uint>>()));

    std::fstream fs;
    fs.open(getFilename(lv), std::fstream::in | std::fstream::app);
    fs.seekg(0, std::ios::beg);
    fs.clear();
    if (fs.peek() != EOF) {
      last_durations_avg_.insert(make_pair(lv, getAverageOfFirstColumn(fs)));
    }
    fs.clear();
    fs.close();
  }

  // durations that are written to file will be removed from the referenced deque
  void writeDurationsToFile(const LevelVector& lv,
                            std::deque<std::pair<long unsigned int, uint>>& durations,
                            int num_entries = std::numeric_limits<int>::max()) {
    std::fstream fs;
    fs.open(getFilename(lv), std::fstream::out | std::fstream::app);
    // std::cout << "writing data points for " << toString(lv) << " to load model" << std::endl;
    // assert(fs.good());
    for (int i = 0; i < num_entries && !durations.empty(); ++i) {
      std::pair<long unsigned int, uint> duration = durations.front();
      fs << std::to_string(duration.first) << ", " << std::to_string(duration.second) << ", "
         << std::endl;
      durations.pop_front();
    }
    fs.close();
  }

  // write all durations for lv to file
  void writeDurationsToFile(const LevelVector& lv) {
    writeDurationsToFile(lv, durationsOfLevels_->at(lv));
  }

  std::unique_ptr<std::map<LevelVector, std::deque<std::pair<long unsigned int, uint>>>>
      durationsOfLevels_;
  std::map<LevelVector, long unsigned int> last_durations_avg_;
};

inline real LearningLoadModel::eval(const LevelVector& l) {
  // std::cout << "eval loadmodel at " << toString(l) << std::endl;
  real ret(0.0);
  if (durationsOfLevels_->at(l).empty()) {
    // if no data yet, use linear load model
    if (!last_durations_avg_.count(l)) {
      std::cout << "using linear load model for initial guess" << std::endl;
      LinearLoadModel llm = LinearLoadModel();
      ret = llm.eval(l);
    } else {
      // use data from last time
      std::cout << "using recoded average of " << last_durations_avg_[l] << " for " << toString(l)
                << std::endl;
      ret = last_durations_avg_[l];
    }
  } else {  // use simple averaging for now //TODO do fancier things, cache results
    for (auto duration : durationsOfLevels_->at(l)) {
      ret += duration.first;
    }
    ret /= static_cast<double>(durationsOfLevels_->at(l).size());
  }

  return ret;
}

} /* namespace combigrid */

#endif /* LEARNINGLOADMODEL_HPP_ */
