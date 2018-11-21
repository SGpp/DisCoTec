/*
 * LearningLoadModel.hpp
 *
 *      Author: pollinta
 */

#ifndef LEARNINGLOADMODEL_HPP_
#define LEARNINGLOADMODEL_HPP_

#include <string>
#include <deque>     

#include "sgpp/distributedcombigrid/loadmodel/LoadModel.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LinearLoadModel.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/utils/Stats.hpp"

namespace combigrid {
    //TODO include "metadata" in model: nrg.dat, parameters etc. ?
    struct durationInformation{
      // MPI_INT task_id;
      int task_id;
      // MPI_UNSIGNED_LONG duration;
      long unsigned int duration;
      // MPI_UNSIGNED nProcesses;
      uint nProcesses;
      // long int order;
    };

    //this data type needs to be communicated to MPI in every process using it
    MPI_Datatype createMPIDurationType();

    class DurationType{
      public:
        DurationType(){
          durationType_ = createMPIDurationType();
        }
        ~DurationType(){
          MPI_Type_free( &durationType_ );
        }
        MPI_Datatype get() const {
          return durationType_;
        }
      private:
        MPI_Datatype durationType_;
    };

    std::string getFilename(const LevelVector& levelVector);

    unsigned long int getAverageOfFirstColumn(std::istream& fs);

class LearningLoadModel : public LoadModel {

 public:
  LearningLoadModel(std::vector<LevelVector> levelVectors){
    durationsOfLevels_ = std::unique_ptr<std::map<LevelVector,std::deque<std::pair<long unsigned int, uint>>>>\
      (new std::map<LevelVector,std::deque<std::pair<long unsigned int, uint>>>);
    for (const auto & lv : levelVectors){
      durationsOfLevels_->insert(make_pair(lv, std::deque<std::pair<long unsigned int, uint>>()));
    }

    for (const auto & lv: levelVectors){
      std::fstream fs;
      // files_.insert(make_pair(lv, std::unique_ptr<std::fstream>(new std::fstream())));
      fs.open(getFilename(lv), std::fstream::in | std::fstream::app);
      fs.seekg(0, std::ios::beg);
      fs.clear();
      if (fs.peek() != EOF){
        last_durations_avg_.insert(make_pair(lv,getAverageOfFirstColumn(fs)));
      }
      fs.clear();
      fs.close();
    }
  }

  inline real eval(const LevelVector& l);

  ~LearningLoadModel() {
    //upon destruction, write data gathered to file
    for (const auto& item : *durationsOfLevels_){
      std::fstream fs;
      LevelVector lv = item.first;
      fs.open(getFilename(lv), std::fstream::out | std::fstream::app);
      // assert(fs.good());
      for (const auto& duration : durationsOfLevels_->at(lv)){
        std::cout << "writing data point for " << toString(lv) << " to load model" << std::endl;
        fs << std::to_string(duration.first) << ", " << std::to_string(duration.second) << ", " << std::endl;
      }
      // std::cout << "open2: " << std::to_string(fs.is_open()) << std::endl;
      fs.close();
    }
  }

  void addDataPoint(durationInformation dI, LevelVector l){
    // std::cout << "adding data point for " << toString(l) << " to load model" << std::endl;
    addDataPoint(dI.duration, dI.nProcesses, l);
  }

  void addDataPoint(long unsigned int duration, uint nProcesses, LevelVector l){ //TODO write to file intermediately for long simulations
    std::pair<long unsigned int, uint> value = std::make_pair(duration, nProcesses);
    durationsOfLevels_->at(l).push_back(value);
  }
 
 private:
  std::unique_ptr<std::map<LevelVector,std::deque<std::pair<long unsigned int, uint>>>> durationsOfLevels_;
  // std::map<LevelVector, std::unique_ptr<std::fstream>> files_;
  std::map<LevelVector, long unsigned int> last_durations_avg_;
};


inline real LearningLoadModel::eval(const LevelVector& l) {
  // std::cout << "eval loadmodel at " << toString(l) << std::endl;
  real ret(0.0); 
  if (durationsOfLevels_->at(l).empty()){
    //if no data yet, use linear load model
    if ( ! last_durations_avg_.count(l) ){
      std::cout << "using linear load model for initial guess" << std::endl;
      LinearLoadModel llm = LinearLoadModel();
      ret = llm.eval(l);
    }else{
      // use data from last time
      std::cout << "using recoded average of " << last_durations_avg_[l] << " for " << toString(l) << std::endl;
      ret = last_durations_avg_[l];
    }
  }else{// use simple averaging for now //TODO do fancier things
    for (auto duration : durationsOfLevels_->at(l)){
      ret += duration.first;
    }
    ret /= static_cast<double>(durationsOfLevels_->at(l).size());
  }

  return ret;
}

} /* namespace combigrid */

#endif /* LEARNINGLOADMODEL_HPP_ */
