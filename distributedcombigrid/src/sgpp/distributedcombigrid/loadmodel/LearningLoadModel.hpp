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
    //TODO include "metadata" in model: nrg.dat, parameters etc.
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

class LearningLoadModel : public LoadModel {

 public:
  LearningLoadModel(std::vector<LevelVector> levelVectors){
    durationsOfLevels_ = std::unique_ptr<std::map<LevelVector,std::deque<std::pair<long unsigned int, uint>>>>\
      (new std::map<LevelVector,std::deque<std::pair<long unsigned int, uint>>>);
    for (const auto & lv : levelVectors){
      durationsOfLevels_->insert(make_pair(lv, std::deque<std::pair<long unsigned int, uint>>()));
    }

    for (const auto & lv: levelVectors){
      files_.insert(make_pair(lv, std::unique_ptr<std::fstream>(new std::fstream())));
      files_[lv]->open(getFilename(lv), std::fstream::in | std::fstream::out | std::fstream::app);
    }
  }

  inline real eval(const LevelVector& l);

  ~LearningLoadModel() {
    //upon destruction, write data gathered to file //TODO
    for (const auto& item : files_){
      LevelVector lv = item.first;
      for (const auto& duration : durationsOfLevels_->at(lv)){
        *(files_[lv]) << std::to_string(duration.first) << ", " << std::to_string(duration.second) << ", " << std::endl;
      }
      // std::cout << "open2: " << std::to_string(files_[lv]->is_open()) << std::endl;
      files_[lv]->close();
    }
  }

  void addDataPoint(durationInformation dI, LevelVector l){
    addDataPoint(dI.duration, dI.nProcesses, l);
  }

  void addDataPoint(long unsigned int duration, uint nProcesses, LevelVector l){
    std::pair<long unsigned int, uint> value = std::make_pair(duration, nProcesses);
    durationsOfLevels_->at(l).push_back(value);
  }
 
 private:
  std::unique_ptr<std::map<LevelVector,std::deque<std::pair<long unsigned int, uint>>>> durationsOfLevels_;
  std::map<LevelVector, std::unique_ptr<std::fstream>> files_;
};


inline real LearningLoadModel::eval(const LevelVector& l) {
  real ret(0.0);
  //if no data yet, use linear load model //TODO read data from last time
  if (durationsOfLevels_->at(l).empty()){
    LinearLoadModel llm = LinearLoadModel();
    ret = llm.eval(l);
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
