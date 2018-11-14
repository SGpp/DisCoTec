/*
 * LearningLoadModel.hpp
 *
 *      Author: pollinta
 */

#ifndef LEARNINGLOADMODEL_HPP_
#define LEARNINGLOADMODEL_HPP_

// #include <iostream>
#include <string>
#include <chrono>
// #include <memory>

#include "sgpp/distributedcombigrid/loadmodel/LoadModel.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LinearLoadModel.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/utils/Stats.hpp"

namespace combigrid {
    //TODO include "metadata" in model: nrg.dat, parameters etc.
    struct durationInformation{
      // MPI_LONG duration;
      long int duration; //todo make larger data type?
      // MPI_UNSIGNED nProcesses;
      uint nProcesses;
      // long int order;
    };

    //this data type needs to be communicated to MPI in every process using it
    MPI_Datatype createMPIDurationType();


class LearningLoadModel : public LoadModel {

 public:
  LearningLoadModel(std::vector<LevelVector> levelVectors){
  }

  inline real eval(const LevelVector& l);

  virtual ~LearningLoadModel() = default;
 
 private:
};


inline real LearningLoadModel::eval(const LevelVector& l) {
  real ret(0.0);
  //if no data yet, use linear load model //TODO
  LinearLoadModel llm = LinearLoadModel();
  ret = llm.eval(l);

  return ret;
}

} /* namespace combigrid */

#endif /* LEARNINGLOADMODEL_HPP_ */
