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
  namespace durationsFile{  

    //TODO include "metadata" in model: nrg.dat, parameters etc.
    struct durationInformation{
      // MPI_LONG duration;
      long int duration; //todo make larger data type
      // MPI_UNSIGNED nProcesses;
      uint nProcesses;
      // long int order;
    };

    //this data type needs to be communicated to MPI in every process using it
    MPI_Datatype createMPIDurationType();

    std::string getFilename(const LevelVector& levelVector);

    // each of the performance files will be opened two times: once to write, from the worker process
    class DurationsWriteFile{
    public:
      DurationsWriteFile(const LevelVector& levelVector) : DurationsWriteFile(getFilename(levelVector)) {}

      void write( durationInformation* buf, int bufsize = 1){
        MPI_File_write(fh_, buf, bufsize, durationType_, MPI_STATUS_IGNORE); 
      }

      void write(const Stats::Event event, size_t nProcesses);
      
      ~DurationsWriteFile(){
        MPI_File_close(&fh_); 
      }

    private:
      DurationsWriteFile(std::string filename){
        durationType_ = createMPIDurationType();
        MPI_File_open(
          MPI_COMM_SELF, filename.c_str(), 
          MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND,// |  MPI_MODE_SEQUENTIAL, 
          MPI_INFO_NULL, &fh_
        ); 
      }
      MPI_File fh_;
      MPI_Datatype durationType_;
    };

    // to read, in the manager process
    class DurationsReadFile{
    public:
      DurationsReadFile(const LevelVector& levelVector) : DurationsReadFile(getFilename(levelVector)) {}

      std::vector<durationInformation> readFromBeginning(size_t numberOfItems);
      
      ~DurationsReadFile(){
        MPI_File_close(&fh_); 
      }
      
    private:   
      DurationsReadFile(std::string filename){
        durationType_ = createMPIDurationType();
        MPI_File_open(MPI_COMM_SELF, filename.c_str(), 
          MPI_MODE_RDONLY, //MPI_MODE_CREATE|
          MPI_INFO_NULL, &fh_
        ); 
      }
      MPI_File fh_;
      MPI_Datatype durationType_;
    };
  }//namespace durationsFile

using namespace durationsFile;

class LearningLoadModel : public LoadModel {

 typedef std::chrono::high_resolution_clock::time_point::duration time_d;

 public:
  LearningLoadModel(std::vector<LevelVector> levelVectors){
    for(LevelVector& l : levelVectors){
      files_.emplace(std::make_pair(l, std::unique_ptr<DurationsReadFile>(new DurationsReadFile(l))));
    }
    numberOfEntriesExpected_ = 0;
  }

  inline real eval(const LevelVector& l);

  virtual ~LearningLoadModel() = default;
 
 private:
  size_t numberOfEntriesExpected_;
  std::map<LevelVector, std::unique_ptr<DurationsReadFile>> files_;
  // std::map<LevelVector, std::vector<long int>> durations_;
};

//using simple averaging for now //TODO
inline real LearningLoadModel::eval(const LevelVector& l) {
  real ret(0.0);
  std::vector<durationInformation> col = files_[l]->readFromBeginning(numberOfEntriesExpected_); //TODO make more efficient, store results in memory
  assert(col[0].duration != 0);
  // if no data yet, use linear load model
  // if(col.empty()){
  //   LinearLoadModel llm = LinearLoadModel();
  //   ret = llm.eval(l);
  // }
  // else{
    std::for_each(col.begin(), col.end(), [&] (durationInformation n) {
      ret += n.duration;
    });
    ret /= static_cast<double>(col.size());
  // }
  return ret;
}

} /* namespace combigrid */

#endif /* LEARNINGLOADMODEL_HPP_ */
