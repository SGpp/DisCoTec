/*
 * LearningLoadModel.hpp
 *
 *  Created on: Oct 9, 2013
 *      Author: heenemo
 */

#ifndef LEARNINGLOADMODEL_HPP_
#define LEARNINGLOADMODEL_HPP_

#include <iostream>
#include <string>
#include <chrono>
#include <memory>

#include "sgpp/distributedcombigrid/loadmodel/LoadModel.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LinearLoadModel.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/utils/Stats.hpp"
#include "CSVFile.hpp"

#ifdef USE_HDF5
// cf. http://micro.stanford.edu/wiki/Install_HDF5:
// "Note that, as of HDF5 version 1.9 there is no way to use the MPI version and the 
// C++ interfaces together. The C interface can be used from C++ anyway."
// #include "H5Cpp.h"
#warning "HDF5 data write: not actually implemented!"
#include "H5public.h"
#endif

namespace combigrid {

class LearningLoadModel : public LoadModel {
 typedef std::chrono::high_resolution_clock::time_point::duration time_d;
 public:
  LearningLoadModel(){
  };

  inline real eval(const LevelVector& l);

  virtual ~LearningLoadModel() = default;

  std::map<LevelVector, std::unique_ptr<csvfile>> files_;

  void addDataPoint(const LevelVector& l, const Stats::Event event, size_t nProcesses){ //TODO include "metadata" in model: nrg.dat, parameters etc.
    createExtensibleData(l);
    long int duration = (event.end - event.start).count(); 
    // long int order = event.end.time_since_epoch().count();

    durationInformation info = {duration, nProcesses};

    #ifdef USE_HDF5
        dataset->write(info, mtype_);
    #else //def USE_HDF5
        (*files_[l]) << (info);
    #endif //def USE_HDF5
  }

  struct durationInformation{
        long int duration;
        size_t nProcesses;
        // long int order;
        friend csvfile& operator <<(csvfile& cf, durationInformation const& d)
        {
            return cf << (d.duration) << d.nProcesses << endrow;
        }
  };

private:

#ifdef USE_HDF5
  // creating an array of structures with this structure
  // cf https://gist.github.com/YukiSakamoto/6319458
  std::unique_ptr<H5::H5File> file_;
  H5::CompType mtype_(sizeof(durationInformation));
  H5::DataSpace space_;
  std::unique_ptr<H5::DataSet> dataset_;

  void getOrCreateH5File(std::string file_name){
      if( std::ifstream(file_name)){
        file_(H5Fcreate( file_name, H5F_ACC_RDWR, H5P_DEFAULT, H5P_DEFAULT ));
      }
      else{
        file_(H5Fcreate( file_name, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT ));
      }
  }

  void prepareMtype(){
    // defining the datatype to pass HDF55
    mtype_.insertMember("level vector", HOFFSET(durationInformation, l), H5::PredType::NATIVE_INT); //TODO
    mtype_.insertMember("duration", HOFFSET(durationInformation, sex), H5::PredType::C_S1);
    mtype_.insertMember("numberOfProcesses", HOFFSET(durationInformation, height), H5::PredType::NATIVE_INT);
    // mtype_.insertMember("order", HOFFSET(durationInformation, order), H5::PredType::NATIVE_FLOAT);
  }

  void addDataSetToFile(std::string datasetName = std::to_string(time(0))){
    prepareMtype();
    // preparation of a dataset

    int rank = 1;
    hsize_t dims[rank] = {0, 0};
    hsize_t max_dims[rank] = {H5S_UNLIMITED, H5S_UNLIMITED};
    space_ = H5Screate_simple(ndims, dims, max_dims);

    // cf. https://stackoverflow.com/questions/15379399/writing-appending-arrays-of-float-to-the-only-dataset-in-hdf5-file-in-c
    // for more efficient extensibility of data set
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(plist, H5D_CHUNKED);
    hsize_t chunk_dims[rank] = {1, 100};
    H5Pset_chunk(plist, ndims, chunk_dims);

    dataset(new H5::DataSet(file->createDataSet(datasetName, mtype, space)));
  }
#else //def USE_HDF5

  void createExtensibleData(const LevelVector& l){
    if (files_.find(l)==files_.end())
      files_.emplace(std::make_pair(l, std::unique_ptr<csvfile>(new csvfile("./learnloadmodel" + std::to_string(time(0)) + "_" + toString(l) + ".csv"))));
  }
#endif //def USE_HDF5
};



//using simple averaging for now //TODO
//inline real LearningLoadModel::eval(const LevelVector& l) const {
inline real LearningLoadModel::eval(const LevelVector& l) {
  real ret(0.0);
  std::vector<long int> col = files_[l]->readColumn<long int>(0); //TODO make more efficient
  // if no data yet, use linear load model
  if(col.empty()){
    LinearLoadModel llm = LinearLoadModel();
    ret = llm.eval(l);
  }
  else{
    std::for_each(col.begin(), col.end(), [&] (int n) {
        ret += n;
    });
    ret /= static_cast<double>(col.size());
  }
  return ret;
}

} /* namespace combigrid */

#endif /* LEARNINGLOADMODEL_HPP_ */
