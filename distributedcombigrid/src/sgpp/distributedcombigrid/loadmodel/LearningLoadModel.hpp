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
    createExtensibleData();
  };

  inline real eval(const LevelVector& l) const;

  virtual ~LearningLoadModel() = default;


  void addDataPoint(const LevelVector& l, const Stats::Event event, size_t nProcesses);
  struct durationInformation{
        const LevelVector l;
        long int duration;
        size_t nProcesses;
        long int order;
        friend csvfile& operator <<(csvfile& cf, durationInformation const& d)
        {
            return cf << d.l << (d.duration) << d.nProcesses << d.order << endrow;
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
    mtype_.insertMember("order", HOFFSET(durationInformation, order), H5::PredType::NATIVE_FLOAT);
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

  std::unique_ptr<csvfile> file_;

  void createExtensibleData(){
    file_ = std::unique_ptr<csvfile>(new csvfile("learnloadmodel" + std::to_string(time(0)) + ".csv"));
    (*file_) << "level vector" << "duration" << "nProcesses" << "order" << endrow;
  }
#endif //def USE_HDF5


};

void LearningLoadModel::addDataPoint(const LevelVector& l, const Stats::Event event, size_t nProcesses){ //TODO "metadata": nrg.dat etc.
    long int duration = (event.end - event.start).count(); 
    long int order = event.end.time_since_epoch().count();

    durationInformation info = {l, duration, nProcesses, order};

    #ifdef USE_HDF5
        dataset->write(info, mtype_);
    #else //def USE_HDF5
        (*file_) << (info);
    #endif //def USE_HDF5
}

//using linear load model's eval function for now //TODO
//inline real LearningLoadModel::eval(const LevelVector& l) const {
inline real LearningLoadModel::eval(const LevelVector& l) const {
  real ret(1.0);

  for (size_t i = 0; i < l.size(); ++i) {
    ret *= std::pow(real(2.0), static_cast<real>(l[i]));
  }

  return ret;
}

} /* namespace combigrid */

#endif /* LEARNINGLOADMODEL_HPP_ */
