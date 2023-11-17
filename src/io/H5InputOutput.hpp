#pragma once
#include <numeric>

#include "utils/Types.hpp"

#ifdef DISCOTEC_USE_HIGHFIVE
// highfive is a C++ hdf5 wrapper, available in spack (-> configure with right boost and mpi
// versions)
#include <highfive/H5File.hpp>
#endif

namespace combigrid {
namespace h5io {

template <typename T>
void writeValuesToH5File(
    const T& values, const std::string& fileName, const std::string& groupName,
    const std::string& dataSetName,
    combigrid::real simulationTime = std::numeric_limits<combigrid::real>::quiet_NaN()) {
#ifdef DISCOTEC_USE_HIGHFIVE
  // check if file already exists, if no, create, if yes, overwrite
  HighFive::File h5_file(fileName, HighFive::File::OpenOrCreate | HighFive::File::ReadWrite |
                                       HighFive::File::Overwrite);

  HighFive::Group group;
  if (h5_file.exist(groupName)) {
    group = h5_file.getGroup(groupName);
  } else {
    group = h5_file.createGroup(groupName);
  }

  HighFive::DataSet dataset =
      group.createDataSet<CombiDataType>(dataSetName, HighFive::DataSpace::From(values));
  dataset.write(values);

  if (!std::isnan(simulationTime)) {
    HighFive::Attribute aTime = dataset.createAttribute<combigrid::real>(
        "simulation_time", HighFive::DataSpace::From(simulationTime));
    aTime.write(simulationTime);
  }
#else  // if not compiled with hdf5
  throw std::runtime_error("requesting hdf5 write but built without hdf5 support");
#endif
}

template <typename T>
void readValuesFromH5File(T& values, const std::string& fileName) {
#ifdef DISCOTEC_USE_HIGHFIVE
  HighFive::File h5_file(fileName, HighFive::File::ReadOnly);

  // assume one group in file
  HighFive::Group group = h5_file.getGroup(h5_file.listObjectNames()[0]);

  // assume one dataset in group
  auto dataset = group.getDataSet(group.listObjectNames()[0]);

  dataset.read(values);

#else  // if not compiled with hdf5
  throw std::runtime_error("requesting hdf5 read but built without hdf5 support");
#endif
}

// some instantiations
void readH5Coordinates(std::vector<std::vector<real>>& coordinates, const std::string& saveFilePath);

void readH5Values(std::vector<real>& values, const std::string& saveFilePath);

}  // namespace h5io
}  // namespace combigrid