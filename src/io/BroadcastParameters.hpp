#pragma once

#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "mpi/MPISystem.hpp"
#include "mpi/MPIUtils.hpp"

namespace combigrid {
namespace broadcastParameters {

inline std::string getFileContentsFromRankZero(const std::string& fileName, MPI_Comm comm) {
  // only one rank reads inputs and broadcasts to others
  auto mpiRank = getCommRank(comm);
  std::string contentString;
  int contentStringSize = 0;
  if (mpiRank == 0) {
    // read in content file
    std::ifstream ifs(fileName);
    contentString.assign((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
    contentStringSize = contentString.size();
  }
  // avoid: copies of strings in this function
  // MPIUtils::broadcastClass<std::string>(contentString, 0, comm);

  // broadcast string size
  MPI_Bcast(&contentStringSize, 1, MPI_INT, 0, comm);
  contentString.resize(contentStringSize);

  // broadcast string
  MPI_Bcast(&contentString[0], contentString.size(), MPI_CHAR, 0, comm);
  return contentString;
}

inline boost::property_tree::ptree getParametersFromRankZero(const std::string& parameterFileName,
                                                             MPI_Comm comm) {
  std::string parameterString = getFileContentsFromRankZero(parameterFileName, comm);
  boost::iostreams::stream parameterStream(
      boost::iostreams::array_source(&parameterString[0], parameterString.size()));

  boost::property_tree::ptree cfg;
  boost::property_tree::ini_parser::read_ini(parameterStream, cfg);
  return cfg;
}

inline boost::property_tree::ptree getJsonFromRankZero(const std::string& jsonFileName,
                                                       MPI_Comm comm) {
  // beware: the json files can get some MB large
  std::string jsonString = getFileContentsFromRankZero(jsonFileName, comm);
  boost::iostreams::stream jsonStream(
      boost::iostreams::array_source(&jsonString[0], jsonString.size()));
  boost::property_tree::ptree cfg;
  boost::property_tree::json_parser::read_json(jsonStream, cfg);
  return cfg;
}

}  // namespace broadcastParameters
}  // namespace combigrid