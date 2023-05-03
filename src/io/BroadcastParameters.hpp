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

boost::property_tree::ptree getParametersFromRankZero(const std::string& parameterFileName,
                                                      MPI_Comm comm) {
  // only one rank reads inputs and broadcasts to others
  auto mpiRank = getCommRank(comm);
  std::string parameterString;
  int parameterStringSize = 0;
  if (mpiRank == 0) {
    // read in parameter file
    std::ifstream ifs(parameterFileName);
    parameterString.assign((std::istreambuf_iterator<char>(ifs)),
                           (std::istreambuf_iterator<char>()));
    parameterStringSize = parameterString.size();
  }
  // avoid: copies of strings in this function
  // MPIUtils::broadcastClass<std::string>(parameterString, 0, comm);

  // broadcast string size
  MPI_Bcast(&parameterStringSize, 1, MPI_INT, 0, comm);
  parameterString.resize(parameterStringSize);

  // broadcast string
  MPI_Bcast(&parameterString[0], parameterString.size(), MPI_CHAR, 0, comm);

  boost::iostreams::stream parameterStream(
      boost::iostreams::array_source(&parameterString[0], parameterString.size()));
  boost::property_tree::ptree cfg;
  boost::property_tree::ini_parser::read_ini(parameterStream, cfg);
  return cfg;
}

}  // namespace broadcastParameters
}  // namespace combigrid