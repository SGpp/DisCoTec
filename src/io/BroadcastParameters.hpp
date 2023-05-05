#pragma once

#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <string>
#include <vector>

#include "mpi/MPISystem.hpp"

namespace combigrid {
namespace broadcastParameters {

std::string getFileContentsFromRankZero(const std::string& fileName, const MPI_Comm& comm);

boost::property_tree::ptree getParametersFromRankZero(const std::string& parameterFileName,
                                                      const MPI_Comm& comm);

boost::property_tree::ptree getJsonFromRankZero(const std::string& jsonFileName, MPI_Comm comm);

std::vector<std::vector<real>> getCoordinatesFromRankZero(const std::string& h5FileName,
                                                          const MPI_Comm& comm);

}  // namespace broadcastParameters
}  // namespace combigrid