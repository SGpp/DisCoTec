#include "io/BroadcastParameters.hpp"

#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "io/H5InputOutput.hpp"
#include "mpi/MPISystem.hpp"
#include "mpi/MPIUtils.hpp"
#include "utils/MonteCarlo.hpp"

namespace combigrid {
namespace broadcastParameters {

std::string getFileContentsFromRankZero(const std::string& fileName, const MPI_Comm& comm) {
  // only one rank reads inputs and broadcasts to others
  auto mpiRank = getCommRank(comm);
  std::string contentString;
  int contentStringSize = 0;
  if (mpiRank == 0) {
    std::ifstream ifs(fileName);
    contentString.assign((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
    contentStringSize = contentString.size();
  }
  MPIUtils::broadcastContainer(contentString, 0, comm);
  return contentString;
}

boost::property_tree::ptree getParametersFromRankZero(const std::string& parameterFileName,
                                                      const MPI_Comm& comm) {
  std::string parameterString = getFileContentsFromRankZero(parameterFileName, comm);
  boost::iostreams::stream parameterStream(
      boost::iostreams::array_source(&parameterString[0], parameterString.size()));

  boost::property_tree::ptree cfg;
  boost::property_tree::ini_parser::read_ini(parameterStream, cfg);
  return cfg;
}

boost::property_tree::ptree getJsonFromRankZero(const std::string& jsonFileName, MPI_Comm comm) {
  // beware: the json files can get some MB large
  std::string jsonString = getFileContentsFromRankZero(jsonFileName, comm);
  boost::iostreams::stream jsonStream(
      boost::iostreams::array_source(&jsonString[0], jsonString.size()));
  boost::property_tree::ptree cfg;
  boost::property_tree::json_parser::read_json(jsonStream, cfg);
  return cfg;
}

std::vector<std::vector<real>> getCoordinatesFromRankZero(const std::string& h5FileName,
                                                          const MPI_Comm& comm) {
  std::vector<std::vector<real>> coordinates;
  std::vector<real> coordinatesSerialized;
  int coordinateSizes[2] = {0, 0};
  auto mpiRank = getCommRank(comm);
  if (mpiRank == 0) {
    h5io::readH5Coordinates(coordinates, h5FileName);
    coordinateSizes[0] = coordinates.size();
    assert(coordinateSizes[0] > 0);
    coordinateSizes[1] = coordinates[0].size();
    assert(coordinateSizes[1] > 0);
    coordinatesSerialized = combigrid::serializeInterpolationCoords(coordinates);
    assert(coordinatesSerialized.size() == coordinateSizes[0] * coordinateSizes[1]);
  }

  // broadcast vector sizes
  MPI_Bcast(&coordinateSizes, 2, MPI_INT, 0, comm);
  coordinatesSerialized.resize(coordinateSizes[0] * coordinateSizes[1]);

  // broadcast serialized vector
  MPI_Bcast(&coordinatesSerialized[0], coordinatesSerialized.size(),
            abstraction::getMPIDatatype(abstraction::getabstractionDataType<real>()), 0, comm);

  // split vector into coordinates
  coordinates = combigrid::deserializeInterpolationCoords(coordinatesSerialized,
                                                          static_cast<DimType>(coordinateSizes[1]));
  return coordinates;
}

}  // namespace broadcastParameters
}  // namespace combigrid