#include "io/H5InputOutput.hpp"

namespace combigrid {
namespace h5io {

// some instantiations
void readH5Coordinates(std::vector<std::vector<real>>& coordinates, std::string saveFilePath) {
  return h5io::readValuesFromH5File(coordinates, saveFilePath);
}

void readH5Values(std::vector<real>& values, std::string saveFilePath) {
  return h5io::readValuesFromH5File(values, saveFilePath);
}

}  // namespace h5io
}  // namespace combigrid