#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/algorithm/string.hpp>
#include <string>
#include "../../../src/sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "../../../src/sgpp/distributedcombigrid/third_level/NetworkUtils.hpp"

class Params
{
  private:
    std::string                         _brokerURL        = "localhost";
    u_int                               _numSystems       = 0;
    u_int                               _dimension        = 0;
    u_int                               _numCombinations  = 0;
    std::vector<std::string>            _systemNames;
    std::vector<combigrid::LevelVector> _commonSubspaces;


  public:
    Params() = delete;
    Params(const std::string& filename);

    void readParameterFile(const std::string& filename);

    std::string                         getbrokerURL()    const;
    u_int                               getNumSystems()   const;
    std::vector<combigrid::LevelVector> getCommonLevels() const;
    std::vector<std::string>            getSystemNames()  const;
};
