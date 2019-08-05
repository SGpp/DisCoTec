#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/algorithm/string.hpp>
#include <string>
#include "sgpp/distributedcombigrid/third_level/NetworkUtils.hpp"

class Params
{
  private:
    std::string                         brokerURL_       = "localhost";
    u_int                               dataPort_        = 0;
    u_int                               numSystems_      = 0;
    u_int                               dimension_       = 0;
    u_int                               numCombinations_ = 0;
    std::vector<std::string>            systemNames_;

  public:
    Params();

    void loadFromFile(const std::string& filename);

    std::string                         getbrokerURL()    const;
    u_int                               getDataPort()     const;
    u_int                               getNumSystems()   const;
    std::vector<std::string>            getSystemNames()  const;
};
