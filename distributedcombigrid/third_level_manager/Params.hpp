#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/algorithm/string.hpp>
#include <string>
#include "sgpp/distributedcombigrid/third_level/NetworkUtils.hpp"

class Params
{
  private:
    u_short port_           = 0;
    u_int   numSystems_     = 0;
    size_t  chunksize_      = 131072;
    bool    loadedFromFile_ = false;

   public:
    Params();

    void loadFromFile(const std::string& filename);
    void loadFromCmd(int argc, char* argv[]);

    bool areLoaded() const;
    u_short getPort() const;
    u_int getNumSystems() const;
    size_t getChunksize() const;
};
