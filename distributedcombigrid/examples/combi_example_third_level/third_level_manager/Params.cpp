#include "Params.hpp"

Params::Params(const std::string& filename)
{
  readParameterFile(filename);
}

void Params::readParameterFile(const std::string& filename)
{
  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini(filename, pt);
  _brokerURL       = pt.get<std::string>("RabbitMQ.url");
  _dimension       = pt.get<u_int>("DistributedCombi.dimension");
  _numCombinations = pt.get<u_int>("DistributedCombi.numCombinations");

  std::string systemsCSV = pt.get<std::string>("DistributedCombi.systems");
  boost::split(_systemNames, systemsCSV, boost::is_any_of(","));
  _numSystems = _systemNames.size();

  // read list of common subspaces given by l11,l12 l21,l22 ... ln1,ln2 (for 2D)
  std::string commonSubspaces_str = pt.get<std::string>("DistributedCombi.commonSubspaces");
  std::vector<std::string> levels_str;
  boost::split(levels_str, commonSubspaces_str, boost::is_any_of(" "));
  for (size_t i = 0; i < levels_str.size(); i++)  {
    combigrid::LevelVector level(_dimension);
    std::istringstream iss(levels_str[i]);
    std::string token;
    u_int dim = 0;
    while ( !std::getline(iss, token, ',') && (dim < _dimension))
    {
      level[dim] = std::stoi(token);
      dim++;
    }
    _commonSubspaces.push_back(level);
  }
}

std::string Params::getbrokerURL() const
{
  return _brokerURL;
}

u_int Params::getNumSystems() const
{
  return _numSystems;
}

std::vector<combigrid::LevelVector> Params::getCommonLevels() const
{
  return _commonSubspaces;
}

std::vector<std::string> Params::getSystemNames() const
{
  return _systemNames;
}
