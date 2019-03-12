#include "Params.hpp"

Params::Params()
{
}

void Params::loadFromFile(const std::string& filename)
{
  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini(filename, pt);
  dataPort_        = pt.get<u_int>("Connection.dataPort");
  brokerURL_       = pt.get<std::string>("Connection.RabbitMQUrl");
  dimension_       = pt.get<u_int>("DistributedCombi.dimension");
  numCombinations_ = pt.get<u_int>("DistributedCombi.numCombinations");

  std::string systemsCSV = pt.get<std::string>("DistributedCombi.systems");
  boost::split(systemNames_, systemsCSV, boost::is_any_of(","));
  numSystems_ = systemNames_.size();

  // read list of common subspaces given for 2D by l11,l12 l21,l22 ... ln1,ln2
  std::string commonSubspaces_str = pt.get<std::string>("DistributedCombi.commonSubspaces");
  std::vector<std::string> levels_str;
  boost::split(levels_str, commonSubspaces_str, boost::is_any_of(" "));
  for (size_t i = 0; i < levels_str.size(); i++)  {
    combigrid::LevelVector level(dimension_);
    std::istringstream iss(levels_str[i]);
    std::string token;
    u_int dim = 0;
    while ( !std::getline(iss, token, ',') && (dim < dimension_))
    {
      level[dim] = std::stoi(token);
      dim++;
    }
    commonSubspaces_.push_back(level);
  }
}

std::string Params::getbrokerURL() const
{
  return brokerURL_;
}

u_int Params::getNumSystems() const
{
  return numSystems_;
}

std::vector<combigrid::LevelVector> Params::getCommonLevels() const
{
  return commonSubspaces_;
}

std::vector<std::string> Params::getSystemNames() const
{
  return systemNames_;
}


u_int Params::getDataPort() const
{
  return dataPort_;
}
