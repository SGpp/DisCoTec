#include "Params.hpp"

Params::Params()
{
}

void Params::loadFromFile(const std::string& filename)
{
  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini(filename, pt);
  port_           = pt.get<u_short>("Connection.port");
  numSystems_     = pt.get<u_int>("DistributedCombi.numSystems");
  loadedFromFile_ = true;
}

u_int Params::getNumSystems() const
{
  return numSystems_;
}

u_short Params::getPort() const
{
  return port_;
}

bool Params::areLoaded() const
{
  return loadedFromFile_;
}
