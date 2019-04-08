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
}

std::string Params::getbrokerURL() const
{
  return brokerURL_;
}

u_int Params::getNumSystems() const
{
  return numSystems_;
}

std::vector<std::string> Params::getSystemNames() const
{
  return systemNames_;
}


u_int Params::getDataPort() const
{
  return dataPort_;
}
