#include "Params.hpp"
#include <boost/program_options.hpp>

Params::Params()
{
}

void Params::loadFromFile(const std::string& filename)
{
  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini(filename, pt);
  port_           = pt.get<u_short>("Connection.port");
  numSystems_     = pt.get<u_int>("DistributedCombi.numSystems");
  chunksize_      = pt.get<size_t>("Connection.chunksize", 131072);
  loadedFromFile_ = true;
}

void Params::loadFromCmd(int argc, char* argv[])
{
  boost::program_options::options_description desc("Allowed options");
  desc.add_options()
          ("port", boost::program_options::value<u_short>()->required(), "port for connections")
          ("chunksize", boost::program_options::value<size_t>()->implicit_value(131072), "size of chunks to recv and send")
          ("numSystems", boost::program_options::value<u_int>()->implicit_value(2), "numer of systems that are going to connect");

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
  boost::program_options::notify(vm);
  port_           = vm["port"].as<u_short>();
  if(vm.count("numSystems")) {
    numSystems_     = vm["numSystems"].as<u_int>();
  } else {
    numSystems_ = 2;
  }
  if(vm.count("chunksize")) {
    chunksize_      = vm["chunksize"].as<size_t>();
  } else {
    chunksize_ = 131072;
  }
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

size_t Params::getChunksize() const { return chunksize_; }
