/*
 * combi_example.cpp
 *
 *  Created on: Sep 23, 2015
 *      Author: heenemo
 */
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/serialization/export.hpp>
#include <string>
#include <vector>

// compulsory includes for basic functionality
#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"
#include "sgpp/distributedcombigrid/combischeme/CombiThirdLevelScheme.hpp"
#include "sgpp/distributedcombigrid/fault_tolerance/FaultCriterion.hpp"
#include "sgpp/distributedcombigrid/fault_tolerance/StaticFaults.hpp"
#include "sgpp/distributedcombigrid/fault_tolerance/WeibullFaults.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LinearLoadModel.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupManager.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupWorker.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessManager.hpp"
#include "sgpp/distributedcombigrid/task/Task.hpp"
#include "sgpp/distributedcombigrid/utils/MonteCarlo.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
// include user specific task. this is the interface to your application

// to allow using test tasks
#define BOOST_CHECK

#include "TaskAdvection.hpp"

using namespace combigrid;

// this is necessary for correct function of task serialization
BOOST_CLASS_EXPORT(TaskAdvection)
BOOST_CLASS_EXPORT(StaticFaults)
BOOST_CLASS_EXPORT(WeibullFaults)
BOOST_CLASS_EXPORT(FaultCriterion)


int main(int argc, char** argv) {
  
// read in parameter file
  std::string paramfile = "ctparam";
  if (argc > 1) paramfile = argv[1];
  boost::property_tree::ptree cfg;
  boost::property_tree::ini_parser::read_ini(paramfile, cfg);

  Stats::initialize();

    std::cout << "DisCoTec! reading " << paramfile << std::endl;
    DimType dim = cfg.get<DimType>("ct.dim");
    combigrid::real dt;
    size_t nsteps, ncombi;
    ncombi = cfg.get<size_t>("ct.ncombi");
    dt = cfg.get<combigrid::real>("application.dt");
    nsteps = cfg.get<size_t>("application.nsteps");


    std::vector<size_t> numValuesToTry{1000, 10000, 100000, 1000000, 10000000};
    std::vector<CombiDataType> l2Values;
    for (auto& numValues : numValuesToTry) {
      l2Values.clear();
      for (int i = 0; i < 20; ++i) {
        Stats::startEvent("manager calculate errors");
	auto interpolationCoords = montecarlo::getRandomCoordinates(numValues, dim);
        // calculate monte carlo errors
        TestFn initialFunction;
        real l0Reference = 0., l1Reference = 0., l2Reference = 0.;
        for (size_t i = 0; i < interpolationCoords.size(); ++i) {
          auto analyticalSln =
              initialFunction(interpolationCoords[i], static_cast<double>(ncombi * nsteps) * dt);
          l0Reference = std::max(analyticalSln, l0Reference);
          l1Reference += analyticalSln;
          l2Reference += std::pow(analyticalSln, 2);
        }
	l1Reference /= numValues;
        l2Reference /= numValues;

        Stats::stopEvent("manager calculate errors");
	std::cout << "Monte carlo values on " << numValues << " points are " << l0Reference << ", "
                  << l1Reference << ", and " << l2Reference << " in total." << std::endl;
	l2Values.push_back(l2Reference);
      }
	double sum = std::accumulate(l2Values.begin(), l2Values.end(), 0.0);
	double mean = sum / l2Values.size();

	double sq_sum = std::inner_product(l2Values.begin(), l2Values.end(), l2Values.begin(), 0.0);
	double variance = sq_sum / l2Values.size() - mean * mean;
	std::cout << "Monte carlo variance on " << numValues << " points is " << variance << std::endl;
    }
  Stats::finalize();
  return 0;
}
