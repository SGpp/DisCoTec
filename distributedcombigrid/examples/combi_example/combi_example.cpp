/*
 * combi_example.cpp
 *
 *  Created on: Sep 23, 2015
 *      Author: heenemo
 */
#include <mpi.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <numeric>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/serialization/export.hpp>

// compulsory includes for basic functionality
#include "sgpp/distributedcombigrid/task/Task.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LinearLoadModel.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupManager.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupWorker.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessManager.hpp"

// include user specific task. this is the interface to your application
#include "TaskExample.hpp"

using namespace combigrid;
using boost::property_tree::ptree;

// this is necessary for correct function of task serialization
BOOST_CLASS_EXPORT(TaskExample)

template <typename T>
std::vector<T> get_as_vector(ptree const& pt, ptree::key_type const& key)
{
  std::vector<T> result;
  result.reserve( pt.size() );

  for(const auto& child : pt.get_child(key)) {
    const ptree& childPt = child.second;
    result.push_back( childPt.get_value<T>() );
  }
  return result;
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  /* when using timers (TIMING is defined in Stats), the Stats class must be
   * initialized at the beginning of the program. (and finalized in the end)
   */
  Stats::initialize();

  // read in parameter file
  boost::property_tree::ptree cfg;
  boost::property_tree::json_parser::read_json("ctparam", cfg);

  // number of process groups and number of processes per group
  size_t ngroup = cfg.get<size_t>("manager.ngroup");
  std::vector<size_t> nprocs = get_as_vector<size_t>( cfg, "manager.nprocs" );

  // divide the MPI processes into process group and initialize the
  // corresponding communicators
  theMPISystem()->init( ngroup, nprocs );

  // this code is only executed by the manager process
  WORLD_MANAGER_EXCLUSIVE_SECTION {
    /* create an abstraction of the process groups for the manager's view
     * a pgroup is identified by the ID in gcomm
     */
    ProcessGroupManagerContainer pgroups;
    for (size_t i = 0; i < ngroup; ++i) {
      int pgroupRootID(i);
      pgroups.emplace_back(
          std::make_shared< ProcessGroupManager > ( pgroupRootID )
                          );
    }

    // create load model
    LoadModel* loadmodel = new LinearLoadModel();

    /* read in parameters from ctparam */
    DimType dim = cfg.get<DimType>("ct.dim");
    LevelVector lmin(dim), lmax(dim), leval(dim);
    std::unordered_map<size_t, IndexVector> pByNProcs;
    combigrid::real dt;
    size_t nsteps, ncombi;
    cfg.get<std::string>("ct.lmin") >> lmin;
    cfg.get<std::string>("ct.lmax") >> lmax;
    cfg.get<std::string>("ct.leval") >> leval;
    {
      const ptree& ps = cfg.get_child("ct.p");
      for(const auto& pConfig : ps) {
        IndexVector p( dim );
        pConfig.second.get_value<std::string>() >> p;
        IndexType procCount = std::accumulate( p.begin(), p.end(), 1, std::multiplies<IndexType>{} );
        pByNProcs[procCount] = std::move( p );
      }
    }
    ncombi = cfg.get<size_t>("ct.ncombi");
    dt = cfg.get<combigrid::real>("application.dt");
    nsteps = cfg.get<size_t>("application.nsteps");

    // todo: read in boundary vector from ctparam
    std::vector<bool> boundary(dim, true);

    // check whether parallelization vector p agrees with nprocs
    for(const auto& nproc : nprocs) {
      assert(pByNProcs.find(nproc) != pByNProcs.end() && "Need parallelization for every proc size");
    }

    /* generate a list of levelvectors and coefficients
     * CombiMinMaxScheme will create a classical combination scheme.
     * however, you could also read in a list of levelvectors and coefficients
     * from a file */
    CombiMinMaxScheme combischeme(dim, lmin, lmax);
    combischeme.createAdaptiveCombischeme();
    std::vector<LevelVector> levels = combischeme.getCombiSpaces();
    std::vector<combigrid::real> coeffs = combischeme.getCoeffs();

    // output combination scheme
    std::cout << "lmin = " << lmin << std::endl;
    std::cout << "lmax = " << lmax << std::endl;
    std::cout << "CombiScheme: " << std::endl;
    std::cout << combischeme << std::endl;

    // create Tasks
    TaskContainer tasks;
    std::vector<int> taskIDs;
    for (size_t i = 0; i < levels.size(); i++) {
      Task* t = new TaskExample(dim, levels[i], boundary, coeffs[i],
                                loadmodel, dt, nsteps, pByNProcs);
      tasks.push_back(t);
      taskIDs.push_back( t->getID() );
    }

    // create combiparameters
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs );

    // create abstraction for Manager
    ProcessManager manager(pgroups, tasks, params);

    // the combiparameters are sent to all process groups before the
    // computations start
    manager.updateCombiParameters();

    std::cout << "set up component grids and run until first combination point"
              << std::endl;
    /* distribute task according to load model and start computation for
     * the first time */

    Stats::startEvent("manager run first");
    manager.runfirst();
    Stats::stopEvent("manager run first");

    std::ofstream myfile;
    myfile.open("out/solution.dat");

    for (size_t i = 0; i < ncombi; ++i) {
      Stats::startEvent("combine");
      manager.combine();
      Stats::stopEvent("combine");

      // evaluate solution
      FullGrid<CombiDataType> fg_eval(dim, leval, boundary);
      Stats::startEvent("eval");
      manager.gridEval(fg_eval);
      Stats::stopEvent("eval");

      // write solution to file
      Stats::startEvent("manager write solution");
      std::vector<double> coords(dim, 0.0);
      for (int i = 0; i < fg_eval.getNrElements(); i++) {
        if (i % fg_eval.length(0) == 0 && i > 0) {
          myfile << std::endl;
        }

        fg_eval.getCoords(i, coords);
        myfile << coords[0] << "\t" << coords[1] << "\t"
               << fg_eval.getElementVector()[i] << std::endl;
      }
      myfile << std::endl << std::endl;
      Stats::stopEvent("manager write solution");

      std::cout << "run until combination point " << i+1 << std::endl;

      // run tasks for next time interval
      Stats::startEvent("manager run");
      manager.runnext();
      Stats::stopEvent("manager run");
    }

    myfile.close();

    // send exit signal to workers in order to enable a clean program termination
    manager.exit();
  }

  // this code is only execute by the worker processes
  else {
    // create abstraction of the process group from the worker's view
    ProcessGroupWorker pgroup;

    // wait for instructions from manager
    SignalType signal = -1;

    while (signal != EXIT)
      signal = pgroup.wait();
  }

  Stats::finalize();

  /* write stats to json file for postprocessing */
  Stats::write( "timers.json" );

  MPI_Finalize();

  return 0;
}
