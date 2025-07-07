/*
 * combi_example.cpp
 *
 *  Created on: Sep 23, 2015
 *      Author: heenemo
 */
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include <boost/serialization/export.hpp>
#include <string>
#include <vector>

// compulsory includes for basic functionality
#include "discotec/fault_tolerance/FaultCriterion.hpp"
#include "discotec/fault_tolerance/StaticFaults.hpp"
#include "discotec/fault_tolerance/WeibullFaults.hpp"
#include "discotec/fullgrid/FullGrid.hpp"
#include "discotec/combischeme/CombiMinMaxScheme.hpp"
#include "discotec/io/BroadcastParameters.hpp"
#include "discotec/loadmodel/LinearLoadModel.hpp"
#include "discotec/manager/CombiParameters.hpp"
#include "discotec/manager/ProcessGroupManager.hpp"
#include "discotec/manager/ProcessGroupWorker.hpp"
#include "discotec/manager/ProcessManager.hpp"
#include "discotec/Task.hpp"
#include "discotec/utils/Types.hpp"
// include user specific task. this is the interface to your application
#include "TaskExample.hpp"

using namespace combigrid;

// this is necessary for correct function of task serialization
#include "discotec/utils/BoostExports.hpp"
BOOST_CLASS_EXPORT(TaskExample)

int main(int argc, char** argv) {
  [[maybe_unused]] auto mpiOnOff = MpiOnOff(&argc, &argv);

  /* when using timers (TIMING is defined in Stats), the Stats class must be
   * initialized at the beginning of the program. (and finalized in the end)
   */
  Stats::initialize();

  // only one rank reads inputs and broadcasts to others
  std::string paramfile = "ctparam";
  if (argc > 1) paramfile = argv[1];
  boost::property_tree::ptree cfg =
      broadcastParameters::getParametersFromRankZero(paramfile, MPI_COMM_WORLD);

  // number of process groups and number of processes per group
  size_t ngroup = cfg.get<size_t>("manager.ngroup");
  size_t nprocs = cfg.get<size_t>("manager.nprocs");

  // divide the MPI processes into process group and initialize the
  // corresponding communicators
  theMPISystem()->init(ngroup, nprocs);

  // this code is only executed by the manager process
  WORLD_MANAGER_EXCLUSIVE_SECTION {
    /* create an abstraction of the process groups for the manager's view
     * a pgroup is identified by the ID in gcomm
     */
    ProcessGroupManagerContainer pgroups;
    for (size_t i = 0; i < ngroup; ++i) {
      int pgroupRootID(i);
      pgroups.emplace_back(std::make_shared<ProcessGroupManager>(pgroupRootID));
    }

    // create load model
    std::unique_ptr<LoadModel> loadmodel = std::unique_ptr<LoadModel>(new LinearLoadModel());


    /* read in parameters from ctparam */
    DimType dim = cfg.get<DimType>("ct.dim");
    LevelVector lmin(dim), lmax(dim), leval(dim);
    std::vector<int> p(dim);
    combigrid::real dt;
    size_t nsteps, ncombi;
    cfg.get<std::string>("ct.lmin") >> lmin;
    cfg.get<std::string>("ct.lmax") >> lmax;
    cfg.get<std::string>("ct.leval") >> leval;
    cfg.get<std::string>("ct.p") >> p;
    ncombi = cfg.get<size_t>("ct.ncombi");
    dt = cfg.get<combigrid::real>("application.dt");
    nsteps = cfg.get<size_t>("application.nsteps");

    // todo: read in boundary vector from ctparam
    std::vector<BoundaryType> boundary(dim, 2);

    // check whether parallelization vector p agrees with nprocs
    IndexType checkProcs = 1;
    for (auto k : p) checkProcs *= k;
    assert(checkProcs == IndexType(nprocs));

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
    std::vector<size_t> taskIDs;
    for (size_t i = 0; i < levels.size(); i++) {
      Task* t = new TaskExample(levels[i], boundary, coeffs[i], loadmodel.get(), dt, nsteps, p);
      tasks.push_back(t);
      taskIDs.push_back(t->getID());
    }

    // create combiparameters
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs, ncombi, 1,
                           CombinationVariant::sparseGridReduce, p);

    // create abstraction for Manager
    ProcessManager manager(pgroups, tasks, params, std::move(loadmodel));

    std::cout << "set up component grids and run until first combination point" << std::endl;

    /* distribute task according to load model and start computation for
     * the first time */
    Stats::startEvent("manager run first");
    manager.runfirst();
    Stats::stopEvent("manager run first");

    for (size_t i = 0; i < ncombi; ++i) {

      // start = MPI_Wtime();

      Stats::startEvent("combine");
      manager.combine();
      Stats::stopEvent("combine");
      // finish = MPI_Wtime();
      // std::cout << "combination " << i << " took: " << finish-start << " seconds" << std::endl;

      // evaluate solution and
      // write solution to file
      std::string filename("out/solution_" + std::to_string(i) + ".out");
      Stats::startEvent("manager write solution");
      manager.parallelEval(leval, filename, 0);
      Stats::stopEvent("manager write solution");

      std::cout << "run until combination point " << i + 1 << std::endl;

      // run tasks for next time interval
      // start = MPI_Wtime();
      Stats::startEvent("manager run");
      manager.runnext();
      Stats::stopEvent("manager run");
      // finish = MPI_Wtime();
      // std::cout << "calculation " << i << " took: " << finish-start << " seconds" << std::endl;
    }

    Stats::startEvent("manager write vtk");
    for (auto group : pgroups)
      manager.writeCombigridsToVTKPlotFile(group);
    Stats::stopEvent("manager write vtk");

    // send exit signal to workers in order to enable a clean program termination
    manager.exit();
  }

  // this code is only execute by the worker processes
  else {
    // create abstraction of the process group from the worker's view
    ProcessGroupWorker pgroup;

    // wait for instructions from manager
    SignalType signal = -1;

    while (signal != EXIT) signal = pgroup.wait();
  }

  Stats::finalize();

  /* write stats to json file for postprocessing */
  Stats::write("timers.json");

  return 0;
}
