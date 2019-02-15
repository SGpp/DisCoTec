/*
 * combi_example.cpp
 *
 *  Created on: Sep 23, 2015
 *      Author: heenemo
 */
#include <mpi.h>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/serialization/export.hpp>
#include <string>
#include <vector>
#include <algorithm>

// compulsory includes for basic functionality
#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LinearLoadModel.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupManager.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupWorker.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessManager.hpp"
#include "sgpp/distributedcombigrid/task/Task.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
// include user specific task. this is the interface to your application
#include "TaskExample.hpp"

using namespace combigrid;

// this is necessary for correct function of task serialization
BOOST_CLASS_EXPORT(TaskExample)

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  /* when using timers (TIMING is defined in Stats), the Stats class must be
   * initialized at the beginning of the program. (and finalized in the end)
   */
  Stats::initialize();

  // read in parameter file
  boost::property_tree::ptree cfg;
  if (argc > 1) {
    boost::property_tree::ini_parser::read_ini(argv[1], cfg);
  } else {
    std::cout << "Usage:\n ./combi_example parameterFile";
    return 0;
  }

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
    LoadModel* loadmodel = new LinearLoadModel();

    /* read in parameters from ctparam */
    DimType dim = cfg.get<DimType>("ct.dim");
    LevelVector lmin(dim), lmax(dim), leval(dim);
    IndexVector p(dim);
    combigrid::real dt;
    size_t nsteps, ncombi;
    cfg.get<std::string>("ct.lmin") >> lmin;
    cfg.get<std::string>("ct.lmax") >> lmax;
    cfg.get<std::string>("ct.leval") >> leval;
    cfg.get<std::string>("ct.p") >> p;
    ncombi = cfg.get<size_t>("ct.ncombi");
    dt = cfg.get<combigrid::real>("application.dt");
    nsteps = cfg.get<size_t>("application.nsteps");

    // read in third Level specific parameters
    std::string thirdLevelHost = cfg.get<std::string>("thirdLevel.host");
    std::string systemName = cfg.get<std::string>("thirdLevel.systemName");
    unsigned short thirdLevelDataPort = cfg.get<unsigned short>("thirdLevel.dataPort");

    // todo: read in boundary vector from ctparam
    std::vector<bool> boundary(dim, true);

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
    std::vector<LevelVector> fullScheme = combischeme.getCombiSpaces();
    std::vector<combigrid::real> fullSchemeCoeffs = combischeme.getCoeffs();

    // For example purpose we just split the levels in half, and assign each
    // half to a system
    auto mid = fullScheme.begin() + ((fullScheme.size()-1) / 2);
    auto midC = fullSchemeCoeffs.begin() + ((fullSchemeCoeffs.size()-1) / 2);
    std::vector<LevelVector> lowerHalf(fullScheme.begin(), mid);
    std::vector<combigrid::real> lowerCoeffs(fullSchemeCoeffs.begin(), midC);
    std::vector<LevelVector> upperHalf(mid+1, fullScheme.end());
    std::vector<combigrid::real> upperCoeffs(midC+1, fullSchemeCoeffs.end());
    assert( !lowerHalf.empty() && !upperHalf.empty() );

    // compute common subspaces:
    // therefore we compute the component wise maximum level which is contained
    // in both sets
    LevelVector maxLevel(dim);
    for (size_t i = 0; i < dim; i++) {
      int lowerMax = 0;
      int upperMax = 0;
      for (size_t j = 0; j < lowerHalf.size(); j++) {
        if (lowerMax < lowerHalf[j][i])
          lowerMax = lowerHalf[j][i];
      }
      for (size_t j = 0; j < upperHalf.size(); j++) {
        if (upperMax < upperHalf[j][i])
          upperMax = upperHalf[j][i];
      }
      maxLevel[i] = std::min(lowerMax, upperMax);
    }

    // by creating a dummy sparse grid with this level, we can extract the subspaces
    SGrid<real> sg(dim, maxLevel, maxLevel, boundary);

    std::vector<LevelVector> commonSubspaces;
    for (size_t subspaceID = 0; subspaceID < sg.getSize(); ++subspaceID) {
      const LevelVector& subL = sg.getLevelVector(subspaceID);
      commonSubspaces.push_back(subL);
    }

    // output combination scheme
    //std::cout << "lmin = " << lmin << std::endl;
    //std::cout << "lmax = " << lmax << std::endl;
    //std::cout << "CombiScheme: " << std::endl;
    //std::cout << combischeme << std::endl;

    // output combi scheme
    std::cout << "UpperHalf:" << std::endl;
    for (int i = 0; i < upperHalf.size(); i++)
      std::cout << upperHalf[i] << ", " << std::endl;
    std::cout << std::endl;
    std::cout << "LowerHalf:" << std::endl;
    for (int i = 0; i < lowerHalf.size(); i++)
      std::cout << lowerHalf[i] << ", " << std::endl;
    std::cout << std::endl;

    // output common subspaces
    std::cout << "Num common subspaces: " << commonSubspaces.size() << std::endl;
    std::cout << "CommonSubspaces:" << std::endl;
    for (int i = 0; i < commonSubspaces.size(); i++)
      std::cout << commonSubspaces[i] << ", " << std::endl;

    // create Tasks
    TaskContainer tasks;
    std::vector<int> taskIDs;
    CombiParameters params;

    std::vector<LevelVector>& levels = fullScheme;
    std::vector<combigrid::real>& coeffs = fullSchemeCoeffs;
    if (systemName == "neon") {
      levels = upperHalf;
      coeffs = upperCoeffs;
    }
    if (systemName == "helium") {
      levels = lowerHalf;
      coeffs = lowerCoeffs;
    }

    for (size_t i = 0; i < upperHalf.size(); i++) {
      Task* t = new TaskExample(dim, upperHalf[i], boundary, upperCoeffs[i], loadmodel, dt, nsteps, p);
      tasks.push_back(t);
      taskIDs.push_back(t->getID());
    }

    assert(!tasks.empty());

    // create combiparameters
    params = CombiParameters(dim, lmin, lmax, boundary, upperHalf, upperCoeffs, taskIDs, thirdLevelHost, thirdLevelDataPort, systemName, commonSubspaces);

    // create abstraction for Manager
    ProcessManager manager(pgroups, tasks, params);

    // the combiparameters are sent to all process groups before the
    // computations start
    manager.updateCombiParameters();

    std::cout << "set up component grids and run until first combination point" << std::endl;

    /* distribute task according to load model and start computation for
     * the first time */
    Stats::startEvent("manager run first");
    manager.runfirst();
    Stats::stopEvent("manager run first");

    for (size_t i = 0; i < ncombi; ++i) {
      Stats::startEvent("combineThirdLevel");
      manager.combineThirdLevel();
      Stats::stopEvent("combineThirdLevel");

      // evaluate solution and
      // write solution to file
      std::string filename("out/solution_" + std::to_string(ncombi) + ".dat");
      Stats::startEvent("manager write solution");
      manager.parallelEval(leval, filename, 0);
      Stats::stopEvent("manager write solution");

      std::cout << "run until combination point " << i + 1 << std::endl;

      // run tasks for next time interval
      Stats::startEvent("manager run");
      manager.runnext();
      Stats::stopEvent("manager run");
    }

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

  MPI_Finalize();

  return 0;
}
