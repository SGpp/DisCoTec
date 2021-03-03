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
#include "sgpp/distributedcombigrid/utils/Types.hpp"
// include user specific task. this is the interface to your application

// to allow using test tasks
#define BOOST_CHECK

#include "TaskCount.hpp"
#include "TaskConstParaboloid.hpp"
#include "TaskAdvection.hpp"

using namespace combigrid;

// this is necessary for correct function of task serialization
BOOST_CLASS_EXPORT(TaskAdvection)
BOOST_CLASS_EXPORT(TaskConstParaboloid)
BOOST_CLASS_EXPORT(TaskCount)
BOOST_CLASS_EXPORT(StaticFaults)
BOOST_CLASS_EXPORT(WeibullFaults)
BOOST_CLASS_EXPORT(FaultCriterion)


namespace shellCommand {
  // cf. https://stackoverflow.com/questions/478898/how-do-i-execute-a-command-and-get-the-output-of-the-command-within-c-using-po
  void exec(const char* cmd) {
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    sleep(2);         // wait for 2 seconds before closing
  }
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  /* when using timers (TIMING is defined in Stats), the Stats class must be
   * initialized at the beginning of the program. (and finalized in the end)
   */
  Stats::initialize();

  // read in parameter file
  std::string paramfile = "ctparam";
  if (argc > 1) paramfile = argv[1];
  boost::property_tree::ptree cfg;
  boost::property_tree::ini_parser::read_ini(paramfile, cfg);

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

    // read in third level parameters if available
    std::string thirdLevelHost, thirdLevelSSHCommand = "";
    unsigned int systemNumber = 0, numSystems = 1;
    unsigned short thirdLevelPort = 0;
    bool hasThirdLevel = static_cast<bool>(cfg.get_child_optional("thirdLevel"));
    if (hasThirdLevel) {
      std::cout << "Using third-level parallelism" << std::endl;
      thirdLevelHost = cfg.get<std::string>("thirdLevel.host");
      systemNumber = cfg.get<unsigned int>("thirdLevel.systemNumber");
      numSystems = cfg.get<unsigned int>("thirdLevel.numSystems");
      thirdLevelPort = cfg.get<unsigned short>("thirdLevel.port");
      thirdLevelSSHCommand = cfg.get<std::string>("thirdLevel.sshCommand");
    }

    // todo: read in boundary vector from ctparam
    std::vector<bool> boundary(dim, true);

    // check whether parallelization vector p agrees with nprocs
    IndexType checkProcs = 1;
    for (auto k : p) checkProcs *= k;
    assert(checkProcs == IndexType(nprocs));

    // set up the ssh tunnel for third level communication, if necessary
    // todo: if this works, move to ProcessManager::setUpThirdLevel
    if (thirdLevelSSHCommand != "") {
      shellCommand::exec(thirdLevelSSHCommand.c_str());
      std::cout << thirdLevelSSHCommand << " returned " << std::endl;
    }

    /* generate a list of levelvectors and coefficients
     * CombiMinMaxScheme will create a classical combination scheme.
     * however, you could also read in a list of levelvectors and coefficients
     * from a file */
    CombiMinMaxScheme combischeme(dim, lmin, lmax);
    combischeme.createAdaptiveCombischeme();
    std::vector<LevelVector> fullLevels = combischeme.getCombiSpaces();
    std::vector<combigrid::real> fullCoeffs = combischeme.getCoeffs();

    // split scheme and assign each half to a system
    std::vector<LevelVector> levels;
    std::vector<combigrid::real> coeffs;
    CombiThirdLevelScheme::createThirdLevelScheme(fullLevels, fullCoeffs, boundary, systemNumber,
                                                  numSystems, levels, coeffs);

    // create load model
    std::unique_ptr<LoadModel> loadmodel = std::unique_ptr<LoadModel>(new LinearLoadModel());
    // std::unique_ptr<LoadModel> loadmodel = std::unique_ptr<LoadModel>(new AnisotropyLoadModel());
    // std::unique_ptr<LoadModel> loadmodel = std::unique_ptr<LoadModel>(new
    // LearningLoadModel(levels));

    // output combination scheme
    std::cout << "lmin = " << lmin << std::endl;
    std::cout << "lmax = " << lmax << std::endl;
    std::cout << "CombiScheme: " << std::endl;
    for (const LevelVector& level : levels) std::cout << level << std::endl;

    // create Tasks
    TaskContainer tasks;
    std::vector<int> taskIDs;
    for (size_t i = 0; i < levels.size(); i++) {
      Task* t =
          new TaskAdvection(dim, levels[i], boundary, coeffs[i], loadmodel.get(), dt, nsteps, p);
      // Task* t = new TaskConstParaboloid(levels[i], boundary, coeffs[i], loadmodel);
      // Task* t = new TaskCount(dim, levels[i], boundary, coeffs[i], loadmodel.get());

      static_assert(!isGENE);

      tasks.push_back(t);
      taskIDs.push_back(t->getID());
    }

    // create combiparameters
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs, ncombi, 1, p,
                           std::vector<IndexType>(dim, 0), std::vector<IndexType>(dim, 0), thirdLevelHost,
                           thirdLevelPort, 0);

    // create abstraction for Manager
    ProcessManager manager(pgroups, tasks, params, std::move(loadmodel));

    manager.updateCombiParameters();
    std::cout << "set up component grids and run until first combination point" << std::endl;

    /* distribute task according to load model and start computation for
     * the first time */
    Stats::startEvent("manager run first");
    manager.runfirst();
    Stats::stopEvent("manager run first");

    // exchange subspace sizes to unify the dsgs in the third level case
    if (hasThirdLevel) {
      Stats::startEvent("manager unify subspace sizes with remote");
      manager.unifySubspaceSizesThirdLevel(),
          Stats::startEvent("manager unify subspace sizes with remote");
    }

    double start, finish;

    for (size_t i = 1; i < ncombi; ++i) {
      start = MPI_Wtime();

      Stats::startEvent("combine");
      if (hasThirdLevel) {
        manager.combineThirdLevel();
      } else {
        manager.combine();
      }
      Stats::stopEvent("combine");
      finish = MPI_Wtime();
      std::cout << "combination " << i << " took: " << finish - start << " seconds" << std::endl;

      // evaluate solution and
      // write solution to file
      std::string filename("out/solution_" + std::to_string(i) + ".raw");
      Stats::startEvent("manager write solution");
      manager.parallelEval(leval, filename, 0);
      Stats::stopEvent("manager write solution");

      std::cout << manager.parallelEvalNorm(leval, 0) << std::endl;

      std::cout << "run until combination point " << i + 1 << std::endl;

      // run tasks for next time interval
      start = MPI_Wtime();
      Stats::startEvent("manager run");
      manager.runnext();
      Stats::stopEvent("manager run");
      finish = MPI_Wtime();
      std::cout << "calculation " << i << " took: " << finish - start << " seconds" << std::endl;

      // run currently sets the dsgs back to zero
      // std::cout << manager.parallelEvalNorm(leval, 0) << std::endl;
    }

    Stats::startEvent("combine");
    if (hasThirdLevel) {
      manager.combineThirdLevel();
    } else {
      manager.combine();
    }
    Stats::stopEvent("combine");

    // evaluate solution and
    // write solution to file
    std::string filename("out/solution_" + std::to_string(ncombi) + ".raw");
    Stats::startEvent("manager write solution");
    manager.parallelEval(leval, filename, 0);
    Stats::stopEvent("manager write solution");

    std::cout << manager.getLpNorms(0) << std::endl;
    std::cout << manager.getLpNorms(1) << std::endl;
    std::cout << manager.getLpNorms(2) << std::endl;
    std::cout << "eval norms " << manager.parallelEvalNorm(leval, 0) << std::endl;

    auto analytical = manager.evalAnalyticalOnDFG(leval, 0);
    std::cout << "analytical " << analytical << std::endl;
    auto error = manager.evalErrorOnDFG(leval, 0);
    std::cout << "errors " << error << std::endl;

    std::cout << "relative errors ";
    for (size_t i=0; i < 3 ; ++i){
      std::cout << error[i]/analytical[i] << " ";
    }
    std::cout << std::endl;

    // TODO pollinta: for a massively parallel setting, this needs to be done differently
    FullGrid<CombiDataType> fg_eval(dim, leval, boundary);
    manager.gridEval(fg_eval);

    // exact solution
    TestFn f;
    // TestFnCount<CombiDataType> f;
    FullGrid<CombiDataType> fg_exact(dim, leval, boundary);
    fg_exact.createFullGrid();
    for (IndexType li = 0; li < fg_exact.getNrElements(); ++li) {
      std::vector<double> coords(dim);
      fg_exact.getCoords(li, coords);
      // fg_exact.getData()[li] = f(coords, ncombi);
      fg_exact.getData()[li] = f(coords, static_cast<double>(ncombi * nsteps) * dt);
    }

    // calculate error
    FullGrid<CombiDataType> fg_error(fg_exact);
    fg_error.add(fg_eval, -1);

    printf("Norms: %f %f\n", fg_eval.getlpNorm(0), fg_eval.getlpNorm(2));

    printf("Error: %f \n", fg_error.getlpNorm(0) / fg_exact.getlpNorm(0));
    printf("Error2: %f \n", fg_error.getlpNorm(2) / fg_exact.getlpNorm(2));

    // send exit signal to workers in order to enable a clean program termination
    manager.exit();
  }

  // this code is only executed by the worker processes
  else {
    // create abstraction of the process group from the worker's view
    ProcessGroupWorker pgroup;

    // wait for instructions from manager
    SignalType signal = -1;

    while (signal != EXIT){
      signal = pgroup.wait();
    } 
  }

  Stats::finalize();

  /* write stats to json file for postprocessing */
  Stats::write("timers.json");

  MPI_Finalize();

  return 0;
}
