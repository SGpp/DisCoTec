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
// #include "sgpp/distributedcombigrid/combischeme/CombiThirdLevelScheme.hpp"
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
#include "sgpp/distributedcombigrid/utils/BoostExports.hpp"
BOOST_CLASS_EXPORT(TaskAdvection)


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


void managerMonteCarlo(ProcessManager& manager, DimType dim, double time, bool hasThirdLevel) {
      // Stats::startEvent("manager get norms");
      // std::cout << manager.getLpNorms(0) << std::endl;
      // std::cout << manager.getLpNorms(1) << std::endl;
      // std::cout << manager.getLpNorms(2) << std::endl;
      // std::cout << "eval norms " << manager.parallelEvalNorm(leval, 0) << std::endl;

      // auto analytical = manager.evalAnalyticalOnDFG(leval, 0);
      // std::cout << "analytical " << analytical << std::endl;
      // auto error = manager.evalErrorOnDFG(leval, 0);
      // std::cout << "errors " << error << std::endl;

      // std::cout << "relative errors ";
      // for (size_t i=0; i < 3 ; ++i){
      //   std::cout << error[i]/analytical[i] << " ";
      // }
      // std::cout << std::endl;
      // Stats::stopEvent("manager get norms");
      
      // 100000 was tested to be sufficient for the 6D blob,
      // but output three times just to make sure
      std::vector<size_t> numValuesToTry{100000};
      for (auto& numValues : numValuesToTry) {
        for (int i = 0; i < 3; ++i) {
          Stats::startEvent("manager monte carlo");
          // third-level monte carlo interpolation
          std::vector<std::vector<real>> interpolationCoords;
          std::vector<CombiDataType> values;
          if (hasThirdLevel) {
            // manager.monteCarloThirdLevel(numValues, interpolationCoords, values);
          } else {
            interpolationCoords = montecarlo::getRandomCoordinates(numValues, dim);
            values = manager.interpolateValues(interpolationCoords);
          }
          Stats::stopEvent("manager monte carlo");
  
          Stats::startEvent("manager calculate errors");
          // calculate monte carlo errors
          TestFn initialFunction;
          real l0Error = 0., l1Error = 0., l2Error = 0.,
               l0Reference = 0., l1Reference = 0., l2Reference = 0.;
          for (size_t i = 0; i < interpolationCoords.size(); ++i) {
            auto analyticalSln =
                initialFunction(interpolationCoords[i], time);
            l0Reference = std::max(analyticalSln, l0Reference);
            l1Reference += analyticalSln;
            l2Reference += std::pow(analyticalSln, 2);
            auto difference = std::abs(analyticalSln - values[i]);
            l0Error = std::max(difference, l0Error);
            l1Error += difference;
            l2Error += std::pow(difference, 2);
          }
          // make them relative errors
          l0Error = l0Error / l0Reference;
          l1Error = l1Error / l1Reference;
          l2Error = std::sqrt(l2Error) / std::sqrt(l2Reference);
          Stats::stopEvent("manager calculate errors");
  
          std::cout << "Monte carlo errors on " << numValues << " points are \n" << time
		    << ", " << l0Error << ", " << l1Error << ", " << l2Error << " " << std::endl;
        }
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
    std::string ctschemeFile = cfg.get<std::string>("ct.ctscheme", "");
    std::string basis = cfg.get<std::string>("ct.basis", "hat");
    dt = cfg.get<combigrid::real>("application.dt");
    nsteps = cfg.get<size_t>("application.nsteps");
    bool evalMCError = cfg.get<bool>("application.mcerror", false);
    std::string output_id = basis + "_" + std::to_string(dim) + "D_" + std::to_string(lmin[0]) +
                            "-" + std::to_string(lmax[0]);

    // todo: read in boundary vector from ctparam
    std::vector<BoundaryType> boundary(dim, 2);

    // check whether parallelization vector p agrees with nprocs
    IndexType checkProcs = 1;
    for (auto k : p) checkProcs *= k;
    if (checkProcs != IndexType(nprocs)) {
      throw std::invalid_argument("specified processes do not match! " +
                                  std::to_string(checkProcs) + " vs " + std::to_string(nprocs));
    }

    CombiMinMaxScheme combischeme(dim, lmin, lmax);
    combischeme.createAdaptiveCombischeme();
    std::vector<LevelVector> levels = combischeme.getCombiSpaces();
    std::vector<combigrid::real> coeffs = combischeme.getCoeffs();

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
    std::vector<size_t> taskIDs;
    for (size_t i = 0; i < levels.size(); i++) {
      Task* t =
          new TaskAdvection(dim, levels[i], boundary, coeffs[i], loadmodel.get(), dt, nsteps, p);
      // Task* t = new TaskConstParaboloid(levels[i], boundary, coeffs[i], loadmodel);
      // Task* t = new TaskCount(dim, levels[i], boundary, coeffs[i], loadmodel.get());

      static_assert(!isGENE, "isGENE");

      tasks.push_back(t);
      taskIDs.push_back(t->getID());
    }

    // create combiparameters
    auto reduceCombinationDimsLmax = std::vector<IndexType>(dim, 0);
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs, ncombi);
    setCombiParametersHierarchicalBasesUniform(params, basis);
    params.setParallelization(p);

    // create abstraction for Manager
    ProcessManager manager(pgroups, tasks, params, std::move(loadmodel));

    manager.updateCombiParameters();
    std::cout << "set up component grids and run until first combination point" << std::endl;

    /* distribute task according to load model and start computation for
     * the first time */
    Stats::startEvent("manager run first");
    manager.runfirst();
    Stats::stopEvent("manager run first");

    // double start, finish;

    std::cout << "time maxnorm l1norm " << std::endl;
    std::cout << std::setprecision(std::numeric_limits<real>::digits10 + 1);

    for (size_t i = 1; i < ncombi; ++i) {
      // start = MPI_Wtime();
      if (tasks.size() > 1) {
        Stats::startEvent("manager combine");
        manager.combine();
        // manager.waitAllFinished();
        Stats::stopEvent("manager combine");
      }
      if (evalMCError && i%1000 == 0) {
      	managerMonteCarlo(manager, dim, static_cast<double>(i * nsteps) * dt, false);
      }
      // write out field at middle of simulation
      if (i == ncombi/2 && dim < 4) {
        std::string filename("out/solution_" + output_id + "_" + std::to_string(i) + ".raw");
        Stats::startEvent("manager write solution");
        manager.parallelEval(leval, filename, 0);
        Stats::stopEvent("manager write solution");
      }
      // finish = MPI_Wtime();
      // std::cout << "combination " << i << " took: " << finish - start << " seconds" << std::endl;

      // evaluate solution and
      // // write solution to file
      // // every ten times
      // if (i%(ncombi/10) == 0) {
      //   std::string filename("out/solution_" + std::to_string(i) + ".raw");
      //   Stats::startEvent("manager write solution");
      //   manager.parallelEval(leval, filename, 0);
      //   Stats::stopEvent("manager write solution");
      //   std::string filename2("out/solution_" + std::to_string(i) + ".vtk");
      //   Stats::startEvent("manager write solution");
      //   manager.parallelEval(leval, filename2, 0);
      //   Stats::stopEvent("manager write solution");
      // }

      // std::cout << manager.parallelEvalNorm(leval, 0) << std::endl;
      // auto error = manager.evalErrorOnDFG(leval, 0);
      // std::cout << "errors " << error << std::endl;

      if (i % 100 == 0) {
        Stats::startEvent("manager get norms");
        // std::cout << "Max norms " << manager.getLpNorms(0) << std::endl;
        // std::cout << "L1 norms " << manager.getLpNorms(1) << std::endl;
        std::cout <<  " " << i * dt << " " << manager.getLpNorms(0)[0] << " " << manager.getLpNorms(1)[0] << std::endl;
        Stats::stopEvent("manager get norms");
      }
      // run tasks for next time interval
      // start = MPI_Wtime();
      Stats::startEvent("manager run");
      manager.runnext();
      Stats::stopEvent("manager run");
      // finish = MPI_Wtime();
      // std::cout << "calculation " << i << " took: " << finish - start << " seconds" << std::endl;

      // run currently sets the dsgs back to zero
      // std::cout << manager.parallelEvalNorm(leval, 0) << std::endl;
    }

    if (tasks.size() > 1) {
      Stats::startEvent("combine");
      manager.combine();
      Stats::stopEvent("combine");
    }

    if (evalMCError) {
      managerMonteCarlo(manager, dim, static_cast<double>(ncombi * nsteps) * dt, false);
    }

    // evaluate solution and
    // write solution to file
    if (dim < 4) {
      std::string filename("out/solution_" + output_id + "_" + std::to_string(ncombi) + ".raw");
      Stats::startEvent("manager write solution");
      manager.parallelEval(leval, filename, 0);
      Stats::stopEvent("manager write solution");
    }
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
