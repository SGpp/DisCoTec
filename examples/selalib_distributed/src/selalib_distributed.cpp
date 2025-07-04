/*
 * SelalibTask.hpp
 *
 *  Created on: Nov 19, 2020
 *      Author: obersteiner
 */

#include <mpi.h>

#include <boost/serialization/export.hpp>
#include <vector>

// compulsory includes for basic functionality
#include "../../../include/discotec/combischeme/CombiMinMaxScheme.hpp"
#include "../../../include/discotec/fault_tolerance/FaultCriterion.hpp"
#include "../../../include/discotec/fault_tolerance/StaticFaults.hpp"
#include "../../../include/discotec/fault_tolerance/WeibullFaults.hpp"
#include "../../../include/discotec/fullgrid/FullGrid.hpp"
#include "../../../include/discotec/io/BroadcastParameters.hpp"
#include "../../../include/discotec/loadmodel/LinearLoadModel.hpp"
#include "../../../include/discotec/manager/CombiParameters.hpp"
#include "../../../include/discotec/manager/ProcessGroupManager.hpp"
#include "../../../include/discotec/manager/ProcessGroupWorker.hpp"
#include "../../../include/discotec/manager/ProcessManager.hpp"
#include "../../../include/discotec/Task.hpp"
#include "../../../include/discotec/utils/Stats.hpp"
#include "../../../include/discotec/utils/Types.hpp"

// include user specific task. this is the interface to your application
#include "SelalibFileTools.hpp"
#include "SelalibTask.hpp"

using namespace combigrid;

// this is necessary for correct function of task serialization
#include "../../../include/discotec/utils/BoostExports.hpp"
BOOST_CLASS_EXPORT(SelalibTask)

void initMpiSelalibStyle(int argc, char** argv) {
  sll_s_allocate_collective();
  int ignore;
#ifdef _OPENMP
  int mpi_mode = MPI_THREAD_MULTIPLE;
#else
  int mpi_mode = MPI_THREAD_SINGLE;
#endif
  MPI_Init_thread(&argc, &argv, mpi_mode, &ignore);
}

int main(int argc, char** argv) {
  std::chrono::high_resolution_clock::time_point init_time =
      std::chrono::high_resolution_clock::now();
  Stats::initialize();

  if (ENABLE_FT) {
    assert(false &&
           "this example is not adapted to fault tolerance; look at gene_distributed if you want "
           "to implement that");
  } else {
    initMpiSelalibStyle(argc, argv);
  }

  // only one rank reads parameter file and broadcasts to others
  std::string paramfile = "ctparam";
  if (argc > 1) paramfile = argv[1];
  boost::property_tree::ptree cfg =
      broadcastParameters::getParametersFromRankZero(paramfile, MPI_COMM_WORLD);

  // number of process groups and number of processes per group
  size_t ngroup = cfg.get<size_t>("manager.ngroup");
  size_t nprocs = cfg.get<size_t>("manager.nprocs");

  // divide the MPI processes into process group and initialize the
  // corresponding communicators
  theMPISystem()->initWorldReusable(MPI_COMM_WORLD, ngroup, nprocs, true, true);

  WORLD_MANAGER_EXCLUSIVE_SECTION {
    Stats::startEvent("manager initialization");
    // manager code
    /* create an abstraction of the process groups for the manager's view
     * a pgroup is identified by the ID in gcomm
     */
    ProcessGroupManagerContainer pgroups;
    // create vector containing the different process groups
    for (size_t i = 0; i < ngroup; ++i) {
      // todo: order of ranks in new group?
      int pgroupRootID(i);
      pgroups.emplace_back(std::make_shared<ProcessGroupManager>(pgroupRootID));
    }

    /* generate a list of levelvectors and coefficients
     * CombiTS_CT will generate a valid combination. however, you could
     * also read in a list of levelvectors and coefficients from a file */
    DimType dim = cfg.get<DimType>("ct.dim");
    LevelVector lmin(dim), lmax(dim), leval(dim), leval2(dim), reduceCombinationDimsLmin(dim),
        reduceCombinationDimsLmax(dim);
    std::vector<int> p(dim), resolution(dim);
    std::vector<bool> hierarchizationDims(dim);
    std::vector<BoundaryType> boundary(dim, 1);
    combigrid::real dt;
    // time inteveral of 1 combination
    // only necessary if number of timesteps varies for each grid
    // otherwise set very high and use ntimesteps to adjust combiinterval
    // combigrid::real combitime;
    // read combination parameters
    size_t nsteps, ncombi;
    cfg.get<std::string>("ct.lmin") >> lmin;  // minimal level vector for each grid
    cfg.get<std::string>("ct.lmax") >> lmax;  // maximum level vector -> level vector of target grid
    cfg.get<std::string>("ct.leval") >> leval;    // level vector of final output
    cfg.get<std::string>("ct.reduceCombinationDimsLmin") >> reduceCombinationDimsLmin;
    cfg.get<std::string>("ct.reduceCombinationDimsLmax") >> reduceCombinationDimsLmax;
    cfg.get<std::string>("ct.p") >> p;  // parallelization of domain (how many procs per dimension)
    cfg.get<std::string>("ct.hierarchization_dims") >>
        hierarchizationDims;                // which dimension should be hierarchized
    std::string basis = cfg.get<std::string>("ct.basis", "hat_periodic");
    ncombi = cfg.get<size_t>("ct.ncombi");  // number of combinations
    std::string basename = cfg.get<std::string>("preproc.basename");
    dt = cfg.get<combigrid::real>(
        "application.dt");  // timestep
    nsteps = cfg.get<size_t>("application.nsteps");  // number of timesteps between combinations
    bool haveDiagnosticsTask = cfg.get<bool>("application.haveDiagnosticsTask", false);
    bool checkpointRestart = cfg.get<bool>("application.checkpoint_restart", false);
    std::string nameDiagnostics = cfg.get<std::string>("application.name_diagnostics" ,"vp_B2_3d3v");
    bool haveResolution = static_cast<bool>(cfg.get_child_optional("ct.resolution"));
    if (haveResolution) {
      cfg.get<std::string>("ct.resolution") >> resolution; // resolution of single grid in the full grid case
      assert(lmax == lmin);
    }

    std::string fg_file_path = cfg.get<std::string>("ct.fg_file_path");

    cfg.get<std::string>("ct.lmin") >> lmin;

    // check parallelization vector p agrees with nprocs
    IndexType checkProcs = 1;

    for (auto k : p) checkProcs *= k;

    assert(checkProcs == IndexType(nprocs));

    // read combi levels and spaces from file as long as min max scheme does
    // not work properly
    std::vector<LevelVector> levels;
    std::vector<combigrid::real> coeffs;
    // std::vector<int> fileTaskIDs;

    CombiMinMaxScheme combischeme(dim, lmin, lmax);
    combischeme.createClassicalCombischeme();
    // combischeme.makeFaultTolerant();
    levels = combischeme.getCombiSpaces();
    coeffs = combischeme.getCoeffs();

    // create load model
    // std::unique_ptr<LoadModel> loadmodel = std::unique_ptr<LoadModel>(new
    // LearningLoadModel(levels));
    std::unique_ptr<LoadModel> loadmodel = std::unique_ptr<LoadModel>(new LinearLoadModel());

    // output of combination setup
    if (haveResolution) {
      std::cout << "resolution = " << resolution << std::endl;
    } else {
      std::cout << "lmin = " << lmin << std::endl;
      std::cout << "lmax = " << lmax << std::endl;
      std::cout << "reduceCombinationDimsLmin = " << reduceCombinationDimsLmin << std::endl;
      std::cout << "reduceCombinationDimsLmax = " << reduceCombinationDimsLmax << std::endl;
      std::cout << "boundary = " << boundary << std::endl;
      std::cout << "hierarchization_dims = " << hierarchizationDims << std::endl;
      [[maybe_unused]] auto numDOF = printCombiDegreesOfFreedom(combischeme.getCombiSpaces(), boundary);
      std::cout << "CombiScheme: " << std::endl;
      for (size_t i = 0; i < levels.size(); ++i) {
        std::cout << "\t" << levels[i] << " " << coeffs[i] << std::endl;
      }
    }

    if (checkpointRestart){
      // change the restart parameter in all existing parameter files in the folders to true
      setCheckpointRestart(basename, levels);
    } else {
      // using a very high diagnostics interval -> write no diagnostics in the component grids
      [[maybe_unused]] size_t veryHighNumber = 2147483647; // =2^31-1
      [[maybe_unused]] size_t sometimes = 100;
      [[maybe_unused]] size_t always = 1;
      // create necessary folders and files to run each task in a separate folder
      std::vector<size_t> taskNumbers(levels.size());  // only necessary in case of static task assignment
      std::iota(taskNumbers.begin(), taskNumbers.end(), 0);
      createTaskFolders(basename, levels, taskNumbers, p, nsteps, dt, always, nameDiagnostics);
      if (haveResolution) {
        assert(levels.size() == 1);
        assert(!haveDiagnosticsTask);
        adaptParameterFileFirstFolder(basename, resolution, p, nsteps, dt, always, nameDiagnostics);
        // if haveResolution: set infty coefficient, does not need to be combined anyways
        coeffs[0] = std::numeric_limits<combigrid::real>::max();
      }
    }

    // create Tasks
    TaskContainer tasks;
    std::vector<size_t> taskIDs;

    std::string baseFolder = "./" + basename;
    // initialize individual tasks (component grids)
    for (size_t i = 0; i < levels.size(); i++) {
      // path to task folder
      std::string path = baseFolder + std::to_string(i);
      Task* t =
          new SelalibTask(levels[i], boundary, coeffs[i], loadmodel.get(), path, dt, nsteps, p);
      tasks.push_back(t);
      taskIDs.push_back(t->getID());
    }
    Task* levalTask;
    if (haveDiagnosticsTask) {
      assert(reduceCombinationDimsLmax == LevelVector(dim, 0));
      // initialize diagnostics task
      if (checkpointRestart) {
        // change the restart parameter
        setCheckpointRestart(basename, std::vector<LevelVector>(1, leval), "_leval_");
      } else {
        // create necessary folder and files to have diagnostics task in a different folder
        std::vector<size_t> taskNumbers(1, 0); 
        createTaskFolders(basename + "_leval_", std::vector<LevelVector>(1, leval), taskNumbers, p, nsteps, dt,
                          1, nameDiagnostics);
      }
      std::string path = baseFolder + "_leval_0";
      auto levalCoefficient = 0.;
      levalTask =
          new SelalibTask(leval, boundary, levalCoefficient, loadmodel.get(), path, dt, nsteps, p);
      tasks.push_back(levalTask);
      taskIDs.push_back(levalTask->getID());
      levels.push_back(leval);
      coeffs.push_back(levalCoefficient);
    }

    // create combiparamters
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, hierarchizationDims, taskIDs,
                           ncombi, 1, CombinationVariant::sparseGridReduce,
                           reduceCombinationDimsLmin, reduceCombinationDimsLmax, 64, false);
    setCombiParametersHierarchicalBasesUniform(params, basis);
    params.setParallelization(p);

    // create Manager with process groups
    ProcessManager manager(pgroups, tasks, params, std::move(loadmodel));

    // combiparameters need to be set before starting the computation
    Stats::startEvent("update combi parameters");
    manager.updateCombiParameters();
    Stats::stopEvent("update combi parameters");

    std::chrono::high_resolution_clock::time_point after_init_time =
        std::chrono::high_resolution_clock::now();
    std::cout << "initialization took "
              << std::chrono::duration_cast<std::chrono::microseconds>(after_init_time - init_time)
                     .count()
              << "us " << std::endl;
    Stats::stopEvent("manager initialization");

    // start computation
    // we perform ncombi many combinations with
    // fixed stepsize or simulation time between each combination
    for (size_t i = 0; i < ncombi; ++i) {
      if (i == 0) {
        /* distribute task according to load model and start computation for
         * the first time */
        Stats::startEvent("manager run first");
        manager.runfirst(!haveResolution);
        Stats::stopEvent("manager run first");
      } else {
        // run tasks for next time interval
        Stats::startEvent("manager run");
        manager.runnext();
        Stats::stopEvent("manager run");
      }

      assert(!ENABLE_FT);
      if(!haveResolution) {
        // combine grids
        Stats::startEvent("manager combine");
        manager.combine();
        Stats::stopEvent("manager combine");

        if (haveDiagnosticsTask){
          Stats::startEvent("manager diag selalib");
          manager.doDiagnostics(levalTask->getID());
          Stats::stopEvent("manager diag selalib");
        }
        if (i%100 == 0) {
          manager.writeSparseGridMinMaxCoefficients(fg_file_path + "_selalib_sg_" + std::to_string(i));
        }
      }
    }

    // std::cout << "Computation finished evaluating on target grid! \n";

    if(!haveResolution) {
      manager.writeSparseGridMinMaxCoefficients(fg_file_path + "_selalib_sg");
      std::cout << manager.getLpNorms(0) << std::endl;
      std::cout << manager.getLpNorms(1) << std::endl;
      std::cout << manager.getLpNorms(2) << std::endl;
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

    while (signal != EXIT) {
      signal = pgroup.wait();
      // std::cout << "Worker worked on signal " << signal << "! \n";
    }
  }

  // finalize timing evaluations
  Stats::finalize();
  /* write stats to json file for postprocessing */
  Stats::write("timers.json");
  MPI_Finalize();

  std::chrono::high_resolution_clock::time_point end_time =
      std::chrono::high_resolution_clock::now();
  std::cout << "total "
            << std::chrono::duration_cast<std::chrono::microseconds>(end_time - init_time).count()
            << " " << std::endl;

  return 0;
}
