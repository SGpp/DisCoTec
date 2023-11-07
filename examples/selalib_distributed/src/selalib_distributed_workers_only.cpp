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
#include "combischeme/CombiMinMaxScheme.hpp"
#include "fault_tolerance/FaultCriterion.hpp"
#include "fault_tolerance/StaticFaults.hpp"
#include "fault_tolerance/WeibullFaults.hpp"
#include "fullgrid/FullGrid.hpp"
#include "io/BroadcastParameters.hpp"
#include "loadmodel/LinearLoadModel.hpp"
#include "manager/CombiParameters.hpp"
#include "manager/ProcessGroupManager.hpp"
#include "manager/ProcessGroupWorker.hpp"
#include "manager/ProcessManager.hpp"
#include "task/Task.hpp"
#include "utils/Stats.hpp"
#include "utils/Types.hpp"

// include user specific task. this is the interface to your application
#include "SelalibFileTools.hpp"
#include "SelalibTask.hpp"

using namespace combigrid;

// this is necessary for correct function of task serialization
#include "utils/BoostExports.hpp"
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
  theMPISystem()->initWorldReusable(MPI_COMM_WORLD, ngroup, nprocs, false, true);

  DimType dim = cfg.get<DimType>("ct.dim");
  LevelVector lmin(dim), lmax(dim), reduceCombinationDimsLmin(dim),
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
  cfg.get<std::string>("ct.lmin") >> lmin;    // minimal level vector for each grid
  cfg.get<std::string>("ct.lmax") >> lmax;    // maximum level vector -> level vector of target grid
  cfg.get<std::string>("ct.reduceCombinationDimsLmin") >> reduceCombinationDimsLmin;
  cfg.get<std::string>("ct.reduceCombinationDimsLmax") >> reduceCombinationDimsLmax;
  cfg.get<std::string>("ct.p") >> p;  // parallelization of domain (how many procs per dimension)
  cfg.get<std::string>("ct.hierarchization_dims") >>
      hierarchizationDims;  // which dimension should be hierarchized
  std::string basis = cfg.get<std::string>("ct.basis", "hat_periodic");
  ncombi = cfg.get<size_t>("ct.ncombi");  // number of combinations
  std::string basename = cfg.get<std::string>("preproc.basename");
  dt = cfg.get<combigrid::real>("application.dt");  // timestep
  nsteps = cfg.get<size_t>("application.nsteps");   // number of timesteps between combinations
  bool haveDiagnosticsTask = cfg.get<bool>("application.haveDiagnosticsTask", false);
  bool checkpointRestart = cfg.get<bool>("application.checkpoint_restart", false);
  std::string nameDiagnostics = cfg.get<std::string>("application.name_diagnostics", "vp_B2_3d3v");
  bool haveResolution = static_cast<bool>(cfg.get_child_optional("ct.resolution"));
  if (haveResolution) {
    cfg.get<std::string>("ct.resolution") >>
        resolution;  // resolution of single grid in the full grid case
    assert(lmax == lmin);
  }

  std::string fg_file_path = cfg.get<std::string>("ct.fg_file_path");

  // check parallelization vector p agrees with nprocs
  IndexType checkProcs = 1;
  for (auto k : p) checkProcs *= k;
  assert(checkProcs == IndexType(nprocs));

  std::vector<LevelVector> levels;
  std::vector<combigrid::real> coeffs;
  std::vector<size_t> taskNumbers;  // only used in case of static task assignment
  const auto& pgroupNumber = theMPISystem()->getProcessGroupNumber();
  // read in CT scheme
  {
    std::unique_ptr<CombiMinMaxScheme> scheme(new CombiMinMaxScheme(dim, lmin, lmax));
    scheme->createClassicalCombischeme();
    size_t totalNumTasks =
        combigrid::getLoadBalancedLevels(*scheme, pgroupNumber, theMPISystem()->getNumGroups(),
                                         boundary, levels, coeffs, taskNumbers);

    assert(!levels.empty());
    assert(levels.size() == coeffs.size());
    assert(levels.size() == taskNumbers.size());

    MASTER_EXCLUSIVE_SECTION {
      std::cout << getTimeStamp() << " Process group " << pgroupNumber << " will run "
                << levels.size() << " of " << totalNumTasks << " tasks." << std::endl;
      printCombiDegreesOfFreedom(levels, boundary);
    }
  }

  // create load model
  // std::unique_ptr<LoadModel> loadmodel = std::unique_ptr<LoadModel>(new
  // LearningLoadModel(levels));
  std::unique_ptr<LoadModel> loadmodel = std::unique_ptr<LoadModel>(new LinearLoadModel());

  // output of combination setup
  MIDDLE_PROCESS_EXCLUSIVE_SECTION {
    if (haveResolution) {
      std::cout << "resolution = " << resolution << std::endl;
    } else {
      std::cout << "lmin = " << lmin << std::endl;
      std::cout << "lmax = " << lmax << std::endl;
      std::cout << "reduceCombinationDimsLmin = " << reduceCombinationDimsLmin << std::endl;
      std::cout << "reduceCombinationDimsLmax = " << reduceCombinationDimsLmax << std::endl;
      std::cout << "boundary = " << boundary << std::endl;
      std::cout << "hierarchization_dims = " << hierarchizationDims << std::endl;
      std::cout << "CombiScheme: " << std::endl;
      for (size_t i = 0; i < levels.size(); ++i) {
        std::cout << "\t" << levels[i] << " " << coeffs[i] << std::endl;
      }
    }
  }
  MASTER_EXCLUSIVE_SECTION {
    if (checkpointRestart) {
      // change the restart parameter in all existing parameter files in the folders to true
      setCheckpointRestart(basename, levels);
    } else {
      // using a very high diagnostics interval -> write no diagnostics in the component grids
      [[maybe_unused]] size_t veryHighNumber = 2147483647;  // =2^31-1
      [[maybe_unused]] size_t sometimes = 100;
      [[maybe_unused]] size_t always = 1;
      // create necessary folders and files to run each task in a separate folder
      createTaskFolders(basename, levels, taskNumbers, p, nsteps, dt, veryHighNumber,
                        nameDiagnostics);
      if (haveResolution) {
        assert(levels.size() == 1);
        assert(!haveDiagnosticsTask);
        adaptParameterFileFirstFolder(basename, resolution, p, nsteps, dt, veryHighNumber,
                                      nameDiagnostics);
        // if haveResolution: set infty coefficient, does not need to be combined anyways
        coeffs[0] = std::numeric_limits<combigrid::real>::max();
      }
    }
  }
  {
    ProcessGroupWorker worker;

    // create combiparamters
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, hierarchizationDims,
                           taskNumbers, ncombi, 1, CombinationVariant::subspaceReduce,
                           reduceCombinationDimsLmin, reduceCombinationDimsLmax, 128, false);
    setCombiParametersHierarchicalBasesUniform(params, basis);
    params.setParallelization(p);

    worker.setCombiParameters(std::move(params));

    // create Tasks

    std::string baseFolder = "./" + basename;
    // initialize individual tasks (component grids)
    for (size_t i = 0; i < levels.size(); i++) {
      // path to task folder
      std::string path = baseFolder + std::to_string(taskNumbers[i]);
      worker.initializeTask(std::unique_ptr<Task>(
          new SelalibTask(levels[i], boundary, coeffs[i], loadmodel.get(), path, dt, nsteps, p)));
    }
    // worker.initializeAllTasks<SelalibTask>(levels, coeffs, taskNumbers, loadmodel.get(), dt,
    // path,
    //
    assert(!haveDiagnosticsTask);

    worker.initCombinedDSGVector();
    worker.zeroDsgsData();

    MASTER_EXCLUSIVE_SECTION {
      std::cout << getTimeStamp() << "group " << theMPISystem()->getProcessGroupNumber()
                << " sparse grid, will allocate "
                << static_cast<real>(worker.getCombinedDSGVector()[0]->getAccumulatedDataSize() *
                                     sizeof(CombiDataType)) /
                       1e6
                << " MB" << std::endl;
    }
    MPI_Barrier(theMPISystem()->getWorldComm());
    std::chrono::high_resolution_clock::time_point after_init_time =
        std::chrono::high_resolution_clock::now();
    MIDDLE_PROCESS_EXCLUSIVE_SECTION
    std::cout << getTimeStamp() << "initialization took "
              << std::chrono::duration_cast<std::chrono::microseconds>(after_init_time - init_time)
                     .count()
              << "us " << std::endl;

    // start computation
    // we perform ncombi many combinations with
    // fixed stepsize or simulation time between each combination
    for (size_t i = 0; i < ncombi; ++i) {
      // run tasks for next time interval
      worker.runAllTasks();
      MIDDLE_PROCESS_EXCLUSIVE_SECTION
      std::cout << getTimeStamp() << "run took " << Stats::getDuration("run") / 1000.0
                << "seconds " << std::endl;
      static_assert(!ENABLE_FT);
      if (!haveResolution) {
        // combine grids
        Stats::startEvent("worker combine");
        worker.combineAtOnce(true);
        Stats::stopEvent("worker combine");
        MIDDLE_PROCESS_EXCLUSIVE_SECTION
        std::cout << getTimeStamp() << "worker combine took "
                  << Stats::getDuration("worker combine") / 1000.0 << "seconds " << std::endl;

        if (i % 100 == 0) {
          worker.writeSparseGridMinMaxCoefficients(fg_file_path + "_selalib_sg_" +
                                                   std::to_string(i));
          Stats::writePartial("stats_worker.json");
        }
      }
    }

    if (!haveResolution) {
      worker.writeSparseGridMinMaxCoefficients(fg_file_path + "_selalib_sg");
      std::cout << worker.getLpNorms(0) << std::endl;
      std::cout << worker.getLpNorms(1) << std::endl;
      std::cout << worker.getLpNorms(2) << std::endl;
    }
  }

  MPI_Barrier(theMPISystem()->getWorldComm());
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
