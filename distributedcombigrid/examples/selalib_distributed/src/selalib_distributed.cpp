/*
 * SelalibTask.hpp
 *
 *  Created on: Nov 19, 2020
 *      Author: obersteiner
 */

#include <mpi.h>

#include <boost/filesystem.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/serialization/export.hpp>
#include <vector>

// compulsory includes for basic functionality
#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"
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
#include "sgpp/distributedcombigrid/utils/Stats.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"

// include user specific task. this is the interface to your application
#include "SelalibTask.hpp"

using namespace combigrid;
namespace fs = boost::filesystem;  // for file operations

// this is necessary for correct function of task serialization
#include "sgpp/distributedcombigrid/utils/BoostExports.hpp"
BOOST_CLASS_EXPORT(SelalibTask)

// helper funtion to read a bool vector from string
inline std::vector<bool>& operator>>(std::string str, std::vector<bool>& vec) {
  std::vector<std::string> strs;
  boost::split(strs, str, boost::is_any_of(" "));

  assert(vec.size() == strs.size());

  for (size_t i = 0; i < strs.size(); ++i) vec[i] = boost::lexical_cast<bool>(strs[i]);

  return vec;
}

// helper function to output bool vector
inline std::ostream& operator<<(std::ostream& os, const std::vector<bool>& l) {
  os << "[";

  for (size_t i = 0; i < l.size(); ++i) os << l[i] << " ";

  os << "]";

  return os;
}

// cf https://stackoverflow.com/questions/37931691/replace-a-word-in-text-file-using-c
std::string getFile(std::ifstream& is) {
  std::string contents;
  // Here is one way to read the whole file
  for (char ch; is.get(ch); contents.push_back(ch)) {
  }
  return contents;
}

// cf https://stackoverflow.com/questions/5878775/how-to-find-and-replace-string/5878802
std::string replaceFirstOccurrence(std::string& s, const std::string& toReplace,
                                   const std::string& replaceWith) {
  std::size_t pos = s.find(toReplace);
  if (pos == std::string::npos) return s;
  return s.replace(pos, toReplace.length(), replaceWith);
}


void setCheckpointRestart(std::string basename, std::vector<LevelVector> levels, std::string suffix = "") {
  std::string baseFolder = "./" + basename;
  for (size_t i = 0; i < levels.size(); i++) {
    // path to task folder
    std::string taskFolder = baseFolder + suffix + std::to_string(i);
    // assert that checkpoint is there
    std::string checkpointString = taskFolder + "/distribution_end-0000.h5";
    if (!fs::exists(checkpointString)) {
      throw std::runtime_error("No checkpoint to re-start from " + checkpointString);
    }
    {
      // adapt each parameter file
      std::ifstream inputFileStream(taskFolder + "/param.nml", std::ifstream::in);
      auto contents = getFile(inputFileStream);
      std::string newRestartInfo = "restart = .true.\n        restart_filename = \"distribution_end-0000.h5\"";
      contents = replaceFirstOccurrence(contents, "restart = .false.", newRestartInfo);
      std::ofstream outputFileStream(taskFolder + "/param_new.nml");
      outputFileStream << contents;
    }
    std::rename((taskFolder + "/param_new.nml").c_str(), (taskFolder + "/param.nml").c_str());
  }
}

bool createTaskFolders(std::string basename, std::vector<LevelVector> levels, IndexVector p,
                       size_t nsteps, double dt, size_t n_diagnostics, std::string suffix = "") {
  std::string baseFolder = "./" + basename;
  std::string templateFolder = "./template";
  assert(fs::exists(templateFolder));
  for (size_t i = 0; i < levels.size(); i++) {
    // path to task folder
    std::string taskFolder = baseFolder + suffix + std::to_string(i);
    // copy template directory
    if (!fs::exists(taskFolder) && !fs::create_directory(taskFolder)) {
      throw std::runtime_error("Cannot create destination directory " + taskFolder);
    }
    // adapt each parameter file
    std::ifstream inputFileStream(templateFolder + "/param.nml", std::ifstream::in);
    auto contents = getFile(inputFileStream);
    contents = replaceFirstOccurrence(contents, "$nsteps", std::to_string(nsteps));
    contents = replaceFirstOccurrence(contents, "$dt", std::to_string(dt));
    contents = replaceFirstOccurrence(contents, "$n_diagnostics", std::to_string(n_diagnostics));
    for (DimType d = 0; d < levels[0].size(); ++d) {
      contents =
          replaceFirstOccurrence(contents, "$nx" + std::to_string(d + 1),
                                 std::to_string(static_cast<LevelType>(std::pow(2, levels[i][d]))));
      contents =
          replaceFirstOccurrence(contents, "$p" + std::to_string(d + 1), std::to_string(p[d]));
    }
    std::ofstream outputFileStream(taskFolder + "/param.nml");
    outputFileStream << contents;
    // copy all other files
    for (const auto& dirEnt : fs::recursive_directory_iterator{templateFolder}) {
      const auto& path = dirEnt.path();
      auto relativePathStr = path.string();
      boost::replace_first(relativePathStr, templateFolder, "");
      fs::copy(path, taskFolder + "/" + relativePathStr, fs::copy_options::skip_existing);
    }
  }
  return true;
}

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
    IndexVector p(dim);
    std::vector<bool> boundary(dim), hierarchizationDims(dim);
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
    // cfg.get<std::string>("ct.leval2") >> leval2;  // level vector of second final output
    cfg.get<std::string>("ct.reduceCombinationDimsLmin") >> reduceCombinationDimsLmin;
    cfg.get<std::string>("ct.reduceCombinationDimsLmax") >> reduceCombinationDimsLmax;
    cfg.get<std::string>("ct.p") >> p;  // parallelization of domain (how many procs per dimension)
    cfg.get<std::string>("ct.boundary") >> boundary;  // which dimension have boundary points
    cfg.get<std::string>("ct.hierarchization_dims") >>
        hierarchizationDims;                // which dimension should be hierarchized
    std::string basis = cfg.get<std::string>("ct.basis", "hat");
    ncombi = cfg.get<size_t>("ct.ncombi");  // number of combinations
    std::string basename = cfg.get<std::string>("preproc.basename");
    dt = cfg.get<combigrid::real>(
        "application.dt");  // timestep
    nsteps = cfg.get<size_t>("application.nsteps");  // number of timesteps between combinations
    bool haveDiagnosticsTask = cfg.get<bool>("application.haveDiagnosticsTask", false);
    bool checkpointRestart = cfg.get<bool>("application.checkpoint_restart", false);

    std::string fg_file_path = cfg.get<std::string>("ct.fg_file_path");
    std::string fg_file_path2 = cfg.get<std::string>("ct.fg_file_path2");

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
    combischeme.createAdaptiveCombischeme();
    combischeme.makeFaultTolerant();
    levels = combischeme.getCombiSpaces();
    coeffs = combischeme.getCoeffs();

    // create load model
    // std::unique_ptr<LoadModel> loadmodel = std::unique_ptr<LoadModel>(new
    // LearningLoadModel(levels));
    std::unique_ptr<LoadModel> loadmodel = std::unique_ptr<LoadModel>(new LinearLoadModel());

    // output of combination setup
    std::cout << "lmin = " << lmin << std::endl;
    std::cout << "lmax = " << lmax << std::endl;
    std::cout << "reduceCombinationDimsLmin = " << reduceCombinationDimsLmin << std::endl;
    std::cout << "reduceCombinationDimsLmax = " << reduceCombinationDimsLmax << std::endl;
    std::cout << "boundary = " << boundary << std::endl;
    std::cout << "hierarchization_dims = " << hierarchizationDims << std::endl;
    auto numDOF = printCombiDegreesOfFreedom(combischeme.getCombiSpaces());
    std::cout << "CombiScheme: " << std::endl;
    for (size_t i = 0; i < levels.size(); ++i) {
      std::cout << "\t" << levels[i] << " " << coeffs[i] << std::endl;
    }

    if (checkpointRestart){
      // change the restart parameter in all existing parameter files in the folders to true
      setCheckpointRestart(basename, levels);
    } else {
      // using a very high diagnostics interval -> write no diagnostics in the component grids
      size_t veryHighNumber = 2147483647; // =2^31-1
      size_t sometimes = 100;
      size_t always = 1;
      // create necessary folders and files to run each task in a separate folder
      createTaskFolders(basename, levels, p, nsteps, dt, always);
    }

    // create Tasks
    TaskContainer tasks;
    std::vector<size_t> taskIDs;

    std::string baseFolder = "./" + basename;
    // initialize individual tasks (component grids)
    for (size_t i = 0; i < levels.size(); i++) {
      // path to task folder
      std::string path = baseFolder + std::to_string(i);
      Task* t = new SelalibTask(dim, levels[i], boundary, coeffs[i], loadmodel.get(), path, dt,
                                nsteps, p);
      tasks.push_back(t);
      taskIDs.push_back(t->getID());
    }
    Task* levalTask;
    if (haveDiagnosticsTask) {
      assert(reduceCombinationDimsLmax == LevelVector(dim, 0));
      // initialize diagnostics task
      if (checkpointRestart){
        // change the restart parameter
        setCheckpointRestart(basename, std::vector<LevelVector>(1, leval), "_leval_");
      } else {
        // create necessary folder and files to have diagnostics task in a different folder
        createTaskFolders(basename + "_leval_", std::vector<LevelVector>(1, leval), p, nsteps, dt, 1);
      }
      std::string path = baseFolder + "_leval_0";
      auto levalCoefficient = 0.;
      levalTask = new SelalibTask(dim, leval, boundary, levalCoefficient, loadmodel.get(), path, dt,
                                  nsteps, p);
      tasks.push_back(levalTask);
      taskIDs.push_back(levalTask->getID());
      levels.push_back(leval);
      coeffs.push_back(levalCoefficient);
    }

    // create combiparamters
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, hierarchizationDims, taskIDs,
                           ncombi, 1, reduceCombinationDimsLmin, reduceCombinationDimsLmax, false);
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
        manager.runfirst();
        Stats::stopEvent("manager run first");
      } else {
        // run tasks for next time interval
        Stats::startEvent("manager run");
        manager.runnext();
        Stats::stopEvent("manager run");
      }

      assert(!ENABLE_FT);
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

    std::cout << manager.getLpNorms(0) << std::endl;
    std::cout << manager.getLpNorms(1) << std::endl;
    std::cout << manager.getLpNorms(2) << std::endl;

    // evaluate solution on the grid defined by leval
    //(basically an interpolation of the sparse grid to fullgrid with resolution leval)
    // Stats::startEvent("manager parallel eval");
    // manager.parallelEval(leval, fg_file_path, 0);
    // Stats::stopEvent("manager parallel eval");
    // Stats::startEvent("manager parallel eval norm");
    // auto combinedNormLeval = manager.parallelEvalNorm(leval, 0);
    // Stats::stopEvent("manager parallel eval norm");

    // std::cout << "Norm is " << combinedNormLeval << " \n";

    // // evaluate solution on the grid defined by leval2
    // Stats::startEvent("manager parallel eval 2");
    // manager.parallelEval(leval2, fg_file_path2, 0);
    // Stats::stopEvent("manager parallel eval 2");

    // std::cout << "Computation finished evaluating on target grid! \n";

    manager.writeSparseGridMinMaxCoefficients(fg_file_path + "_selalib_sg");

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
