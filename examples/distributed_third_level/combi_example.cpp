/*
 * combi_example.cpp
 *
 *  Created on: Sep 23, 2015
 *      Author: heenemo
 */
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <boost/asio.hpp>
#include <boost/serialization/export.hpp>
#include <string>
#include <vector>

// compulsory includes for basic functionality
#include "combischeme/CombiMinMaxScheme.hpp"
#include "combischeme/CombiThirdLevelScheme.hpp"
#include "fault_tolerance/FaultCriterion.hpp"
#include "fault_tolerance/StaticFaults.hpp"
#include "fault_tolerance/WeibullFaults.hpp"
#include "io/BroadcastParameters.hpp"
#include "loadmodel/LinearLoadModel.hpp"
#include "manager/CombiParameters.hpp"
#include "manager/ProcessGroupManager.hpp"
#include "manager/ProcessGroupWorker.hpp"
#include "manager/ProcessManager.hpp"
#include "task/Task.hpp"
#include "utils/MonteCarlo.hpp"
#include "utils/Types.hpp"
#include "mpi/MPIMemory.hpp"
// include user specific task. this is the interface to your application

// to allow using test tasks
#define BOOST_CHECK

#include "TaskAdvection.hpp"

using namespace combigrid;

// this is necessary for correct function of task serialization
#include "utils/BoostExports.hpp"
BOOST_CLASS_EXPORT(TaskAdvection)

namespace shellCommand {
// cf.
// https://stackoverflow.com/questions/478898/how-do-i-execute-a-command-and-get-the-output-of-the-command-within-c-using-po
void exec(const char* cmd) {
  std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
  if (!pipe) {
    throw std::runtime_error("popen() failed!");
  }
  sleep(2);  // wait for 2 seconds before closing
}
}  // namespace shellCommand

void managerMonteCarlo(ProcessManager& manager, DimType dim, double time, bool hasThirdLevel) {
  // 100000 was tested to be sufficient for the 6D blob,
  // but output three times just to make sure
  std::vector<size_t> numValuesToTry{100000};
  for (auto& numValues : numValuesToTry) {
    for (int i = 0; i < 1; ++i) {
      Stats::startEvent("manager monte carlo");
      // third-level monte carlo interpolation
      std::vector<std::vector<real>> interpolationCoords;
      std::vector<CombiDataType> values;
      if (hasThirdLevel) {
        manager.monteCarloThirdLevel(numValues, interpolationCoords, values);
      } else {
        interpolationCoords = montecarlo::getRandomCoordinates(numValues, dim);
        values = manager.interpolateValues(interpolationCoords);
      }
      Stats::stopEvent("manager monte carlo");

      Stats::startEvent("manager calculate errors");
      // calculate monte carlo errors
      TestFn initialFunction;
      real l0Error = 0., l1Error = 0., l2Error = 0., l0Reference = 0., l1Reference = 0.,
           l2Reference = 0.;
      for (size_t i = 0; i < interpolationCoords.size(); ++i) {
        auto analyticalSln = initialFunction(interpolationCoords[i], time);
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

      std::cout << "Monte carlo errors on " << numValues << " points are \n"
                << time << ", " << l0Error << ", " << l1Error << ", " << l2Error << " "
                << std::endl;
    }
  }
}

int main(int argc, char** argv) {
  [[maybe_unused]] auto mpiOnOff = MpiOnOff(&argc, &argv);

  /* when using timers (TIMING is defined in Stats), the Stats class must be
   * initialized at the beginning of the program. (and finalized in the end)
   */
  Stats::initialize();

  // only one rank reads parameter file and broadcasts to others
  std::string paramfile = "ctparam";
  if (argc > 1) paramfile = argv[1];
  boost::property_tree::ptree cfg =
      broadcastParameters::getParametersFromRankZero(paramfile, MPI_COMM_WORLD);

  // number of process groups and number of processes per group
  size_t ngroup = cfg.get<size_t>("manager.ngroup");
  size_t nprocs = cfg.get<size_t>("manager.nprocs");

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
  uint32_t chunkSizeInMebibyte = cfg.get<uint32_t>("ct.chunkSize", 64);
  std::string ctschemeFile = cfg.get<std::string>("ct.ctscheme", "");
  dt = cfg.get<combigrid::real>("application.dt");
  nsteps = cfg.get<size_t>("application.nsteps");
  bool evalMCError = cfg.get<bool>("application.mcerror", false);

  // read in third level parameters if available
  std::string thirdLevelHost, thirdLevelSSHCommand = "";
  unsigned int systemNumber = 0, numSystems = 1;
  unsigned short thirdLevelPort = 0;
  bool hasThirdLevel = static_cast<bool>(cfg.get_child_optional("thirdLevel"));
  bool extraSparseGrid = true;
  std::vector<real> fractionsOfScheme;
  bool brokerOnSameSystem = false;
  if (hasThirdLevel) {
    std::cout << "Using third-level parallelism" << std::endl;
    thirdLevelHost = cfg.get<std::string>("thirdLevel.host");
    systemNumber = cfg.get<unsigned int>("thirdLevel.systemNumber");
    numSystems = cfg.get<unsigned int>("thirdLevel.numSystems");
    thirdLevelPort = cfg.get<unsigned short>("thirdLevel.port");
    thirdLevelSSHCommand = cfg.get<std::string>("thirdLevel.sshCommand", "");
    extraSparseGrid = cfg.get<bool>("thirdLevel.extraSparseGrid");
    brokerOnSameSystem = static_cast<bool>(cfg.get_child_optional("thirdLevel.brokerOnSameSystem"));
    bool hasFractions = static_cast<bool>(cfg.get_child_optional("thirdLevel.fractionsOfScheme"));
    if (hasFractions) {
      std::string fractionsString = cfg.get<std::string>("thirdLevel.fractionsOfScheme");
      std::vector<std::string> stringVector;
      size_t pos = 0;
      std::string delimiter = " ";
      while ((pos = fractionsString.find(delimiter)) != std::string::npos) {
        stringVector.push_back(fractionsString.substr(0, pos));
        fractionsString.erase(0, pos + delimiter.length());
      }
      if (fractionsString.length() > 0) {
        stringVector.push_back(fractionsString);
      }
      fractionsOfScheme.resize(stringVector.size());
      std::transform(stringVector.begin(), stringVector.end(), fractionsOfScheme.begin(),
                     [](const std::string& val) { return std::stod(val); });
    } else {
      fractionsOfScheme = std::vector<real>(numSystems, 1. / static_cast<real>(numSystems));
    }
  }

  if (brokerOnSameSystem) {
    std::cout << "broker running on same system" << std::endl;
    // split communicator into the one for broker and one for workers + manager
    int globalID;
    MPI_Comm_rank(MPI_COMM_WORLD, &globalID);
    MPI_Comm worldComm;
    int color = 0;
    MPI_Comm_split(MPI_COMM_WORLD, color, globalID, &worldComm);
    theMPISystem()->initWorldReusable(worldComm, ngroup, nprocs);
    if (thirdLevelHost != "localhost") {
      throw std::runtime_error("Broker on same system, but third level host is not localhost");
    }
  } else {
    // divide the MPI processes into process group and initialize the
    // corresponding communicators
    theMPISystem()->init(ngroup, nprocs);
  }

  // todo: read in boundary vector from ctparam
  std::vector<BoundaryType> boundary(dim, 1);
  auto forwardDecomposition = false;

  // check whether parallelization vector p agrees with nprocs
  int checkProcs = 1;
  for (auto k : p) checkProcs *= k;
  if(checkProcs != IndexType(nprocs)){
    throw std::invalid_argument("process group size and parallelization do not match");
  }

  std::vector<LevelVector> levels;
  std::vector<combigrid::real> coeffs;
  std::vector<size_t> taskNumbers; // only used in case of static task assignment
  bool useStaticTaskAssignment = false;
  if (ctschemeFile == "") {
    /* generate a list of levelvectors and coefficients
     * CombiMinMaxScheme will create a classical combination scheme.
     * however, you could also read in a list of levelvectors and coefficients
     * from a file */
    CombiMinMaxScheme combischeme(dim, lmin, lmax);
    combischeme.createAdaptiveCombischeme();
    std::vector<LevelVector> fullLevels = combischeme.getCombiSpaces();
    std::vector<combigrid::real> fullCoeffs = combischeme.getCoeffs();

    // split scheme and assign each fraction to a system
    CombiThirdLevelScheme::createThirdLevelScheme(fullLevels, fullCoeffs, systemNumber,
                                                  numSystems, levels, coeffs, fractionsOfScheme);
    WORLD_MANAGER_EXCLUSIVE_SECTION {
      std::cout << fullLevels.size()
                << " component grids in full combination scheme; this system will run "
                << levels.size() << " of them." << std::endl;
      printCombiDegreesOfFreedom(levels, boundary);
    }
  } else {
    // read in CT scheme, if applicable
    std::unique_ptr<CombiMinMaxSchemeFromFile> scheme(
        new CombiMinMaxSchemeFromFile(dim, lmin, lmax, ctschemeFile));
    const auto& pgNumbers = scheme->getProcessGroupNumbers();
    if (pgNumbers.size() > 0) {
      useStaticTaskAssignment = true;
      const auto& allCoeffs = scheme->getCoeffs();
      const auto& allLevels = scheme->getCombiSpaces();
      const auto [itMin, itMax] = std::minmax_element(pgNumbers.begin(), pgNumbers.end());
      assert(*itMin == 0);  // make sure it starts with 0
      // assert(*itMax == ngroup - 1); // and goes up to the maximum group //TODO
      // filter out only those tasks that belong to "our" process group
      const auto& pgroupNumber = theMPISystem()->getProcessGroupNumber();
      for (size_t taskNo = 0; taskNo < pgNumbers.size(); ++taskNo) {
        if (pgNumbers[taskNo] == pgroupNumber) {
          taskNumbers.push_back(taskNo);
          coeffs.push_back(allCoeffs[taskNo]);
          levels.push_back(allLevels[taskNo]);
        }
      }
      MASTER_EXCLUSIVE_SECTION {
        std::cout << " Process group " << pgroupNumber << " will run " << levels.size() << " of "
                  << pgNumbers.size() << " tasks." << std::endl;
        printCombiDegreesOfFreedom(levels, boundary);
      }
    } else {
      // levels and coeffs are only used in manager
      WORLD_MANAGER_EXCLUSIVE_SECTION {
        std::cout << scheme->getCombiSpaces().size() << " tasks to distribute." << std::endl;
      }
    }
    WORLD_MANAGER_EXCLUSIVE_SECTION{
      assert(coeffs.size() == 0);
      assert(levels.size() == 0);
      coeffs = scheme->getCoeffs();
      levels = scheme->getCombiSpaces();
    }
  }
  // create load model
  std::unique_ptr<LoadModel> loadmodel = std::unique_ptr<LoadModel>(new LinearLoadModel());
  // std::unique_ptr<LoadModel> loadmodel = std::unique_ptr<LoadModel>(new AnisotropyLoadModel());
  // std::unique_ptr<LoadModel> loadmodel = std::unique_ptr<LoadModel>(new
  // LearningLoadModel(levels));

  // this code is only executed by the manager process
  WORLD_MANAGER_EXCLUSIVE_SECTION {
    // set up the ssh tunnel for third level communication, if necessary
    // todo: if this works, move to ProcessManager::setUpThirdLevel

    std::string hostnameInfo = "manager = " + boost::asio::ip::host_name();
    std::cout << hostnameInfo << std::endl;

    if (thirdLevelSSHCommand != "") {
      shellCommand::exec(thirdLevelSSHCommand.c_str());
      std::cout << thirdLevelSSHCommand << " returned " << std::endl;
    }

    // output combination scheme
    std::cout << "lmin = " << lmin << std::endl;
    std::cout << "lmax = " << lmax << std::endl;
    // std::cout << "CombiScheme: " << std::endl;
    // for (const LevelVector& level : levels) std::cout << level << std::endl;

    // create Tasks
    TaskContainer tasks;
    std::vector<size_t> taskIDs;
    if (!useStaticTaskAssignment) {
      // the world manager only needs to have task instances if it needs to distribute tasks
      tasks.reserve(levels.size());
      taskIDs.reserve(levels.size());
      for (size_t i = 0; i < levels.size(); i++) {
        Task* t = new TaskAdvection(levels[i], boundary, coeffs[i], loadmodel.get(), dt, nsteps, p);
        // Task* t = new TaskConstParaboloid(levels[i], boundary, coeffs[i], loadmodel);
        // Task* t = new TaskCount(dim, levels[i], boundary, coeffs[i], loadmodel.get());

        static_assert(!isGENE, "isGENE");

        tasks.push_back(t);
        taskIDs.push_back(t->getID());
      }
    } else {
      // taskIDs.resize(levels.size());
      // std::iota(taskIDs.begin(), taskIDs.end(), 0);
    }

    // create combiparameters
    auto reduceCombinationDimsLmax = LevelVector(dim, 1);
    // lie about ncombi, because default is to not use reduced dims for last combi step,
    // which we don't want here because it makes the sparse grid too large
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs, ncombi * 2, 1,
                           CombinationVariant::sparseGridReduce, p, LevelVector(dim, 0),
                           reduceCombinationDimsLmax, chunkSizeInMebibyte, forwardDecomposition,
                           thirdLevelHost, thirdLevelPort, 0);
    IndexVector minNumPoints(dim), maxNumPoints(dim);
    for (DimType d = 0; d < dim; ++d) {
      minNumPoints[d] = combigrid::getNumDofNodal(lmin[d], boundary[d]);
      maxNumPoints[d] = combigrid::getNumDofNodal(lmax[d], boundary[d]);
    }
    // first, test if decomposition possible for small resolution
    auto decomposition = combigrid::getDefaultDecomposition(minNumPoints, p, forwardDecomposition);
    // then assign the actual used one
    decomposition = combigrid::getDefaultDecomposition(maxNumPoints, p, forwardDecomposition);
    // default decomposition works only for powers of 2!
    params.setDecomposition(decomposition);
    std::cout << "manager: generated parameters" << std::endl;

    ProcessGroupManagerContainer pgroups;
    for (size_t i = 0; i < ngroup; ++i) {
      int pgroupRootID(i);
      pgroups.emplace_back(std::make_shared<ProcessGroupManager>(pgroupRootID));
    }
    // create abstraction for Manager
    ProcessManager manager(pgroups, tasks, params, std::move(loadmodel));
    manager.updateCombiParameters();
    auto durationParams = Stats::getDuration("manager update parameters")/ 1000.0;
    std::cout << "manager: updated parameters in " << durationParams << " seconds" << std::endl;

    /* distribute task according to load model and start computation for
     * the first time */
    Stats::startEvent("manager run first");
    if (useStaticTaskAssignment) {
      manager.runnext();
      std::cout << "manager: initialize sparse grid data structures" << std::endl;
      manager.initDsgus();
    } else {
      manager.runfirst(true);
    }
    Stats::stopEvent("manager run first");
    auto durationInit = Stats::getDuration("manager init dsgus")/ 1000.0;
    auto durationRun = Stats::getDuration("manager run first")/ 1000.0;
    std::cout << "manager: ran solver in " << durationRun << " seconds, of which SG init were " << durationInit << "" << std::endl;

    // exchange subspace sizes to unify the dsgs in the third level case
    if (hasThirdLevel) {
      Stats::startEvent("manager unify subspace sizes with remote");
      std::cout << "manager: unify sparse grid data structures w/ remote" << std::endl;
      manager.unifySubspaceSizesThirdLevel(extraSparseGrid);
      Stats::stopEvent("manager unify subspace sizes with remote");
      auto durationUnify = Stats::getDuration("manager unify subspace sizes with remote") / 1000.0;
      std::cout << "manager: unified SG in " << durationUnify << " seconds" << std::endl;
    } else {
      //Stats::startEvent("manager pretend unify subspace sizes with remote");
      //std::cout << "manager: unify sparse grid data structures w/ remote" << std::endl;
      //manager.pretendUnifySubspaceSizesThirdLevel();
      //Stats::stopEvent("manager pretend unify subspace sizes with remote");
      //auto durationUnify =
      //    Stats::getDuration("manager pretend unify subspace sizes with remote") / 1000.0;
      //std::cout << "manager: unified SG in " << durationUnify << " seconds" << std::endl;
    }

    for (size_t i = 1; i < ncombi; ++i) {

      Stats::startEvent("manager combine");
      if (hasThirdLevel) {
        manager.combineThirdLevel();
      } else {
        //manager.pretendCombineThirdLevelForWorkers();
        manager.combine();
      }
      // manager.waitAllFinished();
      Stats::stopEvent("manager combine");
      if (evalMCError && i % 10 == 0) {
        std::cout << "manager: eval Monte Carlo" << std::endl;
        managerMonteCarlo(manager, dim, static_cast<double>(i * nsteps) * dt, hasThirdLevel);
      }
      auto durationCombine = Stats::getDuration("manager combine")/ 1000.0;
      std::cout << "combination " << i << " took: " << durationCombine << " seconds" << std::endl;

      Stats::startEvent("manager write to disk");
      manager.writeDSGsToDisk("uftp_dsgu_");
      Stats::stopEvent("manager write to disk");
      auto durationWrite = Stats::getDuration("manager write to disk")/ 1000.0;
      std::cout << "write " << i << " took: " << durationWrite << " seconds" << std::endl;

      Stats::startEvent("manager read from disk");
      manager.readDSGsFromDisk("uftp_dsgu_");
      Stats::stopEvent("manager read from disk");
      auto durationRead = Stats::getDuration("manager read from disk")/ 1000.0;
      std::cout << "read " << i << " took: " << durationRead << " seconds" << std::endl;

      // run tasks for next time interval
      Stats::startEvent("manager run");
      manager.runnext();
      // manager.waitAllFinished();
      Stats::stopEvent("manager run");
      durationRun = Stats::getDuration("manager run")/ 1000.0;
      std::cout << "calculation " << i << " took: " << durationRun << " seconds" << std::endl;
    }

    Stats::startEvent("manager combine");
    if (hasThirdLevel) {
      manager.combineThirdLevel();
    } else {
      //manager.pretendCombineThirdLevelForWorkers();
      manager.combine();
    }
    Stats::stopEvent("manager combine");

    if (evalMCError) {
      std::cout << "manager: eval Monte Carlo" << std::endl;
      managerMonteCarlo(manager, dim, static_cast<double>(ncombi * nsteps) * dt, hasThirdLevel);
    }
    std::cout << "manager exit" << std::endl;

    // send exit signal to workers in order to enable a clean program termination
    manager.exit();
  }

  // this code is only executed by the worker processes
  else {
    // create abstraction of the process group from the worker's view
    ProcessGroupWorker pgroup;

    // wait for instructions from manager
    SignalType signal = -1;

    while (signal != EXIT) {
      signal = pgroup.wait();
      // if using static task assignment, we initialize all tasks after the combi parameters are
      // updated
      if (useStaticTaskAssignment) {
        if (signal == UPDATE_COMBI_PARAMETERS) {
          // initialize all "our" tasks
          for (size_t taskIndex = 0; taskIndex < taskNumbers.size(); ++taskIndex) {
            auto task = std::unique_ptr<Task>(new TaskAdvection(
                levels[taskIndex], boundary, coeffs[taskIndex], loadmodel.get(), dt, nsteps, p));
            task->setID(taskNumbers[taskIndex]);
            pgroup.initializeTask(std::move(task));
          }
        }
        // make sure not to initialize them twice
        if (signal == RUN_FIRST) {
          throw std::invalid_argument(
              "trying to initialize tasks from file AND manager distribution");
        }
      }
    }
  }

  Stats::finalize();

  /* write stats to json file for postprocessing */
  Stats::write("timers-" + std::to_string(systemNumber) + ".json");

  return 0;
}
