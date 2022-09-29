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
#include "sgpp/distributedcombigrid/mpi/MPIMemory.hpp"
// include user specific task. this is the interface to your application

// to allow using test tasks
#define BOOST_CHECK

#include "TaskAdvection.hpp"

using namespace combigrid;

// this is necessary for correct function of task serialization
#include "sgpp/distributedcombigrid/utils/BoostExports.hpp"
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
  if (hasThirdLevel) {
    std::cout << "Using third-level parallelism" << std::endl;
    thirdLevelHost = cfg.get<std::string>("thirdLevel.host");
    systemNumber = cfg.get<unsigned int>("thirdLevel.systemNumber");
    numSystems = cfg.get<unsigned int>("thirdLevel.numSystems");
    thirdLevelPort = cfg.get<unsigned short>("thirdLevel.port");
    thirdLevelSSHCommand = cfg.get<std::string>("thirdLevel.sshCommand", "");
    extraSparseGrid = cfg.get<bool>("thirdLevel.extraSparseGrid");
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

  // todo: read in boundary vector from ctparam
  std::vector<bool> boundary(dim, true);
  auto forwardDecomposition = true;

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
    CombiThirdLevelScheme::createThirdLevelScheme(fullLevels, fullCoeffs, boundary, systemNumber,
                                                  numSystems, levels, coeffs, fractionsOfScheme);
    WORLD_MANAGER_EXCLUSIVE_SECTION {
      std::cout << fullLevels.size()
                << " component grids in full combination scheme; this system will run "
                << levels.size() << " of them." << std::endl;
      printCombiDegreesOfFreedom(levels);
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
        printCombiDegreesOfFreedom(levels);
      }
      WORLD_MANAGER_EXCLUSIVE_SECTION{
        coeffs = scheme->getCoeffs();
        levels = scheme->getCombiSpaces();
      }
    } else {
      // levels and coeffs are only used in manager
      WORLD_MANAGER_EXCLUSIVE_SECTION {
        coeffs = scheme->getCoeffs();
        levels = scheme->getCombiSpaces();
        std::cout << levels.size() << " tasks to distribute." << std::endl;
      }
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

    /* create an abstraction of the process groups for the manager's view
     * a pgroup is identified by the ID in gcomm
     */
    ProcessGroupManagerContainer pgroups;
    for (size_t i = 0; i < ngroup; ++i) {
      int pgroupRootID(i);
      pgroups.emplace_back(std::make_shared<ProcessGroupManager>(pgroupRootID));
    }

    // output combination scheme
    std::cout << "lmin = " << lmin << std::endl;
    std::cout << "lmax = " << lmax << std::endl;
    // std::cout << "CombiScheme: " << std::endl;
    // for (const LevelVector& level : levels) std::cout << level << std::endl;

    // create Tasks
    TaskContainer tasks;
    tasks.reserve(levels.size());
    std::vector<size_t> taskIDs;
    taskIDs.reserve(levels.size());
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
    auto reduceCombinationDimsLmax = LevelVector(dim, 1);
    // lie about ncombi, because default is to not use reduced dims for last combi step,
    // which we don't want here because it makes the sparse grid too large
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs, ncombi*2, 1, p,
                           LevelVector(dim, 0), reduceCombinationDimsLmax,
                           forwardDecomposition, thirdLevelHost, thirdLevelPort, 0);
    std::vector<IndexVector> decomposition;
    for (DimType d = 0; d < dim; ++d) {
      if (p[d] > (powerOfTwo[lmin[d]] + (boundary[d] ? +1 : -1))) {
        throw std::runtime_error(
            "change p! not all processes can have points on minimum level with current p.");
      }
      IndexVector di;
      if (p[d] == 1) {
        di = {0};
      } else if (p[d] == 2 || p[d] == 4) {
        // forwardDecomposition for powers of 2!
        assert(forwardDecomposition && boundary[d]);
        di = {0, powerOfTwo[lmax[d]]/p[d] + 1};
      } else if (p[d] == 3 && lmax[d] == 10) {
        // naive
        di = {0, 342, 683};
      } else if (p[d] == 3 && lmax[d] == 18 && lmin[d] == 1) {
        // // naive
        // di = {0, 87382, 174763};
        // // optimal for [1]^6 -- [18]^6
        // di = {0, 80448, 181697};
        // optimal for [1,2,2,2,2,2] -- [18,19,19,19,19,19]
        di = {0,  78849, 183296};
      } else if (p[d] == 5 && lmax[d] == 19 && lmin[d] == 2) {
        //naive
        // di = {0, 104858, 209716, 314573, 419431};
        //optimal for [2]^6--[19]^6
        di = {0, 98304, 196609, 327680, 425985};
      } else if (p[d] == 3 && lmax[d] == 19 && lmin[d] == 2) {
        //optimal for [2]^6--[19]^6
        di = {0, 160737, 363552};
      } else {
        throw std::runtime_error("please implement a decomposition matching p and lmax");
      }
      decomposition.push_back(di);
    }
    params.setDecomposition(decomposition);

    if (useStaticTaskAssignment) {
      // read in CT scheme -- again!
      std::unique_ptr<CombiMinMaxSchemeFromFile> scheme(new CombiMinMaxSchemeFromFile(
          dim, lmin, lmax, ctschemeFile));
      const auto& pgNumbers = scheme->getProcessGroupNumbers();
      for (size_t taskNo = 0; taskNo < tasks.size(); ++taskNo) {
        if (pgNumbers[taskNo] < ngroup) { //TODO remove for production
          pgroups[pgNumbers[taskNo]]->storeTaskReference(tasks[taskNo]);
        }
      }
    }

    // create abstraction for Manager
    ProcessManager manager(pgroups, tasks, params, std::move(loadmodel));
    manager.updateCombiParameters();
    std::cout << "set up component grids and run until first combination point" << std::endl;

    /* distribute task according to load model and start computation for
     * the first time */
    Stats::startEvent("manager run first");
    if (useStaticTaskAssignment) {
      manager.runnext();
      manager.initDsgus();
    } else {
      manager.runfirst();
    }
    Stats::stopEvent("manager run first");

    // exchange subspace sizes to unify the dsgs in the third level case
    if (hasThirdLevel) {
      Stats::startEvent("manager unify subspace sizes with remote");
      manager.unifySubspaceSizesThirdLevel(extraSparseGrid);
      Stats::stopEvent("manager unify subspace sizes with remote");
    }

    double start, finish;

    start = MPI_Wtime();
    for (size_t i = 1; i < ncombi; ++i) {

      Stats::startEvent("manager combine");
      if (hasThirdLevel) {
        manager.combineThirdLevel();
      } else {
        manager.combine();
      }
      // manager.waitAllFinished();
      Stats::stopEvent("manager combine");
      if (evalMCError && i % 10 == 0) {
        managerMonteCarlo(manager, dim, static_cast<double>(i * nsteps) * dt, hasThirdLevel);
      }
      finish = MPI_Wtime();
      std::cout << "combination " << i << " took: " << finish - start << " seconds" << std::endl;
      start = finish;

      // run tasks for next time interval
      // start = MPI_Wtime();
      Stats::startEvent("manager run");
      manager.runnext();
      // manager.waitAllFinished();
      Stats::stopEvent("manager run");
      finish = MPI_Wtime();
      std::cout << "calculation " << i << " took: " << finish - start << " seconds" << std::endl;
      start = finish;
    }

    Stats::startEvent("manager combine");
    if (hasThirdLevel) {
      manager.combineThirdLevel();
    } else {
      manager.combine();
    }
    Stats::stopEvent("manager combine");

    // // evaluate solution and
    // // write solution to file
    // std::string filename("out/solution_" + std::to_string(ncombi) + ".raw");
    // Stats::startEvent("manager write solution");
    // manager.parallelEval(leval, filename, 0);
    // Stats::stopEvent("manager write solution");

    if (evalMCError) {
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

#ifdef DEBUG_OUTPUT
    mpimemory::print_memory_usage_local();
    mpimemory::print_memory_usage_local();
    mpimemory::print_memory_usage_local();
    unsigned long former_vmrss = 0;
    unsigned long former_vmsize = 0;
    unsigned long current_vmrss = 0;
    unsigned long current_vmsize = 0;
    double sumMeasured = 0.;
    double sumExpected = 0.;
#endif  // DEBUG_OUTPUT
    while (signal != EXIT) {
      signal = pgroup.wait();
#ifdef DEBUG_OUTPUT
      MASTER_EXCLUSIVE_SECTION { std::cout << "signal " << signal << std::endl; }
      mpimemory::print_memory_usage_local();
      mpimemory::get_memory_usage_local_kb(&former_vmrss, &former_vmsize);
#endif  // DEBUG_OUTPUT
      // if using static task assignment, we initialize all tasks after the combi parameters are
      // updated
      if (useStaticTaskAssignment) {
        if (signal == UPDATE_COMBI_PARAMETERS) {
          // initialize all "our" tasks
          for (size_t taskIndex = 0; taskIndex < taskNumbers.size(); ++taskIndex) {
            auto task = new TaskAdvection(dim, levels[taskIndex], boundary, coeffs[taskIndex],
                                          loadmodel.get(), dt, nsteps, p);
            if (taskIndex > 99) {
              break;
            }
            task->setID(taskNumbers[taskIndex]);
            pgroup.initializeTaskAndFaults(task);
#ifdef DEBUG_OUTPUT
            mpimemory::print_memory_usage_local();
            mpimemory::get_memory_usage_local_kb(&current_vmrss, &current_vmsize);
            auto measured = static_cast<double>(current_vmrss - former_vmrss) / 1e6;
            sumMeasured += measured;
            auto ctdof = getCombiDegreesOfFreedom(levels[taskIndex]);
            auto expected = ctdof * 2. * 8. / 1e9;
            sumExpected += expected;
            MASTER_EXCLUSIVE_SECTION {
              std::cout << "measured " << measured << " GB, expected " << expected << std::endl;
            }

            if (task->getDistributedFullGrid().getElementVector().max_size() <
                    task->getDistributedFullGrid().getNrLocalElements() ||
                task->getDistributedFullGrid().getElementVector().size() <
                    task->getDistributedFullGrid().getNrLocalElements()) {
              throw std::runtime_error(
                  "max_size went wrong! " +
                  std::to_string(task->getDistributedFullGrid().getNrLocalElements()));
            }
            former_vmrss = current_vmrss;
            former_vmsize = current_vmsize;
#endif  // DEBUG_OUTPUT
          }
#ifdef DEBUG_OUTPUT
          MASTER_EXCLUSIVE_SECTION {
            std::cout << "sumMeasured " << sumMeasured << " GB, sumExpected " << sumExpected
                      << std::endl;
          }
#endif  // DEBUG_OUTPUT
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

  MPI_Finalize();

  return 0;
}
