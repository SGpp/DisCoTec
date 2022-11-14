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
// include user specific task. this is the interface to your application

using namespace combigrid;

#include "ConjointLarge.hpp"
// // declared in ConjointLarge.hpp
// std::vector<long long int> dsguConjointSizes;

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

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  std::string hostnameInfo = "manager = " + boost::asio::ip::host_name();
  std::cout << hostnameInfo << std::endl;

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
  size_t ngroup = 0;
  size_t nprocs = 0;

  /* read in parameters from ctparam */
  DimType dim = cfg.get<DimType>("ct.dim");
  LevelVector lmin(dim), lmax(dim);
  std::vector<int> p(dim);
  cfg.get<std::string>("ct.lmin") >> lmin;
  cfg.get<std::string>("ct.lmax") >> lmax;
  cfg.get<std::string>("ct.p") >> p;
  size_t ncombi= cfg.get<size_t>("ct.ncombi");
  std::string ctschemeFile = cfg.get<std::string>("ct.ctscheme", "");

  // read in third level parameters if available
  std::string thirdLevelHost, thirdLevelSSHCommand = "";
  unsigned int systemNumber = 0, numSystems = 1;
  unsigned short thirdLevelPort = 0;
  bool hasThirdLevel = static_cast<bool>(cfg.get_child_optional("thirdLevel"));
  std::vector<real> fractionsOfScheme;
  bool brokerOnSameSystem = false;
  if (hasThirdLevel) {
    std::cout << "Using third-level parallelism" << std::endl;
    thirdLevelHost = cfg.get<std::string>("thirdLevel.host");
    systemNumber = cfg.get<unsigned int>("thirdLevel.systemNumber");
    numSystems = cfg.get<unsigned int>("thirdLevel.numSystems");
    thirdLevelPort = cfg.get<unsigned short>("thirdLevel.port");
    thirdLevelSSHCommand = cfg.get<std::string>("thirdLevel.sshCommand", "");
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

  std::vector<BoundaryType> boundary(dim, 1);
  auto forwardDecomposition = false;

  std::vector<LevelVector> levels;
  std::vector<combigrid::real> coeffs;
  std::vector<size_t> taskNumbers;  // only used in case of static task assignment
  bool useStaticTaskAssignment = false;
  // 6D -> use precomputed sizes
  if (dim < 6) {
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
        WORLD_MANAGER_EXCLUSIVE_SECTION {
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
      WORLD_MANAGER_EXCLUSIVE_SECTION { printCombiDegreesOfFreedom(levels); }
    }
  }
  // create load model
  std::unique_ptr<LoadModel> loadmodel = std::unique_ptr<LoadModel>(new LinearLoadModel());

  // this code is only executed by the manager process
  WORLD_MANAGER_EXCLUSIVE_SECTION {
    // set up the ssh tunnel for third level communication, if necessary
    // todo: if this works, move to ProcessManager::setUpThirdLevel
    if (thirdLevelSSHCommand != "") {
      shellCommand::exec(thirdLevelSSHCommand.c_str());
      std::cout << thirdLevelSSHCommand << " returned " << std::endl;
    }

    // output combination scheme
    std::cout << "lmin = " << lmin << std::endl;
    std::cout << "lmax = " << lmax << std::endl;
    std::cout << "CombiScheme size: " << levels.size() << std::endl;

    // create combiparameters
    std::vector<size_t> taskIDs(levels.size());
    std::iota(taskIDs.begin(), taskIDs.end(), 0);
    auto reduceCombinationDimsLmax = LevelVector(dim, 1);
    // lie about ncombi, because default is to not use reduced dims for last combi step,
    // which we don't want here because it makes the sparse grid too large
    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs, ncombi*2, 1, p,
                           LevelVector(dim, 0), reduceCombinationDimsLmax,
                           forwardDecomposition, thirdLevelHost, thirdLevelPort, 0);
    // auto sgDOF = printSGDegreesOfFreedomAdaptive(lmin, lmax-reduceCombinationDimsLmax);
    IndexVector minNumPoints(dim), maxNumPoints(dim);
    for (DimType d = 0; d < dim; ++d) {
      minNumPoints[d] = combigrid::getNumDofNodal(lmin[d], boundary[d]);
      maxNumPoints[d] = combigrid::getNumDofNodal(lmax[d], boundary[d]);
    }
    // first, test if decomposition possible for small resolution
    auto decomposition = combigrid::getDefaultDecomposition(minNumPoints, p, forwardDecomposition);
    // then assign the actual used one
    decomposition = combigrid::getDefaultDecomposition(maxNumPoints, p, forwardDecomposition);
    {
      // compute conjoint size, precomputing this for the large scheme
      // (put into separate header)
      //   std::unique_ptr<CombiMinMaxSchemeFromFile> scheme(
      //       new CombiMinMaxSchemeFromFile(dim, lmin, lmax, ctschemeFile));
      //   dsguConjointSizes = getPartitionedNumDOFSGConjoint(*scheme, lmin, lmax, decomposition);
      if (dim < 6) {
        CombiMinMaxScheme scheme(levels, coeffs);
        dsguConjointSizes = getPartitionedNumDOFSGConjoint(scheme, lmin, lmax, decomposition);
      }
    }

    // std::cout << "conjoint" << std::endl;
    // std::cout << dsguConjointSizes << std::endl;
    std::cout << "and that makes a total of DOF " << std::endl;
    std::cout << std::accumulate(dsguConjointSizes.begin(), dsguConjointSizes.end(), 0ll)
              << std::endl;
    std::cout << "distributed over a total of partitions " << std::endl;
    std::cout << dsguConjointSizes.size() << std::endl;

    // create abstraction for Manager
    ProcessGroupManagerContainer pgroups;
    TaskContainer tasks;
    ProcessManager manager(pgroups, tasks, params, std::move(loadmodel));

    manager.updateCombiParameters();
    for (size_t i = 0; i < ncombi; ++i) {
      Stats::startEvent("manager pretend to third-level combine");
      if (hasThirdLevel) {
        auto numErrors = manager.pretendCombineThirdLevelForBroker(dsguConjointSizes, true);
        if (numErrors > 0) {
          throw std::runtime_error("wrong numbers!");
        }
      } else {
        throw std::runtime_error("this was not intended!");
      }
      Stats::stopEvent("manager pretend to third-level combine");
    }

    // send exit signal to workers in order to enable a clean program termination
    manager.exit();
  } else {
    throw std::runtime_error("manager only test should run on one process only!");
  }

  Stats::finalize();

  /* write stats to json file for postprocessing */
  Stats::write("timers-" + std::to_string(systemNumber) + ".json");

  MPI_Finalize();

  return 0;
}
