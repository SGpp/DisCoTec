
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <string>
#include <vector>
#include <algorithm>

#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"
#include "sgpp/distributedcombigrid/combischeme/CombiThirdLevelScheme.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LinearLoadModel.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupManager.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessManager.hpp"
#include "sgpp/distributedcombigrid/utils/MonteCarlo.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"

// this is not a boost test, because it is intended to be actually run on different machines
//

using namespace combigrid;

void mockUpDSGWriteToDisk(std::string filePrefix,
                          const std::vector<long long int>& dsgPartitionSizes) {
  for (size_t partitionIndex = 0; partitionIndex < dsgPartitionSizes.size(); ++partitionIndex) {
    std::string myFilename = filePrefix + std::to_string(partitionIndex);
    std::ofstream ofp(myFilename, std::ios::out | std::ios::binary);
    // generate data on heap
    std::unique_ptr<std::vector<real>> mockUpData(
        new std::vector<real>(dsgPartitionSizes[partitionIndex]));
    combigrid::montecarlo::getNumberSequenceFromSeed(*mockUpData, partitionIndex);
    ofp.write(reinterpret_cast<const char*>(mockUpData->data()),
              dsgPartitionSizes[partitionIndex] * sizeof(real));
    ofp.close();
  }
}

void validateExchangedData(std::string filePrefix,
                           const std::vector<long long int>& dsgPartitionSizes) {
  // generate data on heap
  std::unique_ptr<std::vector<real>> mockUpData(
      new std::vector<real>(dsgPartitionSizes[0]));
  combigrid::montecarlo::getNumberSequenceFromSeed(*mockUpData, 0);
  std::cout << "test first data" << std::endl;
  decltype(mockUpData) mockUpDataInverted(new std::vector<real>());
  mockUpDataInverted->reserve(mockUpData->size());
  std::cout << "test other data" << std::endl;
  std::transform(mockUpData->begin(), mockUpData->end(), std::back_inserter(*mockUpDataInverted),
                 std::negate<real>());
  std::cout << "test transform data" << std::endl;
  std::string myFilename = filePrefix + std::to_string(0);
  std::ifstream ifp(myFilename, std::ios::in | std::ios::binary);
  std::unique_ptr<std::vector<real>> readData(new std::vector<real>(dsgPartitionSizes[0]));
  ifp.read(reinterpret_cast<char*>(readData->data()), dsgPartitionSizes[0] * sizeof(real));
  for (size_t i = 0; i < mockUpData->size(); ++i) {
    if (readData->at(i) + mockUpDataInverted->at(i) != 0.0) {
      throw std::runtime_error("binary data does not match");
    }
  }
}

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

  /* when using timers (TIMING is defined in Stats), the Stats class must be
   * initialized at the beginning of the program. (and finalized in the end)
   */
  Stats::initialize();

  // read in parameter file
  std::string paramfile;
  if (argc > 1) {
    paramfile = argv[1];
  } else {
    throw std::runtime_error("pass parameter file as argument plz");
  }
  boost::property_tree::ptree cfg;
  boost::property_tree::ini_parser::read_ini(paramfile, cfg);

  // number of process groups and number of processes per group
  size_t ngroup = 0;
  size_t nprocs = 0;

  // divide the MPI processes into process group and initialize the
  // corresponding communicators
  theMPISystem()->init(ngroup, nprocs);

  /* read in parameters from ctparam */
  DimType dim = cfg.get<DimType>("ct.dim");
  LevelVector lmin(dim), lmax(dim);
  IndexVector p(dim);
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
  if (hasThirdLevel) {
    std::cout << "Using third-level parallelism" << std::endl;
    thirdLevelHost = cfg.get<std::string>("thirdLevel.host");
    systemNumber = cfg.get<unsigned int>("thirdLevel.systemNumber");
    numSystems = cfg.get<unsigned int>("thirdLevel.numSystems");
    thirdLevelPort = cfg.get<unsigned short>("thirdLevel.port");
    thirdLevelSSHCommand = cfg.get<std::string>("thirdLevel.sshCommand", "");
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

  std::vector<bool> boundary(dim, true);
  auto forwardDecomposition = true;

  std::vector<LevelVector> levels;
  std::vector<combigrid::real> coeffs;
  std::vector<size_t> taskNumbers; // only used in case of static task assignment
  bool useStaticTaskAssignment = false;
  long long int ctDOF = 0;
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
      ctDOF = printCombiDegreesOfFreedom(levels);
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
        ctDOF = printCombiDegreesOfFreedom(levels);
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
    WORLD_MANAGER_EXCLUSIVE_SECTION {
      ctDOF = printCombiDegreesOfFreedom(levels);
    }
  }
  // create load model
  std::unique_ptr<LoadModel> loadmodel = std::unique_ptr<LoadModel>(new LinearLoadModel());

  // this code is only executed by the manager process
  WORLD_MANAGER_EXCLUSIVE_SECTION {
    // set up the ssh tunnel for third level communication, if necessary
    // todo: if this works, move to ProcessManager::setUpThirdLevel
    // if (thirdLevelSSHCommand != "") {
    //   shellCommand::exec(thirdLevelSSHCommand.c_str());
    //   std::cout << thirdLevelSSHCommand << " returned " << std::endl;
    // }

    // output combination scheme
    std::cout << "lmin = " << lmin << std::endl;
    std::cout << "lmax = " << lmax << std::endl;
    // std::cout << "CombiScheme: " << std::endl;
    // for (const LevelVector& level : levels) std::cout << level << std::endl;

    // create combiparameters
    std::vector<size_t> taskIDs(coeffs.size());
    auto reduceCombinationDimsLmax = LevelVector(dim, 1);
    // lie about ncombi, because default is to not use reduced dims for last combi step,
    // which we don't want here because it makes the sparse grid too large
    // CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs, ncombi*2, 1, p,
    //                        LevelVector(dim, 0), reduceCombinationDimsLmax,
    //                        forwardDecomposition, thirdLevelHost, thirdLevelPort, 0);
    //TODO!

    CombiParameters params(dim, lmin, lmax, boundary, levels, coeffs, taskIDs, ncombi*2, 1, p,
                           LevelVector(dim, 0), reduceCombinationDimsLmax,
                           forwardDecomposition);
    // auto sgDOF = printSGDegreesOfFreedomAdaptive(lmin, lmax-reduceCombinationDimsLmax);
    // if(ctDOF < sgDOF) {
    //   throw std::runtime_error("Partitions don't match " + std::to_string(ctDOF) + " " + std::to_string(sgDOF));
    // }
    std::vector<IndexVector> decomposition;
    for (DimType d = 0; d < dim; ++d) {
      if (p[d] > (powerOfTwo[lmin[d]] + (boundary[d] ? +1 : -1))) {
        throw std::runtime_error(
            "change p! not all processes can have points on minimum level with current p.");
      }
      IndexVector di;
      if (p[d] == 1) {
        di = {0};
      } else if (p[d] == 2) {
        // forwardDecomposition for powers of 2!
        assert(forwardDecomposition && boundary[d]);
        di = {0, powerOfTwo[lmax[d]]/p[d] + 1};
      } else if (p[d] == 3 && lmax[d] == 10) {
        // naive
        di = {0, 342, 683};
      } else if (p[d] == 3 && lmax[d] == 18 && lmin[d] == 1) {
        // // naive
        // di = {0, 87382, 174763};
        // optimal
        di = {0, 80450, 181695};
      } else if (p[d] == 5 && lmax[d] == 19 && lmin[d] == 2) {
        //naive
        // di = {0, 104858, 209716, 314573, 419431};
        //optimal
        di = {0, 98304, 196609, 327680, 425985};
      } else {
        throw std::runtime_error("please provide a decomposition matching p and lmax");
      }
      decomposition.push_back(di);
    }
    params.setDecomposition(decomposition);

     // create abstraction for Manager
    ProcessGroupManagerContainer pgroups;
    TaskContainer tasks;
    ProcessManager manager(pgroups, tasks, params, std::move(loadmodel));

    manager.updateCombiParameters();

    // auto dsguSizes =
    //     getPartitionedNumDOFSGAdaptive(lmin, lmax - reduceCombinationDimsLmax, lmax, decomposition);
    // auto sgSumDOFPartitioned = std::accumulate(dsguSizes.begin(), dsguSizes.end(), static_cast<size_t>(0));
    // if(sgSumDOFPartitioned != sgDOF || ctDOF < sgDOF) {
    //   throw std::runtime_error("Partitions don't match" + std::to_string(sgSumDOFPartitioned) + " " + std::to_string(sgDOF));
    // }
    // Stats::startEvent("manager calculate dsguConjointSizes");
    // std::vector<long long int> dsguConjointSizes;
    // {
    //   std::unique_ptr<CombiMinMaxSchemeFromFile> scheme(
    //       new CombiMinMaxSchemeFromFile(dim, lmin, lmax, ctschemeFile));
    //   dsguConjointSizes = getPartitionedNumDOFSGConjoint(*scheme, lmin, lmax, decomposition);
    // }
    // Stats::stopEvent("manager calculate dsguConjointSizes");
    Stats::startEvent("manager calculate dsguConjointSize");
    long long int dsguConjointSize;
    {
      std::cout <<"scheme from file" << std::endl;
      std::unique_ptr<CombiMinMaxSchemeFromFile> scheme(
          new CombiMinMaxSchemeFromFile(dim, lmin, lmax, ctschemeFile));
      std::cout <<"calc conjoint" << std::endl;
      dsguConjointSize = getNumDOFSGConjoint(*scheme, lmin);
      std::cout <<"have dof conjoint" << std::endl;
    }
    Stats::stopEvent("manager calculate dsguConjointSize");

    // std::cout << dsguSizes << std::endl;
    std::cout << "vs conjoint" << std::endl;
    std::cout << dsguConjointSize << std::endl;
    std::cout << "and that makes a total of DOF " << std::endl;
    // std::cout << std::accumulate(dsguSizes.begin(), dsguSizes.end(), 0ll) << std::endl;
    // std::cout << std::accumulate(dsguConjointSizes.begin(), dsguConjointSizes.end(), 0ll) << std::endl;
    // std::cout << "distributed over a total of partitions " << std::endl;
    // std::cout << dsguConjointSizes.size() << std::endl;

    for (size_t i = 1; i < ncombi; ++i) {
      Stats::startEvent("manager pretend to third-level combine");
      // if (hasThirdLevel) {
        // manager.pretendCombineThirdLevel(dsguConjointSizes, false);//TODO
        std::cout << "mock" << std::endl;
        mockUpDSGWriteToDisk("dsg_part_" + std::to_string(systemNumber) + "_", {dsguConjointSize});
        std::cout << "test" << std::endl;
        validateExchangedData("dsg_part_" + std::to_string(systemNumber) + "_", {dsguConjointSize});
        sleep(2);
      // } else {
      //   throw std::runtime_error("this was not intended!");
      // }
      Stats::stopEvent("manager pretend to third-level combine");
    }
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
