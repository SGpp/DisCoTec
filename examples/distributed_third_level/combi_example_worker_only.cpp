// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <boost/asio.hpp>
#include <boost/serialization/export.hpp>
#include <filesystem>
#include <string>
#include <vector>

#include "combischeme/CombiMinMaxScheme.hpp"
#include "combischeme/CombiThirdLevelScheme.hpp"
#include "fault_tolerance/FaultCriterion.hpp"
#include "fault_tolerance/StaticFaults.hpp"
#include "fault_tolerance/WeibullFaults.hpp"
#include "io/BroadcastParameters.hpp"
#include "io/H5InputOutput.hpp"
#include "loadmodel/LinearLoadModel.hpp"
#include "manager/CombiParameters.hpp"
#include "manager/ProcessGroupWorker.hpp"
#include "task/Task.hpp"
#include "utils/MonteCarlo.hpp"
#include "utils/Types.hpp"

// to allow using test tasks
#define BOOST_CHECK

#include "TaskAdvection.hpp"

using namespace combigrid;

// this is necessary for correct function of task serialization
#include "utils/BoostExports.hpp"
BOOST_CLASS_EXPORT(TaskAdvection)

int main(int argc, char** argv) {
  [[maybe_unused]] auto mpiOnOff = MpiOnOff(&argc, &argv);
  /* when using timers (TIMING is defined in Stats), the Stats class must be
   * initialized at the beginning of the program. (and finalized in the end)
   */
  Stats::initialize();
  auto startInit = std::chrono::high_resolution_clock::now();

  // only one rank reads inputs and broadcasts to others
  std::string paramfile = "ctparam";
  if (argc > 1) paramfile = argv[1];
  boost::property_tree::ptree cfg =
      broadcastParameters::getParametersFromRankZero(paramfile, MPI_COMM_WORLD);

  // number of process groups and number of processes per group
  size_t ngroup = cfg.get<size_t>("manager.ngroup");
  size_t nprocs = cfg.get<size_t>("manager.nprocs");

  theMPISystem()->initWorldReusable(MPI_COMM_WORLD, ngroup, nprocs, false, true);
  {
    /* read other parameters from ctparam */
    DimType dim = cfg.get<DimType>("ct.dim");
    LevelVector lmin(dim), lmax(dim);
    std::vector<int> p(dim);
    combigrid::real dt;
    size_t nsteps, ncombi;
    cfg.get<std::string>("ct.lmin") >> lmin;
    cfg.get<std::string>("ct.lmax") >> lmax;
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
    if (hasThirdLevel) {
      systemNumber = cfg.get<unsigned int>("thirdLevel.systemNumber");
      numSystems = cfg.get<unsigned int>("thirdLevel.numSystems");
      assert(numSystems == 2);
      assert(systemNumber < numSystems);
      extraSparseGrid = cfg.get<bool>("thirdLevel.extraSparseGrid", true);
      MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << "running in file-based third level mode"
                                                 << std::endl;
    }

    // periodic boundary conditions
    std::vector<BoundaryType> boundary(dim, 1);
    auto forwardDecomposition = false;

    // check whether parallelization vector p agrees with nprocs
    int checkProcs = 1;
    for (auto k : p) checkProcs *= k;
    if (checkProcs != IndexType(nprocs)) {
      throw std::invalid_argument("process group size and parallelization do not match");
    }

    std::vector<LevelVector> levels;
    std::vector<combigrid::real> coeffs;
    std::vector<size_t> taskNumbers;  // only used in case of static task assignment
    bool useStaticTaskAssignment = false;
    if (ctschemeFile == "") {
      throw std::runtime_error("No CT scheme file specified");
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
          std::cout << getTimeStamp() << " Process group " << pgroupNumber << " will run "
                    << levels.size() << " of " << pgNumbers.size() << " tasks." << std::endl;
          printCombiDegreesOfFreedom(levels, boundary);
        }
      }
    }
    if (!useStaticTaskAssignment) {
      throw std::runtime_error("Dynamic task assignment not to be used here");
    }
    // create load model
    std::unique_ptr<LoadModel> loadmodel = std::unique_ptr<LoadModel>(new LinearLoadModel());

    // create combiparameters
    auto reduceCombinationDimsLmax = LevelVector(dim, 1);
    // lie about ncombi, because default is to not use reduced dims for last combi step,
    // which we don't want here because it makes the sparse grid too large
    CombiParameters params(dim, lmin, lmax, boundary, ncombi * 2, 1,
                           CombinationVariant::chunkedOutgroupSparseGridReduce, p,
                           LevelVector(dim, 0), reduceCombinationDimsLmax, chunkSizeInMebibyte,
                           forwardDecomposition, thirdLevelHost, thirdLevelPort, 0);
    setCombiParametersHierarchicalBasesUniform(params, "hat_periodic");
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
    MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << getTimeStamp() << "generated parameters"
                                               << std::endl;

    // read interpolation coordinates
    std::vector<std::vector<double>> interpolationCoords;
    interpolationCoords.resize(1e5, std::vector<double>(dim, -1.));
    std::string interpolationCoordsFile = "interpolation_coords_" + std::to_string(dim) + "D_" +
                                          std::to_string(interpolationCoords.size()) + ".h5";
    // // if the file does not exist, one rank creates it
    // if (theMPISystem()->getWorldRank() == 0) {
    //   if (!std::filesystem::exists(interpolationCoordsFile)) {
    //     interpolationCoords = montecarlo::getRandomCoordinates(interpolationCoords.size(), dim);
    //     h5io::writeValuesToH5File(interpolationCoords, interpolationCoordsFile, "worker_group",
    //                               "only");
    //   }
    // }
    // get them e.g. with `wget https://darus.uni-stuttgart.de/api/access/datafile/195524` (1e6)
    // or `wget https://darus.uni-stuttgart.de/api/access/datafile/195545` (1e5)
    interpolationCoords = broadcastParameters::getCoordinatesFromRankZero(
        interpolationCoordsFile, theMPISystem()->getWorldComm());

    if (interpolationCoords.size() != 1e5) {
      sleep(1);
      throw std::runtime_error("not enough interpolation coordinates");
    }
    MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << getTimeStamp() << "read interpolation coordinates"
                                               << std::endl;

    ProcessGroupWorker worker;
    worker.setCombiParameters(std::move(params));

    worker.initializeAllTasks<TaskAdvection>(levels, coeffs, taskNumbers, loadmodel.get(), dt,
                                             nsteps, p);
    MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << getTimeStamp() << "worker: initialized tasks "
                                               << std::endl;

    worker.initCombinedDSGVector();
    MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << getTimeStamp() << "worker: initialized SG"
                                               << std::endl;

    // read (extra) sparse grid sizes, as generated with subspace_writer
    // for target scenarios, consider `wget
    // https://darus.uni-stuttgart.de/api/access/datafile/195543` or similar
    if (extraSparseGrid) {
      std::string conjointSubspaceFileName =  // cf. subspace_writer.cpp
          ctschemeFile.substr(
              0, ctschemeFile.length() - std::string("_part0_00008groups.json").length()) +
          "conjoint.sizes";
      worker.reduceExtraSubspaceSizes(conjointSubspaceFileName, true);
    }

    OUTPUT_GROUP_EXCLUSIVE_SECTION {
      MASTER_EXCLUSIVE_SECTION {
        std::cout << getTimeStamp() << "worker: read sizes, will allocate "
                  << static_cast<real>(worker.getCombinedDSGVector()[0]->getAccumulatedDataSize() *
                                       sizeof(CombiDataType)) /
                         1e6
                  << " (but only "
                  << static_cast<real>(combigrid::CombiCom::getGlobalReduceChunkSize<CombiDataType>(
                                           chunkSizeInMebibyte) *
                                       sizeof(CombiDataType)) /
                         1e6
                  << " MB at once)"
                  << " plus "
                  << static_cast<real>(worker.getExtraDSGVector()[0]->getAccumulatedDataSize() *
                                       sizeof(CombiDataType)) /
                         1e6
                  << " MB" << std::endl;
      }
    }
    MPI_Barrier(theMPISystem()->getWorldComm());
    // allocate sparse grids now
    worker.zeroDsgsData();
    MIDDLE_PROCESS_EXCLUSIVE_SECTION {
      auto endInit = std::chrono::high_resolution_clock::now();
      auto durationInit =
          std::chrono::duration_cast<std::chrono::seconds>(endInit - startInit).count();
      std::cout << getTimeStamp() << "initialization took: " << durationInit << " seconds"
                << std::endl;
    }
    MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << getTimeStamp() << "start simulation loop"
                                               << std::endl;
    for (size_t i = 0; i < ncombi; ++i) {
      // run tasks for next time interval
      worker.runAllTasks();
      auto durationRun = Stats::getDuration("run") / 1000.0;
      MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << getTimeStamp() << "calculation " << i
                                                 << " took: " << durationRun << " seconds"
                                                 << std::endl;

      if (evalMCError) {
        Stats::startEvent("write interpolated");
        worker.writeInterpolatedValuesSingleFile(
            interpolationCoords, "worker_interpolated_" + std::to_string(systemNumber));
        Stats::stopEvent("write interpolated");
        OTHER_OUTPUT_GROUP_EXCLUSIVE_SECTION {
          MASTER_EXCLUSIVE_SECTION {
            std::cout << getTimeStamp() << "interpolation " << i
                      << " took: " << Stats::getDuration("write interpolated") / 1000.0
                      << " seconds" << std::endl;
          }
        }
      }

      auto startCombineWrite = std::chrono::high_resolution_clock::now();
      std::string writeSparseGridFile =
          "dsg_" + std::to_string(systemNumber) + "_step" + std::to_string(i);
      std::string writeSparseGridFileToken = writeSparseGridFile + "_token.txt";

      worker.combineSystemWideAndWrite(writeSparseGridFile, writeSparseGridFileToken);
      // everyone writes partial stats
      Stats::writePartial("stats_worker_" + std::to_string(systemNumber) + "_group" +
                              std::to_string(theMPISystem()->getProcessGroupNumber()) + ".json",
                          theMPISystem()->getLocalComm());

      MIDDLE_PROCESS_EXCLUSIVE_SECTION {
        auto endCombineWrite = std::chrono::high_resolution_clock::now();
        auto durationCombineWrite =
            std::chrono::duration_cast<std::chrono::seconds>(endCombineWrite - startCombineWrite)
                .count();
        std::cout << getTimeStamp() << "combination-local/write " << i
                  << " took: " << durationCombineWrite << " seconds" << std::endl;
      }
      auto startCombineRead = std::chrono::high_resolution_clock::now();
      std::string readSparseGridFile;
      if (hasThirdLevel) {
        readSparseGridFile =
            "dsg_" + std::to_string((systemNumber + 1) % 2) + "_step" + std::to_string(i);
        std::string readSparseGridFileToken = readSparseGridFile + "_token.txt";
        worker.combineReadDistributeSystemWide(readSparseGridFile, readSparseGridFileToken, false,
                                               false);

      } else {
        readSparseGridFile = writeSparseGridFile;
        worker.combineReadDistributeSystemWide(readSparseGridFile, writeSparseGridFileToken, true,
                                               false);
      }
      MIDDLE_PROCESS_EXCLUSIVE_SECTION {
        auto endCombineRead = std::chrono::high_resolution_clock::now();
        auto durationCombineRead =
            std::chrono::duration_cast<std::chrono::seconds>(endCombineRead - startCombineRead)
                .count();
        std::cout << getTimeStamp() << "combination-wait/read/reduce " << i
                  << " took: " << durationCombineRead << " seconds ; read " << readSparseGridFile
                  << std::endl;
      }
    }
    // run tasks for last time interval
    worker.runAllTasks();
    auto durationRun = Stats::getDuration("run") / 1000.0;
    MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << getTimeStamp() << "last calculation " << ncombi
                                               << " took: " << durationRun << " seconds"
                                               << std::endl;

    if (evalMCError) {
      Stats::startEvent("write interpolated");
      worker.writeInterpolatedValuesSingleFile(
          interpolationCoords, "worker_interpolated_" + std::to_string(systemNumber));
      Stats::stopEvent("write interpolated");
      OTHER_OUTPUT_GROUP_EXCLUSIVE_SECTION {
        MASTER_EXCLUSIVE_SECTION {
          std::cout << getTimeStamp() << "last interpolation " << ncombi
                    << " took: " << Stats::getDuration("write interpolated") / 1000.0 << " seconds"
                    << std::endl;
        }
      }
    }
    worker.exit();

    Stats::finalize();

    /* write stats to json file for postprocessing */
    Stats::write("timers_system" + std::to_string(systemNumber) + "_group" +
                     std::to_string(theMPISystem()->getProcessGroupNumber()) + ".json",
                 theMPISystem()->getLocalComm());
  }

  return 0;
}
