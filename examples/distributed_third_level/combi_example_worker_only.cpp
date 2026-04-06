// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <boost/asio.hpp>
#include <boost/serialization/export.hpp>
#include <filesystem>
#include <string>
#include <vector>

#include "discotec/io/H5InputOutput.hpp"
#include "discotec/io/ParameterIO.hpp"
#include "discotec/loadmodel/LinearLoadModel.hpp"
#include "discotec/manager/ProcessGroupWorker.hpp"
#include "discotec/utils/MonteCarlo.hpp"

// to allow using test tasks
#define BOOST_CHECK

#include "TaskAdvection.hpp"

using namespace combigrid;

// this is necessary for correct function of task serialization
#include "discotec/utils/BoostExports.hpp"
BOOST_CLASS_EXPORT(TaskAdvection)

int main(int argc, char** argv) {
  [[maybe_unused]] auto mpiOnOff = MpiOnOff(&argc, &argv);
  Stats::initialize();
  auto startInit = std::chrono::high_resolution_clock::now();

  // read ctparam: rank 0 reads and broadcasts
  std::string paramfile = argc > 1 ? argv[1] : "ctparam";
  auto [ngroup, nprocs, cfg] = readParameterFile(paramfile, MPI_COMM_WORLD);

  theMPISystem()->initWorldReusable(MPI_COMM_WORLD, ngroup, nprocs, false, true);

  // build CombiParameters and get task distribution
  auto [params, tasks] = buildCombiParametersFromConfig(cfg);

  // application-specific parameters
  DimType dim = cfg.get<DimType>("ct.dim");
  combigrid::real dt = cfg.get<combigrid::real>("application.dt");
  size_t nsteps = cfg.get<size_t>("application.nsteps");
  size_t ncombi = cfg.get<size_t>("ct.ncombi");
  bool evalMCError = cfg.get<bool>("application.mcerror", false);
  uint16_t numberOfFileParts = cfg.get<uint16_t>("io.numberParts", 1);

  theMPISystem()->initOutputGroupComm(numberOfFileParts);

  // third-level parameters
  bool hasThirdLevel = static_cast<bool>(cfg.get_child_optional("thirdLevel"));
  unsigned int systemNumber = 0, numSystems = 1;
  bool extraSparseGrid = true;
  if (hasThirdLevel) {
    systemNumber = cfg.get<unsigned int>("thirdLevel.systemNumber");
    numSystems = cfg.get<unsigned int>("thirdLevel.numSystems");
    assert(numSystems > 1);
    assert(systemNumber < numSystems);
    extraSparseGrid = cfg.get<bool>("thirdLevel.extraSparseGrid", true);
    MIDDLE_PROCESS_EXCLUSIVE_SECTION
    std::cout << "running in file-based third level mode" << std::endl;
  }

  MASTER_EXCLUSIVE_SECTION {
    auto pgroupNumber = theMPISystem()->getProcessGroupNumber();
    std::cout << getTimeStamp() << " Process group " << pgroupNumber << " will run "
              << tasks.levels.size() << " of " << tasks.totalNumTasks << " tasks." << std::endl;
    printCombiDegreesOfFreedom(tasks.levels, params.getBoundary());
  }

  // set decomposition
  auto& p = params.getParallelization();
  bool forwardDecomposition = cfg.get<bool>("ct.forwardDecomposition", false);
  IndexVector maxNumPoints(dim);
  for (DimType d = 0; d < dim; ++d) {
    maxNumPoints[d] = getNumDofNodal(params.getLMax()[d], params.getBoundary()[d]);
  }
  auto decomposition = getDefaultDecomposition(maxNumPoints, p, forwardDecomposition);
  params.setDecomposition(decomposition);
  MIDDLE_PROCESS_EXCLUSIVE_SECTION
  std::cout << getTimeStamp() << "generated parameters" << std::endl;

  // read interpolation coordinates
  std::vector<std::vector<double>> interpolationCoords;
  interpolationCoords.resize(1e5, std::vector<double>(dim, -1.));
  std::string interpolationCoordsFile = "interpolation_coords_" + std::to_string(dim) + "D_" +
                                        std::to_string(interpolationCoords.size()) + ".h5";
#ifdef DISCOTEC_USE_HIGHFIVE
  if (theMPISystem()->getWorldRank() == 0) {
    if (!std::filesystem::exists(interpolationCoordsFile)) {
      interpolationCoords =
          montecarlo::getRandomCoordinates(static_cast<int>(interpolationCoords.size()), dim);
      // write to a temporary file and rename atomically to avoid races
      // when multiple systems start concurrently in the same directory
      std::string tmpFile = interpolationCoordsFile + ".tmp_sys" + std::to_string(systemNumber);
      h5io::writeValuesToH5File(interpolationCoords, tmpFile, "worker_group", "only");
      std::filesystem::rename(tmpFile, interpolationCoordsFile);
    }
  }
#endif
  interpolationCoords = broadcastParameters::getCoordinatesFromRankZero(
      interpolationCoordsFile, theMPISystem()->getWorldComm());
  if (interpolationCoords.size() != static_cast<size_t>(1e5)) {
    sleep(1);
    throw std::runtime_error("not enough interpolation coordinates");
  }
  MIDDLE_PROCESS_EXCLUSIVE_SECTION
  std::cout << getTimeStamp() << "read interpolation coordinates" << std::endl;

  // set up worker, tasks, and sparse grid
  std::unique_ptr<LoadModel> loadmodel = std::make_unique<LinearLoadModel>();
  ProcessGroupWorker<> worker;
  worker.setCombiParameters(std::move(params));
  worker.initializeAllTasks<TaskAdvection>(tasks.levels, tasks.coeffs, tasks.taskNumbers,
                                           loadmodel.get(), dt, nsteps, p);
  MIDDLE_PROCESS_EXCLUSIVE_SECTION
  std::cout << getTimeStamp() << "worker: initialized tasks" << std::endl;

  worker.initCombinedDSGVector();
  MIDDLE_PROCESS_EXCLUSIVE_SECTION
  std::cout << getTimeStamp() << "worker: initialized SG" << std::endl;

  // read extra sparse grid sizes for third-level
  std::string ctschemeFile = cfg.get<std::string>("ct.ctscheme", "");
  if (extraSparseGrid && !ctschemeFile.empty()) {
    std::string conjointSubspaceFileName =
        ctschemeFile.substr(
            0, ctschemeFile.length() - std::string("_part0_00008groups.json").length()) +
        "conjoint.sizes";
    worker.reduceExtraSubspaceSizes({conjointSubspaceFileName}, true);
  }

  MASTER_EXCLUSIVE_SECTION {
    uint32_t chunkSizeInMebibyte = cfg.get<uint32_t>("ct.chunkSize", 64);
    std::cout << getTimeStamp() << "group " << theMPISystem()->getProcessGroupNumber()
              << ": set sparse grid sizes, will allocate "
              << static_cast<real>(worker.getCombinedDSGVector()[0]->getAccumulatedDataSize() *
                                   sizeof(CombiDataType)) /
                     1e6
              << " MB (but only "
              << static_cast<real>(
                     CombiCom::getGlobalReduceChunkSize<CombiDataType>(chunkSizeInMebibyte) *
                     sizeof(CombiDataType)) /
                     1e6
              << " MB at once)" << std::endl;
  }

  worker.zeroDsgsData();
  MPI_Barrier(theMPISystem()->getWorldComm());

  MIDDLE_PROCESS_EXCLUSIVE_SECTION {
    auto endInit = std::chrono::high_resolution_clock::now();
    auto durationInit =
        std::chrono::duration_cast<std::chrono::seconds>(endInit - startInit).count();
    std::cout << getTimeStamp() << "initialization took: " << durationInit << " seconds"
              << std::endl;
  }

  // simulation loop
  MIDDLE_PROCESS_EXCLUSIVE_SECTION
  std::cout << getTimeStamp() << "start simulation loop" << std::endl;
  for (size_t i = 0; i < ncombi; ++i) {
    MPI_Barrier(theMPISystem()->getWorldComm());
    worker.runAllTasks();
    MIDDLE_PROCESS_EXCLUSIVE_SECTION
    std::cout << getTimeStamp() << "calculation " << i
              << " took: " << static_cast<double>(Stats::getDuration("run")) / 1000.0 << " seconds"
              << std::endl;

    if (evalMCError) {
      Stats::startEvent("write interpolated");
      worker.writeInterpolatedValuesSingleFile(interpolationCoords, "worker_interpolated");
      Stats::stopEvent("write interpolated");
    }

    MPI_Barrier(theMPISystem()->getWorldComm());
    auto startCombine = std::chrono::high_resolution_clock::now();

    if (hasThirdLevel) {
      // third-level file-based exchange with iteration-keyed filenames
      std::string iterStr = std::to_string(i);
      std::string writeSparseGridFile =
          "dsgu_" + std::to_string(systemNumber) + "_i" + iterStr + ".json";
      std::string writeSparseGridFileToken =
          "dsgu_" + std::to_string(systemNumber) + "_i" + iterStr + "_token";
      worker.combineSystemWideAndWrite(writeSparseGridFile, writeSparseGridFileToken);

      std::vector<std::string> readSparseGridFiles, readSparseGridFileTokens;
      for (unsigned int sys = 0; sys < numSystems; ++sys) {
        if (sys != systemNumber) {
          readSparseGridFiles.push_back("dsgu_" + std::to_string(sys) + "_i" + iterStr + ".json");
          readSparseGridFileTokens.push_back("dsgu_" + std::to_string(sys) + "_i" + iterStr +
                                             "_token");
        }
      }
      worker.combineReadDistributeSystemWide(readSparseGridFiles, readSparseGridFileTokens, false,
                                             true);
    } else {
      worker.combineAtOnce();
    }

    auto endCombine = std::chrono::high_resolution_clock::now();
    auto durationCombine =
        static_cast<double>(
            std::chrono::duration_cast<std::chrono::milliseconds>(endCombine - startCombine)
                .count()) /
        1000.0;
    MIDDLE_PROCESS_EXCLUSIVE_SECTION
    std::cout << getTimeStamp() << "combination " << i << " took: " << durationCombine << " seconds"
              << std::endl;
  }

  // last iteration
  worker.runAllTasks();
  MIDDLE_PROCESS_EXCLUSIVE_SECTION
  std::cout << getTimeStamp() << "last calculation " << ncombi
            << " took: " << static_cast<double>(Stats::getDuration("run")) / 1000.0 << " seconds"
            << std::endl;
  if (evalMCError) {
    Stats::startEvent("write interpolated");
    worker.writeInterpolatedValuesSingleFile(interpolationCoords, "worker_interpolated");
    Stats::stopEvent("write interpolated");
  }

  worker.exit();
  Stats::finalize();

  /* write stats to json file for postprocessing */
  Stats::write("timers_" + std::to_string(theMPISystem()->getProcessGroupNumber()) + ".json",
               theMPISystem()->getLocalComm());

  return 0;
}
