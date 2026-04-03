#include <algorithm>
#include <boost/serialization/export.hpp>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

// include user specific task. this is the interface to your application
#include "../distributed_third_level/TaskAdvection.hpp"
#include "discotec/io/H5InputOutput.hpp"
#include "discotec/io/ParameterIO.hpp"
#include "discotec/loadmodel/LinearLoadModel.hpp"
#include "discotec/manager/ProcessGroupWorker.hpp"
#include "discotec/utils/MonteCarlo.hpp"

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
  MIDDLE_PROCESS_EXCLUSIVE_SECTION
  std::cout << getTimeStamp() << "initialized communicators" << std::endl;

  // build CombiParameters and get task distribution
  auto [params, tasks] = buildCombiParametersFromConfig(cfg);

  // application-specific parameters
  combigrid::real dt = cfg.get<combigrid::real>("application.dt");
  size_t nsteps = cfg.get<size_t>("application.nsteps");
  size_t ncombi = cfg.get<size_t>("ct.ncombi");
  bool evalMCError = cfg.get<bool>("application.mcerror", false);
  DimType dim = cfg.get<DimType>("ct.dim");

  MASTER_EXCLUSIVE_SECTION {
    auto pgroupNumber = theMPISystem()->getProcessGroupNumber();
    std::cout << getTimeStamp() << " Process group " << pgroupNumber << " will run "
              << tasks.levels.size() << " of " << tasks.totalNumTasks << " tasks." << std::endl;
    printCombiDegreesOfFreedom(tasks.levels, params.getBoundary());
  }

  // set decomposition (needed for this task type)
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

  // read interpolation coordinates for MC error evaluation
  std::vector<std::vector<double>> interpolationCoords;
  if (evalMCError) {
    interpolationCoords.resize(1e5, std::vector<double>(dim, -1.));
    std::string interpolationCoordsFile = "interpolation_coords_" + std::to_string(dim) + "D_" +
                                          std::to_string(interpolationCoords.size()) + ".h5";
    if (theMPISystem()->getWorldRank() == 0) {
      if (!std::filesystem::exists(interpolationCoordsFile)) {
        interpolationCoords =
            montecarlo::getRandomCoordinates(static_cast<int>(interpolationCoords.size()), dim);
        h5io::writeValuesToH5File(interpolationCoords, interpolationCoordsFile, "worker_group",
                                  "only");
      }
    }
    interpolationCoords = broadcastParameters::getCoordinatesFromRankZero(
        interpolationCoordsFile, theMPISystem()->getWorldComm());
    if (interpolationCoords.size() != static_cast<size_t>(1e5)) {
      sleep(1);
      throw std::runtime_error("not enough interpolation coordinates");
    }
  }

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

  MASTER_EXCLUSIVE_SECTION {
    uint32_t chunkSizeInMebibyte = cfg.get<uint32_t>("ct.chunkSize", 128);
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
    auto durationRun = static_cast<double>(Stats::getDuration("run")) / 1000.0;
    MIDDLE_PROCESS_EXCLUSIVE_SECTION
    std::cout << getTimeStamp() << "calculation " << i << " took: " << durationRun << " seconds"
              << std::endl;

    if (evalMCError) {
      Stats::startEvent("write interpolated");
      worker.writeInterpolatedValuesSingleFile(interpolationCoords, "worker_interpolated");
      Stats::stopEvent("write interpolated");
      OTHER_OUTPUT_GROUP_EXCLUSIVE_SECTION {
        MASTER_EXCLUSIVE_SECTION {
          std::cout << getTimeStamp() << "interpolation " << i << " took: "
                    << static_cast<double>(Stats::getDuration("write interpolated")) / 1000.0
                    << " seconds" << std::endl;
        }
      }
    }

    MPI_Barrier(theMPISystem()->getWorldComm());
    auto startCombine = std::chrono::high_resolution_clock::now();
    worker.combineAtOnce();
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

  // run tasks for last time interval
  worker.runAllTasks();
  auto durationRun = static_cast<double>(Stats::getDuration("run")) / 1000.0;
  MIDDLE_PROCESS_EXCLUSIVE_SECTION
  std::cout << getTimeStamp() << "last calculation " << ncombi << " took: " << durationRun
            << " seconds" << std::endl;
  if (evalMCError) {
    Stats::startEvent("write interpolated");
    worker.writeInterpolatedValuesSingleFile(interpolationCoords, "worker_interpolated");
    Stats::stopEvent("write interpolated");
    OTHER_OUTPUT_GROUP_EXCLUSIVE_SECTION {
      MASTER_EXCLUSIVE_SECTION {
        std::cout << getTimeStamp() << "last interpolation " << ncombi << " took: "
                  << static_cast<double>(Stats::getDuration("write interpolated")) / 1000.0
                  << " seconds" << std::endl;
      }
    }
  }

  worker.exit();

  Stats::finalize();
  Stats::write("timers.json");
  return 0;
}
