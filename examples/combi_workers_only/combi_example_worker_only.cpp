#include <algorithm>
#include <boost/serialization/export.hpp>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

// include user specific task. this is the interface to your application
#include "../distributed_third_level/TaskAdvection.hpp"
#include "combischeme/CombiMinMaxScheme.hpp"
#include "io/BroadcastParameters.hpp"
#include "io/H5InputOutput.hpp"
#include "loadmodel/LinearLoadModel.hpp"
#include "manager/CombiParameters.hpp"
#include "manager/ProcessGroupWorker.hpp"
#include "mpi/MPISystem.hpp"
#include "task/Task.hpp"
#include "utils/MonteCarlo.hpp"
#include "utils/Types.hpp"

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
  {
    // number of process groups and number of processes per group
    size_t ngroup = cfg.get<size_t>("manager.ngroup");
    size_t nprocs = cfg.get<size_t>("manager.nprocs");

    // divide the MPI processes into process group and initialize the
    // corresponding communicators
    theMPISystem()->initWorldReusable(MPI_COMM_WORLD, ngroup, nprocs, false, true);
    MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << getTimeStamp() << "initialized communicators"
                                               << std::endl;

    /* read other parameters from ctparam */
    DimType dim = cfg.get<DimType>("ct.dim");
    LevelVector lmin(dim), lmax(dim);
    std::vector<bool> hierarchizationDims(dim, true);
    std::vector<int> p(dim);
    cfg.get<std::string>("ct.lmin") >> lmin;
    cfg.get<std::string>("ct.lmax") >> lmax;
    cfg.get<std::string>("ct.p") >> p;
    if (cfg.get_child_optional("ct.hierarchization_dims")) {
      cfg.get<std::string>("ct.hierarchization_dims") >>
          hierarchizationDims;  // which dimensions should be hierarchized
    }

    size_t ncombi = cfg.get<size_t>("ct.ncombi");
    uint32_t chunkSizeInMebibyte = cfg.get<uint32_t>("ct.chunkSize", 128);
    std::string basis = cfg.get<std::string>("ct.basis", "hat_periodic");
    std::string ctschemeFile = cfg.get<std::string>("ct.ctscheme");
    combigrid::real dt = cfg.get<combigrid::real>("application.dt");
    size_t nsteps = cfg.get<size_t>("application.nsteps");
    bool evalMCError = cfg.get<bool>("application.mcerror", false);

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

    const auto& pgroupNumber = theMPISystem()->getProcessGroupNumber();
    // read in CT scheme
    {
      std::unique_ptr<CombiMinMaxSchemeFromFile> scheme(
          new CombiMinMaxSchemeFromFile(dim, lmin, lmax, ctschemeFile));
      size_t totalNumTasks = 0;
      if (scheme->getProcessGroupNumbers().size() > 0) {
        totalNumTasks =
            combigrid::getAssignedLevels(*scheme, pgroupNumber, levels, coeffs, taskNumbers);
      } else {
        totalNumTasks =
            combigrid::getLoadBalancedLevels(*scheme, pgroupNumber, theMPISystem()->getNumGroups(),
                                             boundary, levels, coeffs, taskNumbers);
      }
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
    std::unique_ptr<LoadModel> loadmodel = std::unique_ptr<LoadModel>(new LinearLoadModel());

    // create combiparameters
    auto reduceCombinationDimsLmax = LevelVector(dim, 1);
    CombiParameters params(dim, lmin, lmax, boundary, ncombi, 1,
                           CombinationVariant::chunkedOutgroupSparseGridReduce, p,
                           LevelVector(dim, 0), reduceCombinationDimsLmax, chunkSizeInMebibyte,
                           forwardDecomposition);
    setCombiParametersHierarchicalBasesUniform(params, basis);
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
    if (evalMCError) {
      interpolationCoords.resize(1e5, std::vector<double>(dim, -1.));
      std::string interpolationCoordsFile = "interpolation_coords_" + std::to_string(dim) + "D_" +
                                            std::to_string(interpolationCoords.size()) + ".h5";
      // if the file does not exist, one rank creates it
      if (theMPISystem()->getWorldRank() == 0) {
        if (!std::filesystem::exists(interpolationCoordsFile)) {
          interpolationCoords = montecarlo::getRandomCoordinates(interpolationCoords.size(), dim);
          h5io::writeValuesToH5File(interpolationCoords, interpolationCoordsFile, "worker_group",
                                    "only");
        }
      }
      interpolationCoords = broadcastParameters::getCoordinatesFromRankZero(
          interpolationCoordsFile, theMPISystem()->getWorldComm());

      if (interpolationCoords.size() != 1e5) {
        sleep(1);
        throw std::runtime_error("not enough interpolation coordinates");
      }
    }

    ProcessGroupWorker worker;
    worker.setCombiParameters(std::move(params));

    // create Tasks
    worker.initializeAllTasks<TaskAdvection>(levels, coeffs, taskNumbers, loadmodel.get(), dt,
                                             nsteps, p);
    MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << getTimeStamp() << "worker: initialized tasks "
                                               << std::endl;

    worker.initCombinedDSGVector();
    MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << getTimeStamp() << "worker: initialized SG"
                                               << std::endl;

    MASTER_EXCLUSIVE_SECTION {
      std::cout << getTimeStamp() << "group " << theMPISystem()->getProcessGroupNumber()
                << ": set sparse grid sizes, will allocate "
                << static_cast<real>(worker.getCombinedDSGVector()[0]->getAccumulatedDataSize() *
                                     sizeof(CombiDataType)) /
                       1e6
                << " MB (but only "
                << static_cast<real>(combigrid::CombiCom::getGlobalReduceChunkSize<CombiDataType>(
                                         chunkSizeInMebibyte) *
                                     sizeof(CombiDataType)) /
                       1e6
                << " MB at once)" << std::endl;
    }

    // allocate sparse grids
    worker.zeroDsgsData();

    MPI_Barrier(theMPISystem()->getWorldComm());

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
      MPI_Barrier(theMPISystem()->getWorldComm());
      worker.runAllTasks();
      auto durationRun = Stats::getDuration("run") / 1000.0;
      MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << getTimeStamp() << "calculation " << i
                                                 << " took: " << durationRun << " seconds"
                                                 << std::endl;

      if (evalMCError) {
        Stats::startEvent("write interpolated");
        // evaluate these errors by calling
        // `python3 tools/hdf5_interpolation_norms.py worker_interpolated_values_0.h5
        // --solution=advection --coordinates=interpolation_coords_6D_100000.h5`
        worker.writeInterpolatedValuesSingleFile(interpolationCoords, "worker_interpolated");
        Stats::stopEvent("write interpolated");
        OTHER_OUTPUT_GROUP_EXCLUSIVE_SECTION {
          MASTER_EXCLUSIVE_SECTION {
            std::cout << getTimeStamp() << "interpolation " << i
                      << " took: " << Stats::getDuration("write interpolated") / 1000.0
                      << " seconds" << std::endl;
          }
        }
      }

      MPI_Barrier(theMPISystem()->getWorldComm());
      auto startCombine = std::chrono::high_resolution_clock::now();
      worker.combineAtOnce();
      auto endCombine = std::chrono::high_resolution_clock::now();
      auto durationCombine =
          std::chrono::duration_cast<std::chrono::milliseconds>(endCombine - startCombine).count() /
          1000.0;
      MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << getTimeStamp() << "combination " << i
                                                 << " took: " << durationCombine << " seconds"
                                                 << std::endl;
    }

    // run tasks for last time interval
    worker.runAllTasks();
    auto durationRun = Stats::getDuration("run") / 1000.0;
    MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << getTimeStamp() << "last calculation " << ncombi
                                               << " took: " << durationRun << " seconds"
                                               << std::endl;
    if (evalMCError) {
      Stats::startEvent("write interpolated");
      worker.writeInterpolatedValuesSingleFile(interpolationCoords, "worker_interpolated");
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
  }

  Stats::finalize();

  /* write stats to json file for postprocessing */
  Stats::write("timers.json");

  return 0;
}
