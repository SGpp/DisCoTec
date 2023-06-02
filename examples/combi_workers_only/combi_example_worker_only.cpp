// include user specific task. this is the interface to your application
#include <algorithm>
#include <boost/serialization/export.hpp>
#include <chrono>
#include <iostream>
#include <string>
#include <vector>

#include "../distributed_advection/TaskAdvection.hpp"
#include "combischeme/CombiMinMaxScheme.hpp"
#include "io/BroadcastParameters.hpp"
#include "loadmodel/LinearLoadModel.hpp"
#include "manager/CombiParameters.hpp"
#include "manager/ProcessGroupWorker.hpp"
#include "mpi/MPISystem.hpp"
#include "task/Task.hpp"
#include "utils/Types.hpp"

using namespace combigrid;

// this is necessary for correct function of task serialization
#include "utils/BoostExports.hpp"
BOOST_CLASS_EXPORT(TaskAdvection)

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

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
    theMPISystem()->initWorldReusable(MPI_COMM_WORLD, ngroup, nprocs, false);
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
    std::string basis = cfg.get<std::string>("ct.basis", "hat_periodic");
    std::string ctschemeFile = cfg.get<std::string>("ct.ctscheme");
    combigrid::real dt = cfg.get<combigrid::real>("application.dt");
    size_t nsteps = cfg.get<size_t>("application.nsteps");

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
      const auto& pgNumbers = scheme->getProcessGroupNumbers();
      const auto& allCoeffs = scheme->getCoeffs();
      const auto& allLevels = scheme->getCombiSpaces();
      assert(allLevels.size() == allCoeffs.size());
      if (pgNumbers.size() > 0) {
        const auto [itMin, itMax] = std::minmax_element(pgNumbers.begin(), pgNumbers.end());
        assert(*itMin == 0);  // make sure it starts with 0
        // filter out only those tasks that belong to "our" process group
        for (size_t taskNo = 0; taskNo < pgNumbers.size(); ++taskNo) {
          if (pgNumbers[taskNo] == pgroupNumber) {
            taskNumbers.push_back(taskNo);
            coeffs.push_back(allCoeffs[taskNo]);
            levels.push_back(allLevels[taskNo]);
          }
        }
      } else {
        // assume tasks are distributed round-robin (no load balancing!)
        assert(allLevels.size() >= theMPISystem()->getNumGroups());
        for (size_t taskNo = pgroupNumber; taskNo < allLevels.size();
             taskNo += theMPISystem()->getNumGroups()) {
          taskNumbers.push_back(taskNo);
          coeffs.push_back(allCoeffs[taskNo]);
          levels.push_back(allLevels[taskNo]);
        }
      }
      assert(!levels.empty());
      assert(levels.size() == coeffs.size());

      MASTER_EXCLUSIVE_SECTION {
        std::cout << getTimeStamp() << " Process group " << pgroupNumber << " will run "
                  << levels.size() << " of " << allLevels.size() << " tasks." << std::endl;
        printCombiDegreesOfFreedom(levels, boundary);
      }
    }

    // create load model
    std::unique_ptr<LoadModel> loadmodel = std::unique_ptr<LoadModel>(new LinearLoadModel());

    // create combiparameters
    auto reduceCombinationDimsLmax = LevelVector(dim, 1);
    CombiParameters params(dim, lmin, lmax, boundary, ncombi, 1,
                           CombinationVariant::sparseGridReduce, p, LevelVector(dim, 0),
                           reduceCombinationDimsLmax, forwardDecomposition);
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

    ProcessGroupWorker worker;
    worker.setCombiParameters(std::move(params));

    // create Tasks
    worker.initializeAllTasks<TaskAdvection>(levels, coeffs, taskNumbers, loadmodel.get(), dt,
                                             nsteps, p);

    worker.initCombinedDSGVector();

    auto durationInitSG = Stats::getDuration("init dsgus") / 1000.0;
    MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout
        << getTimeStamp() << "worker: initialized SG, registration was " << durationInitSG
        << " seconds" << std::endl;

    OUTPUT_GROUP_EXCLUSIVE_SECTION {
      MASTER_EXCLUSIVE_SECTION {
        std::cout << getTimeStamp() << "worker: set sparse grid sizes, will allocate "
                  << static_cast<real>(worker.getCombinedDSGVector()[0]->getAccumulatedDataSize() *
                                       sizeof(CombiDataType)) /
                         1e6
                  << " MB" << std::endl;
      }
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

      MPI_Barrier(theMPISystem()->getWorldComm());
      worker.combineUniform();
      auto durationCombine = Stats::getDuration("combine") / 1000.0;
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

    // numerical evaluation of Monte Carlo errors -> cf. examples/distributed_advection

    worker.exit();
  }

  Stats::finalize();

  /* write stats to json file for postprocessing */
  Stats::write("timers.json");

  MPI_Finalize();

  return 0;
}
