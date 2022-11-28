#include "sgpp/distributedcombigrid/manager/ProcessGroupWorker.hpp"

#include "boost/lexical_cast.hpp"

#include "sgpp/distributedcombigrid/combicom/CombiCom.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/hierarchization/DistributedHierarchization.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LearningLoadModel.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupSignals.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/mpi/MPIUtils.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/DistributedSparseGrid.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/DistributedSparseGridUniform.hpp"

#include <algorithm>
#include <iostream>
#include <string>
#ifdef HAVE_HIGHFIVE
#include <chrono>
#include <random>
// highfive is a C++ hdf5 wrapper, available in spack (-> configure with right boost and mpi versions)
#include <highfive/H5File.hpp>
#endif
#include "sgpp/distributedcombigrid/mpi_fault_simulator/MPI-FT.h"

namespace combigrid {

ProcessGroupWorker::ProcessGroupWorker()
    : currentTask_(nullptr),
      status_(PROCESS_GROUP_WAIT),
      combinedUniDSGVector_(0),
      combiParameters_(),
      combiParametersSet_(false),
      currentCombi_(0) {
  t_fault_ = -1;
  startTimeIteration_ = (std::chrono::high_resolution_clock::now());
  MASTER_EXCLUSIVE_SECTION {
    std::string fname = "out/all-betas-" + std::to_string(theMPISystem()->getGlobalRank()) + ".txt";
    // betasFile_ = std::ofstream( fname, std::ofstream::out );
  }
}

ProcessGroupWorker::~ProcessGroupWorker() {
  for (auto& task : tasks_) {
    delete task;
    task = nullptr;
  }
}

// Do useful things with the info about how long a task took.
// this gets called whenever a task was run, i.e., signals RUN_FIRST(once), RUN_NEXT(possibly multiple times),
// RECOMPUTE(possibly multiple times), and in ready(possibly multiple times)
void ProcessGroupWorker::processDuration(const Task& t, const Stats::Event e,
                                         unsigned int numProcs) {
  return;
  MASTER_EXCLUSIVE_SECTION {
    DurationInformation info = {t.getID(), Stats::getEventDurationInUsec(e),
	    			                    t.getCurrentTime(), t.getCurrentTimestep(),
				                        theMPISystem()->getWorldRank(), static_cast<unsigned int>(numProcs)};
    MPIUtils::sendClass(&info, theMPISystem()->getManagerRank(), theMPISystem()->getGlobalComm());
  }
}

SignalType ProcessGroupWorker::wait() {
  if (status_ == PROCESS_GROUP_FAIL) {  // in this case worker got reused
    status_ = PROCESS_GROUP_WAIT;
  }
  if (status_ != PROCESS_GROUP_WAIT) {
#ifdef DEBUG_OUTPUT
    int myRank = theMPISystem()->getWorldRank();
    std::cout << "status is " << status_ << "of rank " << myRank << "\n";
    std::cout << "executing next task\n";
#endif
    return RUN_NEXT;
  }
  SignalType signal = -1;

  MASTER_EXCLUSIVE_SECTION {
    // receive signal from manager
    MPI_Recv(&signal, 1, MPI_INT, theMPISystem()->getManagerRank(), TRANSFER_SIGNAL_TAG,
             theMPISystem()->getGlobalComm(), MPI_STATUS_IGNORE);
  }
  // distribute signal to other processes of pgroup
  MPI_Bcast(&signal, 1, MPI_INT, theMPISystem()->getMasterRank(), theMPISystem()->getLocalComm());
#ifdef DEBUG_OUTPUT
  std::cout << theMPISystem()->getWorldRank() << " waits for signal " << signal << " \n";
#endif
  // process signal
  switch (signal) {
    case RUN_FIRST: {
      receiveAndInitializeTaskAndFaults();
      status_ = PROCESS_GROUP_BUSY;

      // execute task
      Stats::startEvent("worker run first");
      currentTask_->run(theMPISystem()->getLocalComm());
      Stats::Event e = Stats::stopEvent("worker run first");
      // std::cout << "from runfirst ";
      processDuration(*currentTask_, e, theMPISystem()->getNumProcs());
    } break;
    case RUN_NEXT: {
      assert(tasks_.size() > 0);
      // // free space for computation
      // deleteDsgsData();

      // reset finished status of all tasks
      if (tasks_.size() != 0) {
        for (size_t i = 0; i < tasks_.size(); ++i) tasks_[i]->setFinished(false);

        status_ = PROCESS_GROUP_BUSY;

        // set currentTask
        currentTask_ = tasks_[0];

        // run first task
        // if isGENE, this is done in GENE's worker_routines.cpp
        if (!isGENE) {
          Stats::startEvent("worker run");
        }

        Stats::Event e = Stats::Event();
        currentTask_->run(theMPISystem()->getLocalComm());
        e.end = std::chrono::high_resolution_clock::now();
        // std::cout << "from runnext ";
        processDuration(*currentTask_, e, theMPISystem()->getNumProcs());

      } else {
        std::cout << "Possible error: No tasks! \n";
      }

    } break;
    case ADD_TASK: {  // add a new task to the process group
      // initalize task and set values to zero
      // the task will get the proper initial solution during the next combine
      // TODO test if this signal works in case of not-GENE
      receiveAndInitializeTaskAndFaults();

      currentTask_->setZero();

      currentTask_->setFinished(true);

      if (isGENE) {
        currentTask_->changeDir(theMPISystem()->getLocalComm());
      }
    } break;
    case RESET_TASKS: {  // deleta all tasks (used in process recovery)
      std::cout << "resetting tasks" << std::endl;

      deleteTasks();
      status_ = PROCESS_GROUP_BUSY;

    } break;
    case EVAL: {
      // receive x

      // loop over all tasks
      // t.eval(x)
    } break;
    case EXIT: {
      // write out tasks that were in use when the computation ended
      // (i.e. after fault tolerance or rescheduling changes)
      MASTER_EXCLUSIVE_SECTION {
        // serialize tasks as string
        std::stringstream tasksStream;
        for (const auto& t : tasks_) {
          tasksStream <<t->getID() << ": " << t->getCoefficient() << t->getLevelVector()  << "; ";
        }
        std::string tasksString = tasksStream.str();
        Stats::setAttribute("tasks: levels", tasksString);
      }
      if (isGENE) {
        if(chdir("../ginstance")){};
      }
      deleteTasks();
    } break;
    case SYNC_TASKS: {
      MASTER_EXCLUSIVE_SECTION {
        for (size_t i = 0; i < tasks_.size(); ++i) {
          Task::send(&tasks_[i], theMPISystem()->getManagerRank(), theMPISystem()->getGlobalComm());
        }
      }
    } break;
    case INIT_DSGUS: {
      Stats::startEvent("worker initialize dsgu");
      initCombinedUniDSGVector();
      Stats::stopEvent("worker initialize dsgu");

    } break;
    case COMBINE: {  // start combination
      Stats::startEvent("worker combine");
      combineUniform();
      currentCombi_++;
      Stats::stopEvent("worker combine");

    } break;
    case COMBINE_LOCAL_AND_GLOBAL: {
      Stats::startEvent("worker combine local");
      combineLocalAndGlobal();
      Stats::stopEvent("worker combine local");

    } break;
    case COMBINE_THIRD_LEVEL: {
      Stats::startEvent("worker combine third level");
      combineThirdLevel();
      currentCombi_++;
      Stats::stopEvent("worker combine third level");
    } break;
    case COMBINE_THIRD_LEVEL_FILE: {
      Stats::startEvent("worker combine third level file");
      combineThirdLevelFileBased();
      currentCombi_++;
      Stats::stopEvent("worker combine third level file");
    } break;
    case WAIT_FOR_TL_COMBI_RESULT: {
      Stats::startEvent("worker wait third level result");
      waitForThirdLevelCombiResult();
      currentCombi_++;
      Stats::stopEvent("worker wait third level result");

    } break;
    case REDUCE_SUBSPACE_SIZES_TL_AND_ALLOCATE_EXTRA_SG: {
      Stats::startEvent("worker unify extra third level");
      reduceSubspaceSizesThirdLevel(true);
      Stats::stopEvent("worker unify extra third level");

    } break;
    case REDUCE_SUBSPACE_SIZES_TL: {
      Stats::startEvent("worker unify sizes third level");
      reduceSubspaceSizesThirdLevel(false);
      Stats::stopEvent("worker unify sizes third level");

    } break;
    case WAIT_FOR_TL_SIZE_UPDATE: {
      Stats::startEvent("worker wait third level size");
      waitForThirdLevelSizeUpdate();
      Stats::stopEvent("worker wait third level size");

    } break;
    case WRITE_DFGS_TO_VTK: {
      Stats::startEvent("worker write vtk all tasks");
      writeVTKPlotFilesOfAllTasks();
      Stats::stopEvent("worker write vtk all tasks");
    } break;
    case WRITE_DSGS_TO_DISK: {
      Stats::startEvent("worker write to disk");
      std::string filenamePrefix;
      MASTER_EXCLUSIVE_SECTION {
        MPIUtils::receiveClass(&filenamePrefix, theMPISystem()->getManagerRank(),
                               theMPISystem()->getGlobalComm());
      }
      MPIUtils::broadcastClass(&filenamePrefix, theMPISystem()->getMasterRank(),
                               theMPISystem()->getLocalComm());
      writeDSGsToDisk(filenamePrefix);
      Stats::stopEvent("worker write to disk");
    } break;
    case READ_DSGS_FROM_DISK: {
      Stats::startEvent("worker read from disk");
      std::string filenamePrefix;
      MASTER_EXCLUSIVE_SECTION {
        MPIUtils::receiveClass(&filenamePrefix, theMPISystem()->getManagerRank(),
                               theMPISystem()->getGlobalComm());
      }
      MPIUtils::broadcastClass(&filenamePrefix, theMPISystem()->getMasterRank(),
                               theMPISystem()->getLocalComm());
      readDSGsFromDisk(filenamePrefix);
      Stats::stopEvent("worker read from disk");
    } break;
    case GRID_EVAL: {  // not supported anymore

      Stats::startEvent("worker eval");
      gridEval();
      Stats::stopEvent("worker eval");

      return signal;

    } break;
    case COMBINE_FG: {
      combineFG();

    } break;
    case UPDATE_COMBI_PARAMETERS: {  // update combiparameters (e.g. in case of faults -> FTCT)

      updateCombiParameters();

    } break;
    case RECOMPUTE: {  // recompute the received task (immediately computes tasks ->
                       // difference to ADD_TASK)
      receiveAndInitializeTaskAndFaults();
      currentTask_->setZero();

      // fill task with combisolution
      if (!isGENE) {
        fillDFGFromDSGU(currentTask_);
      }
      // execute task
      Stats::Event e = Stats::Event();
      currentTask_->run(theMPISystem()->getLocalComm());
      e.end = std::chrono::high_resolution_clock::now();
      // std::cout << "from recompute ";
      processDuration(*currentTask_, e, theMPISystem()->getNumProcs());

    } break;
    case RECOVER_COMM: {  // start recovery in case of faults
      theMPISystem()->recoverCommunicators(true);
      return signal;
    } break;
    case PARALLEL_EVAL: {  // output final grid
      Stats::startEvent("worker parallel eval");
      parallelEval();
      Stats::stopEvent("worker parallel eval");
    } break;
    case DO_DIAGNOSTICS: {  // task-specific diagnostics/post-processing
      Stats::startEvent("worker diagnostics");
      doDiagnostics();
      Stats::stopEvent("worker diagnostics");
    } break;
    case GET_L2_NORM: {  // evaluate norm on dfgs and send
      Stats::startEvent("worker get L2 norm");
      sendLpNorms(2);
      Stats::stopEvent("worker get L2 norm");
    } break;
    case GET_L1_NORM: {  // evaluate norm on dfgs and send
      Stats::startEvent("worker get L1 norm");
      sendLpNorms(1);
      Stats::stopEvent("worker get L1 norm");
    } break;
    case GET_MAX_NORM: {  // evaluate norm on dfgs and send
      Stats::startEvent("worker get max norm");
      sendLpNorms(0);
      Stats::stopEvent("worker get max norm");
    } break;
    case PARALLEL_EVAL_NORM: {  // evaluate norms on new dfg and send
      Stats::startEvent("worker parallel eval norm");
      parallelEvalNorm();
      Stats::stopEvent("worker parallel eval norm");
    } break;
    case EVAL_ANALYTICAL_NORM: {  // evaluate analytical norms on new dfg and send
      Stats::startEvent("worker eval analytical norm");
      evalAnalyticalOnDFG();
      Stats::stopEvent("worker eval analytical norm");
    } break;
    case EVAL_ERROR_NORM: {  // evaluate analytical norms on new dfg and send difference
      Stats::startEvent("worker eval error norm");
      evalErrorOnDFG();
      Stats::stopEvent("worker eval error norm");
    } break;
    case INTERPOLATE_VALUES: {  // interpolate values on given coordinates
      Stats::startEvent("worker interpolate values");
      auto values = interpolateValues();
      // send result
      MASTER_EXCLUSIVE_SECTION {
        MPI_Send(values.data(), values.size(), abstraction::getMPIDatatype(
          abstraction::getabstractionDataType<CombiDataType>()), theMPISystem()->getManagerRank(),
          TRANSFER_INTERPOLATION_TAG, theMPISystem()->getGlobalComm());
      }
      Stats::stopEvent("worker interpolate values");
    } break;
    case WRITE_INTERPOLATED_VALUES_PER_GRID: {  // interpolate values on given coordinates and write values to .h5
      Stats::startEvent("worker write interpolated values");
      writeInterpolatedValuesPerGrid();
      Stats::stopEvent("worker write interpolated values");
    } break;
    case RESCHEDULE_ADD_TASK: {
      assert(currentTask_ == nullptr);

      receiveAndInitializeTaskAndFaults(); // receive and initalize new task
			// now the variable currentTask_ contains the newly received task
      currentTask_->setZero();
      fillDFGFromDSGU(currentTask_);
      currentTask_->setFinished(true);
      currentTask_ = nullptr;
    } break;
    case RESCHEDULE_REMOVE_TASK: {
      assert(currentTask_ == nullptr);

      size_t taskID;
      MASTER_EXCLUSIVE_SECTION {
        MPI_Recv(
            &taskID, 1,
            abstraction::getMPIDatatype(abstraction::getabstractionDataType<decltype(taskID)>()),
            theMPISystem()->getManagerRank(), 0, theMPISystem()->getGlobalComm(),
            MPI_STATUS_IGNORE);
      }
      MPI_Bcast(
          &taskID, 1,
          abstraction::getMPIDatatype(abstraction::getabstractionDataType<decltype(taskID)>()),
          theMPISystem()->getMasterRank(), theMPISystem()->getLocalComm());

      // search for task send to group master and remove
      for (size_t i = 0; i < tasks_.size(); ++i) {
        if (tasks_[i]->getID() == taskID) {
          MASTER_EXCLUSIVE_SECTION {
            // send to group master
            Task::send(&tasks_[i], theMPISystem()->getManagerRank(),
                       theMPISystem()->getGlobalComm());
          }
          delete(tasks_[i]);
          tasks_.erase(tasks_.begin() + i);
          break;  // only one task has the taskID
        }
      }
    } break;
    case WRITE_DSG_MINMAX_COEFFICIENTS: {
      std::string filename;
      MASTER_EXCLUSIVE_SECTION {
        MPIUtils::receiveClass(&filename, theMPISystem()->getManagerRank(),
                              theMPISystem()->getGlobalComm());
      }
      for (size_t i = 0; i < combinedUniDSGVector_.size(); ++i) {
        combinedUniDSGVector_[i]->writeMinMaxCoefficents(filename, i);
      }
		} break;
    default: { assert(false && "signal not implemented"); }
  }
  if (isGENE) {
    // special solution for GENE
    // todo: find better solution and remove this
    if ((signal == RUN_FIRST || signal == RUN_NEXT || signal == RECOMPUTE) &&
        !currentTask_->isFinished() && omitReadySignal)
      return signal;
  }
  // in the general case: send ready signal.
  // if(!omitReadySignal)
  ready();
  if (isGENE) {
    if (signal == ADD_TASK) {  // ready resets currentTask but needs to be set for GENE
      currentTask_ = tasks_.back();
    }
  }
  if (!isGENE && signal == RUN_NEXT) {
    Stats::stopEvent("worker run");
  }
  return signal;
}

void ProcessGroupWorker::decideToKill() {
  // decide if processor was killed during this iteration
  currentTask_->decideToKill();
}

void ProcessGroupWorker::ready() {
  if (ENABLE_FT) {
    // with this barrier the local root but also each other process can detect
    // whether a process in the group has failed
    int globalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);
    // std::cout << "rank " << globalRank << " is ready \n";
    int err = simft::Sim_FT_MPI_Barrier(theMPISystem()->getLocalCommFT());

    if (err == MPI_ERR_PROC_FAILED) {
      status_ = PROCESS_GROUP_FAIL;

      std::cout << "rank " << globalRank << " fault detected" << std::endl;
    }
  }
  if (status_ != PROCESS_GROUP_FAIL) {
    // check if there are unfinished tasks
    // all the tasks that are not the first in their process group will be run in this loop
    for (size_t i = 0; i < tasks_.size(); ++i) {
      if (!tasks_[i]->isFinished()) {
        status_ = PROCESS_GROUP_BUSY;

        // set currentTask
        currentTask_ = tasks_[i];
        // if isGENE, this is done in GENE's worker_routines.cpp
        // if (!isGENE) {
        //   Stats::startEvent("worker run");
        // }
        currentTask_->run(theMPISystem()->getLocalComm());
        Stats::Event e;
        // if (!isGENE) {
        //   Stats::stopEvent("worker run");
        // }

        processDuration(*currentTask_, e, theMPISystem()->getNumProcs());
        if (ENABLE_FT) {
          // with this barrier the local root but also each other process can detect
          // whether a process in the group has failed
          int err = simft::Sim_FT_MPI_Barrier(theMPISystem()->getLocalCommFT());

          if (err == MPI_ERR_PROC_FAILED) {
            status_ = PROCESS_GROUP_FAIL;
            break;
          }
        }
        // merge problem?
        // todo: gene specific voodoo
        if (isGENE && !currentTask_->isFinished()) {
          return;
        }
      }
    }

    // all tasks finished -> group waiting
    if (status_ != PROCESS_GROUP_FAIL) {
      status_ = PROCESS_GROUP_WAIT;
    }
  }

  // send ready status to manager
  MASTER_EXCLUSIVE_SECTION {
    StatusType status = status_;

    if (ENABLE_FT) {
      simft::Sim_FT_MPI_Send(&status, 1, MPI_INT, theMPISystem()->getManagerRank(), TRANSFER_STATUS_TAG,
                             theMPISystem()->getGlobalCommFT());
    } else {
      MPI_Send(&status, 1, MPI_INT, theMPISystem()->getManagerRank(), TRANSFER_STATUS_TAG,
               theMPISystem()->getGlobalComm());
    }
  }

  // reset current task
  currentTask_ = NULL;

  // if failed proc in this group detected the alive procs go into recovery state
  if (ENABLE_FT) {
    if (status_ == PROCESS_GROUP_FAIL) {
      theMPISystem()->recoverCommunicators(false);
      status_ = PROCESS_GROUP_WAIT;
    }
  }
}

/**
 * This method reduces the lmax and lmin vectors of the sparse grid according to the reduction
 * specifications in ctparam. It is taken care of that lmin does not fall below 1 and lmax >= lmin.
 * We do not reduce the levels in the last combination as we do not want to lose any information
 * for the final checkpoint.
 */
void reduceSparseGridCoefficients(LevelVector& lmax, LevelVector& lmin,
                                  IndexType totalNumberOfCombis, IndexType currentCombi,
                                  LevelVector reduceLmin, LevelVector reduceLmax) {
  //checking for valid combi step
  assert(currentCombi >= 0);
  if(!(currentCombi < totalNumberOfCombis)) {
    MASTER_EXCLUSIVE_SECTION {
      std::cout << "combining more often than totalNumberOfCombis -- do this for postprocessing only" << std::endl;
    }
  }

  // this if-clause is currently always true, as initCombinedUniDSGVector is called only once,
  // at the beginning of the computation.
  // Leaving it here, in case the SG subspaces need to be re-initialized at some point.
  if (currentCombi < totalNumberOfCombis - 1) {  // do not reduce in last iteration
    for (size_t i = 0; i < reduceLmin.size(); ++i) {
      assert(reduceLmax[i] >= 0 && reduceLmin[i] >= 0);  // check for valid reduce values
      if (lmin[i] > 1) {
        lmin[i] = std::max(static_cast<LevelType>(1), static_cast<LevelType>(lmin[i] - reduceLmin[i]));
      }
    }
    for (size_t i = 0; i < reduceLmax.size(); ++i) {
      lmax[i] = std::max(lmin[i], static_cast<LevelType>(lmax[i] - reduceLmax[i]));
    }
  }
}

void registerAllSubspacesInDSGU(DistributedSparseGridUniform<CombiDataType>& dsgu,
                                const CombiParameters& combiParameters) {
  // the last level vector should have the highest level sum
  const auto highestLevelSum = levelSum(dsgu.getAllLevelVectors().back());
  for (const auto& level : dsgu.getAllLevelVectors()) {
    if (levelSum(level) == highestLevelSum) {
      const auto& boundary = combiParameters.getBoundary();
      auto dfgDecomposition = combigrid::downsampleDecomposition(
          combiParameters.getDecomposition(), combiParameters.getLMax(), level, boundary);
      auto uniDFG = std::unique_ptr<DistributedFullGrid<CombiDataType>>(
          new DistributedFullGrid<CombiDataType>(
              combiParameters.getDim(), level, dsgu.getCommunicator(), boundary,
              combiParameters.getParallelization(), false, dfgDecomposition));
      uniDFG->registerUniformSG(dsgu);
    } else {
      assert(levelSum(level) < highestLevelSum);
    }
  }
}

/** Initializes the dsgu for each species by setting the subspace sizes of all
 * dfgs in the global reduce comm. After calling, all workers which share the
 * same spatial distribution of the dsgu (those who combine during global
 * reduce) have the same sized subspaces and thus share the same sized dsg.
 *
 * Attention: No data is created here, only subspace sizes are shared.
 */
void ProcessGroupWorker::initCombinedUniDSGVector() {
  if (tasks_.size() == 0) {
    std::cout << "Possible error: task size is 0! \n";
  }
  assert(combiParametersSet_);
  // we assume here that every task has the same number of grids, e.g. species in GENE
  LevelVector lmin = combiParameters_.getLMin();
  LevelVector lmax = combiParameters_.getLMax();

  // the dsg can be smaller than lmax because the highest subspaces do not have
  // to be exchanged
  // todo: use a flag to switch on/off optimized combination

  reduceSparseGridCoefficients(lmax, lmin, combiParameters_.getNumberOfCombinations(),
                               currentCombi_, combiParameters_.getLMinReductionVector(),
                               combiParameters_.getLMaxReductionVector());

#ifdef DEBUG_OUTPUT
  MASTER_EXCLUSIVE_SECTION {
    std::cout << "lmin: " << lmin << std::endl;
    std::cout << "lmax: " << lmax << std::endl;
  }
#endif

  // get all subspaces in the (optimized) combischeme, create dsgs
  combinedUniDSGVector_.resize(static_cast<size_t>(combiParameters_.getNumGrids()));
  for (auto& uniDSG : combinedUniDSGVector_) {
    uniDSG = std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>(
        new DistributedSparseGridUniform<CombiDataType>(combiParameters_.getDim(), lmax, lmin,
                                                        theMPISystem()->getLocalComm()));
    // // this registers all possible subspaces in the DSGU
    // // can be used to test the memory consumption of the "filled" DSGU
    // registerAllSubspacesInDSGU(*uniDSG, combiParameters_);
#ifdef DEBUG_OUTPUT
    MASTER_EXCLUSIVE_SECTION {
      std::cout << "dsg size: " << uniDSG->getRawDataSize() << " * " << sizeof(CombiDataType)
                << std::endl;
    }
#endif  // def DEBUG_OUTPUT
  }

  // register dsgs in all dfgs
  Stats::startEvent("worker register dsgus");
  for (size_t g = 0; g < combinedUniDSGVector_.size(); ++g) {
    for (Task* t : tasks_) {
#ifdef DEBUG_OUTPUT
      MASTER_EXCLUSIVE_SECTION { std::cout << "register task " << t->getID() << std::endl; }
#endif  // def DEBUG_OUTPUT
      DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid(static_cast<int>(g));
      // set subspace sizes local
      dfg.registerUniformSG(*(combinedUniDSGVector_[g]));
    }
    // // we may clear the levels_ member of the sparse grids here to save memory
    // // but only if we need no new full grids initialized from the sparse grids!
    // // ...such as for rescheduling or interpolation (parallelEval/ evalNorm / ...)
    // combinedUniDSGVector_[(size_t) g]->resetLevels();
  }
  Stats::stopEvent("worker register dsgus");

  // global reduce of subspace sizes
  CommunicatorType globalReduceComm = theMPISystem()->getGlobalReduceComm();
  for (auto& uniDSG : combinedUniDSGVector_) {
    uniDSG->reduceSubspaceSizes(globalReduceComm);
#ifdef DEBUG_OUTPUT
    MASTER_EXCLUSIVE_SECTION {
      std::cout << "dsg size: " << uniDSG->getRawDataSize() << " * " << sizeof(CombiDataType)
                << std::endl;
    }
#endif  // def DEBUG_OUTPUT
  }
}

void ProcessGroupWorker::hierarchizeFullGrids() {
  auto numGrids = combiParameters_.getNumGrids();
  // real localMax(0.0);
  //  std::vector<CombiDataType> beforeCombi;
  bool anyNotBoundary =
      std::any_of(combiParameters_.getBoundary().begin(), combiParameters_.getBoundary().end(),
                  [](BoundaryType b) { return b == 0; });
  for (Task* t : tasks_) {
    for (IndexType g = 0; g < numGrids; g++) {
      DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid(static_cast<int>(g));

      // hierarchize dfg
      if (anyNotBoundary) {
        DistributedHierarchization::hierarchize<CombiDataType>(
            dfg, combiParameters_.getHierarchizationDims(),
            combiParameters_.getHierarchicalBases());
      } else {
        DistributedHierarchization::hierarchize<CombiDataType>(
            dfg, combiParameters_.getHierarchizationDims(), combiParameters_.getHierarchicalBases(),
            combiParameters_.getLMin());
      }
    }
  }
}

void ProcessGroupWorker::addFullGridsToUniformSG() {
  assert(combinedUniDSGVector_.size() > 0 &&
         "Initialize dsgu first with "
         "initCombinedUniDSGVector()");
  auto numGrids = combiParameters_.getNumGrids();
  for (Task* t : tasks_) {
    for (IndexType g = 0; g < numGrids; g++) {
      DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid(static_cast<int>(g));

      // lokales reduce auf sg ->
      dfg.addToUniformSG(*combinedUniDSGVector_[g], t->getCoefficient());
#ifdef DEBUG_OUTPUT
      std::cout << "Combination: added task " << t->getID() << " with coefficient "
                << t->getCoefficient() << "\n";
#endif
    }
  }
}

void ProcessGroupWorker::reduceUniformSG() {
  // we assume here that every task has the same number of grids, e.g. species in GENE
  auto numGrids = combiParameters_.getNumGrids();

#ifdef DEBUG_OUTPUT
  MASTER_EXCLUSIVE_SECTION { std::cout << "reduce uniform sg \n"; }
#endif

  for (IndexType g = 0; g < numGrids; g++) {
    CombiCom::distributedGlobalReduce(*combinedUniDSGVector_[g]);
    assert(CombiCom::sumAndCheckSubspaceSizes(*combinedUniDSGVector_[g]));
  }
}

void ProcessGroupWorker::combineLocalAndGlobal() {
  assert(combinedUniDSGVector_.size() > 0 && "Initialize dsgu first with "
                                             "initCombinedUniDSGVector()");
  // assert(combinedUniDSGVector_[0]->isSubspaceDataCreated());

#ifdef DEBUG_OUTPUT
  MASTER_EXCLUSIVE_SECTION { std::cout << "start combining \n"; }
#endif
  zeroDsgsData();

  Stats::startEvent("worker hierarchize");
  hierarchizeFullGrids();
  Stats::stopEvent("worker hierarchize");

#ifdef DEBUG_OUTPUT
  MASTER_EXCLUSIVE_SECTION { std::cout << "mid combining \n"; }
#endif

  Stats::startEvent("worker local reduce");
  addFullGridsToUniformSG();
  Stats::stopEvent("worker local reduce");

#ifdef DEBUG_OUTPUT
  MASTER_EXCLUSIVE_SECTION { std::cout << "almost done combining \n"; }
#endif

  Stats::startEvent("worker global reduce");
  reduceUniformSG();
  Stats::stopEvent("worker global reduce");

#ifdef DEBUG_OUTPUT
  MASTER_EXCLUSIVE_SECTION { std::cout << "end combining \n"; }
#endif
}

void ProcessGroupWorker::combineUniform() {
  combineLocalAndGlobal();
  integrateCombinedSolution();
}

void ProcessGroupWorker::parallelEval() {
  if(uniformDecomposition)
    parallelEvalUniform();
  else
    assert(false && "not yet implemented");
}
// cf https://stackoverflow.com/questions/874134/find-out-if-string-ends-with-another-string-in-c
static bool endsWith(const std::string& str, const std::string& suffix)
{
    return str.size() >= suffix.size() && 0 == str.compare(str.size()-suffix.size(), suffix.size(), suffix);
}

// helper function to output bool vector
inline std::ostream& operator<<(std::ostream& os, const std::vector<bool>& l) {
  os << "[";

  for (size_t i = 0; i < l.size(); ++i)
    os << l[i] << " ";

  os << "]";

  return os;
}

LevelVector ProcessGroupWorker::receiveLevalAndBroadcast(){
  const auto dim = combiParameters_.getDim();

  // receive leval and broadcast to group members
  std::vector<int> tmp(dim);
  MASTER_EXCLUSIVE_SECTION {
    MPI_Recv(&tmp[0], static_cast<int>(dim), MPI_INT, theMPISystem()->getManagerRank(), TRANSFER_LEVAL_TAG,
             theMPISystem()->getGlobalComm(), MPI_STATUS_IGNORE);
  }

  MPI_Bcast(&tmp[0], dim, MPI_INT, theMPISystem()->getMasterRank(), theMPISystem()->getLocalComm());
  LevelVector leval(tmp.begin(), tmp.end());
  return leval;
}

void ProcessGroupWorker::fillDFGFromDSGU(DistributedFullGrid<CombiDataType>& dfg, IndexType g) {
  // fill dfg with hierarchical coefficients from distributed sparse grid
  dfg.extractFromUniformSG(*combinedUniDSGVector_[g]);

  bool anyNotBoundary =
      std::any_of(combiParameters_.getBoundary().begin(), combiParameters_.getBoundary().end(),
                  [](BoundaryType b) { return b == 0; });

  if (anyNotBoundary) {
    DistributedHierarchization::dehierarchizeDFG(dfg, combiParameters_.getHierarchizationDims(),
                                                 combiParameters_.getHierarchicalBases());
  } else {
    DistributedHierarchization::dehierarchizeDFG(dfg, combiParameters_.getHierarchizationDims(),
                                                 combiParameters_.getHierarchicalBases(),
                                                 combiParameters_.getLMin());
  }
}

void ProcessGroupWorker::fillDFGFromDSGU(Task* t) {
  auto numGrids = static_cast<int>(
      combiParameters_
          .getNumGrids());  // we assume here that every task has the same number of grids
  for (int g = 0; g < numGrids; g++) {
    assert(combinedUniDSGVector_[g] != nullptr);
    this->fillDFGFromDSGU(t->getDistributedFullGrid(g), g);
  }
}

void ProcessGroupWorker::parallelEvalUniform() {
  assert(uniformDecomposition);

  assert(combiParametersSet_);
  auto numGrids = combiParameters_.getNumGrids();  // we assume here that every task has the same number of grids

  auto leval = receiveLevalAndBroadcast();
  const auto dim = static_cast<DimType>(leval.size());

  // receive filename and broadcast to group members
  std::string filename;
  MASTER_EXCLUSIVE_SECTION {
    MPIUtils::receiveClass(&filename, theMPISystem()->getManagerRank(),
                           theMPISystem()->getGlobalComm());
  }

  MPIUtils::broadcastClass(&filename, theMPISystem()->getMasterRank(),
                           theMPISystem()->getLocalComm());

  for (IndexType g = 0; g < numGrids; g++) {  // loop over all grids and plot them
    // create dfg
    bool forwardDecomposition = combiParameters_.getForwardDecomposition();
    auto levalDecomposition = combigrid::downsampleDecomposition(
            combiParameters_.getDecomposition(),
            combiParameters_.getLMax(), leval,
            combiParameters_.getBoundary());

    DistributedFullGrid<CombiDataType> dfg(
      dim, leval, theMPISystem()->getLocalComm(), combiParameters_.getBoundary(),
      combiParameters_.getParallelization(), forwardDecomposition, levalDecomposition);
    this->fillDFGFromDSGU(dfg, g);
    // save dfg to file with MPI-IO
    if(endsWith(filename, ".vtk")){
      dfg.writePlotFileVTK(filename.c_str());
    }else{
      std::string fn = filename;
      auto pos = fn.find(".");
      if (pos != std::string::npos){
        // if filename contains ".", insert grid number before that
        fn.insert(pos, "_" + std::to_string(g));
      }
      dfg.writePlotFile(fn.c_str());
    }
  }
}

void ProcessGroupWorker::sendLpNorms(int p) {
  // get Lp norm on every worker; reduce through dfg function
  std::vector<double> lpnorms;
  lpnorms.reserve(tasks_.size());
  for (const auto& t : tasks_) {
    auto lpnorm = t->getDistributedFullGrid().getLpNorm(p);
    lpnorms.push_back(lpnorm);
    // std::cout << t->getID() << " ";
  }
  // send from master to manager
  MASTER_EXCLUSIVE_SECTION {
    MPI_Send(lpnorms.data(), static_cast<int>(lpnorms.size()), MPI_DOUBLE,
             theMPISystem()->getManagerRank(), TRANSFER_NORM_TAG, theMPISystem()->getGlobalComm());
  }
}

void sendEvalNorms(const DistributedFullGrid<CombiDataType>& dfg){
  // get Lp norm on every worker; reduce through dfg function
  for (int p = 0; p < 3; ++p) {
    auto lpnorm = dfg.getLpNorm(p);

    // send from master to manager
    MASTER_EXCLUSIVE_SECTION {
      MPI_Send(&lpnorm, 1, MPI_DOUBLE,
              theMPISystem()->getManagerRank(), TRANSFER_NORM_TAG, theMPISystem()->getGlobalComm());
    }
  }
}

void ProcessGroupWorker::parallelEvalNorm() {
  auto leval = receiveLevalAndBroadcast();
  const auto dim = static_cast<DimType>(leval.size());
  bool forwardDecomposition = combiParameters_.getForwardDecomposition();
  auto levalDecomposition = combigrid::downsampleDecomposition(
          combiParameters_.getDecomposition(),
          combiParameters_.getLMax(), leval,
          combiParameters_.getBoundary());

  DistributedFullGrid<CombiDataType> dfg(
      dim, leval, theMPISystem()->getLocalComm(), combiParameters_.getBoundary(),
      combiParameters_.getParallelization(), forwardDecomposition, levalDecomposition);

  this->fillDFGFromDSGU(dfg, 0);

  sendEvalNorms(dfg);
}

void ProcessGroupWorker::evalAnalyticalOnDFG() {
  auto leval = receiveLevalAndBroadcast();
  const auto dim = static_cast<DimType>(leval.size());
  bool forwardDecomposition = combiParameters_.getForwardDecomposition();
  auto levalDecomposition = combigrid::downsampleDecomposition(
          combiParameters_.getDecomposition(),
          combiParameters_.getLMax(), leval,
          combiParameters_.getBoundary());

  DistributedFullGrid<CombiDataType> dfg(
      dim, leval, theMPISystem()->getLocalComm(), combiParameters_.getBoundary(),
      combiParameters_.getParallelization(), forwardDecomposition, levalDecomposition);

  // interpolate Task's analyticalSolution
  for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
    std::vector<double> coords(leval.size());
    dfg.getCoordsLocal(li, coords);

    dfg.getData()[li] = tasks_[0]->analyticalSolution(coords, 0);
  }

  sendEvalNorms(dfg);
}

void ProcessGroupWorker::evalErrorOnDFG() {
  auto leval = receiveLevalAndBroadcast();
  const auto dim = static_cast<DimType>(leval.size());
  bool forwardDecomposition = combiParameters_.getForwardDecomposition();
  auto levalDecomposition = combigrid::downsampleDecomposition(
          combiParameters_.getDecomposition(),
          combiParameters_.getLMax(), leval,
          combiParameters_.getBoundary());

  DistributedFullGrid<CombiDataType> dfg(
      dim, leval, theMPISystem()->getLocalComm(), combiParameters_.getBoundary(),
      combiParameters_.getParallelization(), forwardDecomposition, levalDecomposition);

  this->fillDFGFromDSGU(dfg, 0);
  // interpolate Task's analyticalSolution
  for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
    std::vector<double> coords(leval.size());
    dfg.getCoordsLocal(li, coords);

    dfg.getData()[li] -= tasks_[0]->analyticalSolution(coords, 0);
  }

  sendEvalNorms(dfg);
}

void ProcessGroupWorker::doDiagnostics() {
  // receive taskID and broadcast
  size_t taskID;
  MASTER_EXCLUSIVE_SECTION {
    MPI_Recv(&taskID, 1,
             abstraction::getMPIDatatype(abstraction::getabstractionDataType<decltype(taskID)>()),
             theMPISystem()->getManagerRank(), 0, theMPISystem()->getGlobalComm(),
             MPI_STATUS_IGNORE);
  }
  MPI_Bcast(&taskID, 1,
            abstraction::getMPIDatatype(abstraction::getabstractionDataType<decltype(taskID)>()),
            theMPISystem()->getMasterRank(), theMPISystem()->getLocalComm());

  // call diagnostics on that Task
  for (auto task : tasks_) {
    if (task->getID() == taskID) {
      std::vector<DistributedSparseGridUniform<CombiDataType>*> dsgsToPassToTask;
      for (auto& dsgPtr : combinedUniDSGVector_) {
        dsgsToPassToTask.push_back(dsgPtr.get());
      }
      task->doDiagnostics(dsgsToPassToTask, combiParameters_.getHierarchizationDims());
      return;
    }
  }
  assert(false && "this taskID is not here");
}

void receiveAndBroadcastInterpolationCoords(std::vector<std::vector<real>>& interpolationCoords, DimType dim) {
  std::vector<real> interpolationCoordsSerial;
  auto realType = abstraction::getMPIDatatype(abstraction::getabstractionDataType<real>());
  int coordsSize;
  MASTER_EXCLUSIVE_SECTION {
    MPI_Status status;
    MPI_Probe(theMPISystem()->getManagerRank(), TRANSFER_INTERPOLATION_TAG, theMPISystem()->getGlobalComm(), &status);
    MPI_Get_count(&status, realType, &coordsSize);

    // resize buffer to appropriate size and receive
    interpolationCoordsSerial.resize(coordsSize);
    int result = MPI_Recv(interpolationCoordsSerial.data(), coordsSize, realType, theMPISystem()->getManagerRank(),
             TRANSFER_INTERPOLATION_TAG, theMPISystem()->getGlobalComm(), &status);
    assert(result == MPI_SUCCESS);
    assert(status.MPI_ERROR == MPI_SUCCESS);
    for (const auto& coord : interpolationCoordsSerial) {
      assert(coord >= 0.0 && coord <= 1.0);
    }
  }
  // broadcast size of vector, and then vector
  MPI_Bcast(&coordsSize, 1, MPI_INT, theMPISystem()->getMasterRank(),
            theMPISystem()->getLocalComm());
  interpolationCoordsSerial.resize(coordsSize);
  MPI_Bcast(interpolationCoordsSerial.data(), coordsSize, realType, theMPISystem()->getMasterRank(),
            theMPISystem()->getLocalComm());
  for (const auto& coord : interpolationCoordsSerial) {
    assert(coord >= 0.0 && coord <= 1.0);
  }

  // split vector into coordinates
  const int dimInt = static_cast<int>(dim);
  auto numCoordinates = coordsSize / dimInt;
  assert( coordsSize % dimInt == 0);
  interpolationCoords.resize(numCoordinates);
  auto it = interpolationCoordsSerial.cbegin();
  for (auto& coord: interpolationCoords) {
    coord.insert(coord.end(), it, it+dimInt);
    it += dimInt;
  }
  assert(it == interpolationCoordsSerial.end());
}

std::vector<CombiDataType> ProcessGroupWorker::interpolateValues() {
  assert(combiParameters_.getNumGrids() == 1 && "interpolate only implemented for 1 species!");
  // receive coordinates and broadcast to group members
  std::vector<std::vector<real>> interpolationCoords;
  receiveAndBroadcastInterpolationCoords(interpolationCoords, combiParameters_.getDim());
  auto numCoordinates = interpolationCoords.size();

  // call interpolation function on tasks and reduce with combination coefficient
  std::vector<CombiDataType> values(numCoordinates, 0.);
  for (Task* t : tasks_){
    auto coeff = t->getCoefficient();
    for (size_t i = 0; i < numCoordinates; ++i) {
      values[i] += t->getDistributedFullGrid().evalLocal(interpolationCoords[i]) * coeff;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, values.data(), static_cast<int>(numCoordinates),
                abstraction::getMPIDatatype(abstraction::getabstractionDataType<CombiDataType>()),
                MPI_SUM, theMPISystem()->getLocalComm());
  return values;
}

void ProcessGroupWorker::writeInterpolatedValuesPerGrid() {
#ifdef HAVE_HIGHFIVE
  assert(combiParameters_.getNumGrids() == 1 && "interpolate only implemented for 1 species!");
  // receive coordinates and broadcast to group members
  std::vector<std::vector<real>> interpolationCoords;
  receiveAndBroadcastInterpolationCoords(interpolationCoords, combiParameters_.getDim());

  // call interpolation function on tasks and write out task-wise
  for (size_t i = 0; i < tasks_.size(); ++i) {
    auto taskVals = tasks_[i]->getDistributedFullGrid().getInterpolatedValues(interpolationCoords);
    if (i % (theMPISystem()->getNumProcs()) == theMPISystem()->getLocalRank()) {
      // generate a rank-local per-run random number
      // std::random_device dev;
      static std::mt19937 rng(std::chrono::high_resolution_clock::now().time_since_epoch().count());
      static std::uniform_int_distribution<std::mt19937::result_type> dist(
          1, std::numeric_limits<size_t>::max());
      static size_t rankLocalRandom = dist(rng);

      std::string saveFilePath = "interpolated_" + std::to_string(tasks_[i]->getID()) + ".h5";
      // check if file already exists, if no, create
      HighFive::File h5_file(saveFilePath,
                             HighFive::File::OpenOrCreate | HighFive::File::ReadWrite);

      // std::vector<std::string> elems = h5_file.listObjectNames();
      // std::cout << elems << std::endl;

      std::string groupName = "run_" + std::to_string(rankLocalRandom);
      HighFive::Group group;
      if (h5_file.exist(groupName)) {
        group = h5_file.getGroup(groupName);
      } else {
        group = h5_file.createGroup(groupName);
      }

      std::string datasetName = "interpolated_" + std::to_string(currentCombi_);
      HighFive::DataSet dataset =
          group.createDataSet<CombiDataType>(datasetName, HighFive::DataSpace::From(taskVals));

      dataset.write(taskVals);

      HighFive::Attribute a2 = dataset.createAttribute<combigrid::real>(
          "simulation_time", HighFive::DataSpace::From(tasks_[i]->getCurrentTime()));
      a2.write(tasks_[i]->getCurrentTime());
    }
  }
#else  // if not compiled with hdf5
  throw std::runtime_error("requesting hdf5 write but built without hdf5 support");
#endif
}

void ProcessGroupWorker::gridEval() {  // not supported anymore
  /* error if no tasks available
   * todo: however, this is not a real problem, we could can create an empty
   * grid an contribute to the reduce operation. at the moment even the dim
   * parameter is stored in the tasks, so if no task available we have no access
   * to this parameter.
   */
  assert(tasks_.size() > 0);

  assert(combiParametersSet_);
  const DimType dim = combiParameters_.getDim();

  LevelVector leval(dim);

  // receive leval
  MASTER_EXCLUSIVE_SECTION {
    // receive size of levelvector = dimensionality
    MPI_Status status;
    int bsize;
    MPI_Probe(theMPISystem()->getManagerRank(), TRANSFER_LEVAL_TAG, theMPISystem()->getGlobalComm(), &status);
    MPI_Get_count(&status, MPI_INT, &bsize);

    assert(bsize == static_cast<int>(dim));

    std::vector<int> tmp(dim);
    MPI_Recv(&tmp[0], bsize, MPI_INT, theMPISystem()->getManagerRank(), TRANSFER_LEVAL_TAG,
             theMPISystem()->getGlobalComm(), MPI_STATUS_IGNORE);
    leval = LevelVector(tmp.begin(), tmp.end());
  }

  assert(combiParametersSet_);
  const std::vector<BoundaryType>& boundary = combiParameters_.getBoundary();
  FullGrid<CombiDataType> fg_red(dim, leval, boundary);

  // create the empty grid on only on localroot
  MASTER_EXCLUSIVE_SECTION { fg_red.createFullGrid(); }

  // collect fg on pgrouproot and reduce
  for (size_t i = 0; i < tasks_.size(); ++i) {
    Task* t = tasks_[i];

    FullGrid<CombiDataType> fg(t->getDim(), t->getLevelVector(), boundary);

    MASTER_EXCLUSIVE_SECTION { fg.createFullGrid(); }

    t->getFullGrid(fg, theMPISystem()->getMasterRank(), theMPISystem()->getLocalComm());

    MASTER_EXCLUSIVE_SECTION { fg_red.add(fg, t->getCoefficient()); }
  }
  // global reduce of f_red
  MASTER_EXCLUSIVE_SECTION {
    CombiCom::FGReduce(fg_red, theMPISystem()->getManagerRank(), theMPISystem()->getGlobalComm());
  }
}

void ProcessGroupWorker::receiveAndInitializeTaskAndFaults(bool mayAlreadyExist /*=true*/) {
  Task* t;

  // local root receives task
  MASTER_EXCLUSIVE_SECTION {
    Task::receive(&t, theMPISystem()->getManagerRank(), theMPISystem()->getGlobalComm());
  }

  // broadcast task to other process of pgroup
  Task::broadcast(&t, theMPISystem()->getMasterRank(), theMPISystem()->getLocalComm());

  if (!mayAlreadyExist) {
    // check if task already exists on this group
    for (auto tmp : tasks_) assert(tmp->getID() != t->getID());
  }

  MPI_Barrier(theMPISystem()->getLocalComm());
  initializeTaskAndFaults(t);
}

void ProcessGroupWorker::initializeTaskAndFaults(Task* t) {
  assert(combiParametersSet_);
  // add task to task storage
  tasks_.push_back(t);

  // set currentTask
  currentTask_ = tasks_.back();

  // initalize task
  auto taskDecomposition = combigrid::downsampleDecomposition(
          combiParameters_.getDecomposition(),
          combiParameters_.getLMax(), currentTask_->getLevelVector(),
          combiParameters_.getBoundary());
  currentTask_->init(theMPISystem()->getLocalComm(), taskDecomposition);
  if (ENABLE_FT) {
    t_fault_ = currentTask_->initFaults(t_fault_, startTimeIteration_);
  }
}

// todo: this is just a temporary function which will drop out some day
// also this function requires a modified fgreduce method which uses allreduce
// instead reduce in manger
void ProcessGroupWorker::combineFG() {
  // gridEval();

  // TODO: Sync back to fullgrids
}

void ProcessGroupWorker::deleteTasks() {
      // freeing tasks
      for (auto tmp : tasks_) delete (tmp);
      tasks_.clear();
}

void ProcessGroupWorker::updateCombiParameters() {
  // local root receives combi parameters
  MASTER_EXCLUSIVE_SECTION {
    MPIUtils::receiveClass(&combiParameters_, theMPISystem()->getManagerRank(), theMPISystem()->getGlobalComm());
    //std::cout << "master received combiparameters \n";
  }

  // broadcast parameters to other process of pgroup
  MPIUtils::broadcastClass(&combiParameters_, theMPISystem()->getMasterRank(), theMPISystem()->getLocalComm());
  //std::cout << "worker received combiparameters \n";

  combiParametersSet_ = true;

  // overwrite local comm with cartesian communicator
  if (!isGENE  && combiParameters_.isParallelizationSet()){
    // cf. https://www.rookiehpc.com/mpi/docs/mpi_cart_create.php
    // get decompositon from combi params
    auto par = combiParameters_.getParallelization();

    // important: note reverse ordering of dims! -- cf DistributedFullGrid //TODO(pollinta) remove reverse ordering
    std::vector<int> dims (par.size());
    if (reverseOrderingDFGPartitions) {
      dims.assign(par.rbegin(), par.rend());
    } else {
      dims.assign(par.begin(), par.end());
    }

    std::vector<int> periods;
    // Make dimensions periodic depending on boundary parameters
    for (const auto& b : combiParameters_.getBoundary()) {
      if (b == 1) {
        periods.push_back(1);

      } else {
        periods.push_back(0);
      }
    }
    assert(combiParameters_.getDim() == static_cast<DimType>(periods.size()));
    assert(periods.size() == dims.size());

    // don't let MPI assign arbitrary ranks
    int reorder = false;

    // Create a communicator given the topology.
    MPI_Comm new_communicator;
    MPI_Cart_create(theMPISystem()->getLocalComm(), combiParameters_.getDim(), dims.data(),
                    periods.data(), reorder, &new_communicator);

    theMPISystem()->storeLocalComm(new_communicator);
  }
}

void ProcessGroupWorker::integrateCombinedSolution() {
  auto numGrids = static_cast<int>(combiParameters_.getNumGrids());
  for (Task* taskToUpdate : tasks_) {
    for (int g = 0; g < numGrids; g++) {
      // fill dfg with hierarchical coefficients from distributed sparse grid
      taskToUpdate->getDistributedFullGrid(g).extractFromUniformSG(*combinedUniDSGVector_[g]);
    }
  }

  bool anyNotBoundary =
      std::any_of(combiParameters_.getBoundary().begin(), combiParameters_.getBoundary().end(),
                  [](BoundaryType b) { return b == 0; });
  Stats::startEvent("worker dehierarchize");
  for (Task* taskToUpdate : tasks_) {
    for (int g = 0; g < numGrids; g++) {
      if (anyNotBoundary) {
        DistributedHierarchization::dehierarchizeDFG(taskToUpdate->getDistributedFullGrid(g),
                                                     combiParameters_.getHierarchizationDims(),
                                                     combiParameters_.getHierarchicalBases());
      } else {
        DistributedHierarchization::dehierarchizeDFG(
            taskToUpdate->getDistributedFullGrid(g), combiParameters_.getHierarchizationDims(),
            combiParameters_.getHierarchicalBases(), combiParameters_.getLMin());
      }
    }
  }
  Stats::stopEvent("worker dehierarchize");
}

template <typename T>
void copyFromCombinedDSGUToExtraDSGU(const DistributedSparseGridUniform<T>& combinedUniDSG,
                                     DistributedSparseGridUniform<T>& extraUniDSG) {
  // copy partial data from uniDSG to extraDSG
  for (decltype(combinedUniDSG.getNumSubspaces()) i = 0; i < combinedUniDSG.getNumSubspaces(); ++i) {
    assert(extraUniDSG.getDataSize(i) == 0 ||
           extraUniDSG.getDataSize(i) == combinedUniDSG.getDataSize(i));
    std::copy_n(combinedUniDSG.getData(i), extraUniDSG.getDataSize(i), extraUniDSG.getData(i));
  }
}

template <typename T>
void copyFromExtraDSGUToCombinedDSGU(const DistributedSparseGridUniform<T>& extraUniDSG,
                                       DistributedSparseGridUniform<T>& combinedUniDSG) {
  // copy partial data from extraDSG back to uniDSG
  for (decltype(combinedUniDSG.getNumSubspaces()) i = 0; i < combinedUniDSG.getNumSubspaces(); ++i) {
    std::copy_n(extraUniDSG.getData(i), extraUniDSG.getDataSize(i), combinedUniDSG.getData(i));
  }
}

void ProcessGroupWorker::combineThirdLevel() {
  assert(combinedUniDSGVector_.size() != 0);
  assert(combiParametersSet_);

  assert(theMPISystem()->getThirdLevelComms().size() == 1 && "init thirdLevel communicator failed");
  const CommunicatorType& managerComm = theMPISystem()->getThirdLevelComms()[0];
  const CommunicatorType& globalReduceComm = theMPISystem()->getGlobalReduceComm();
  const RankType& globalReduceRank = theMPISystem()->getGlobalReduceRank();
  const RankType& manager = theMPISystem()->getThirdLevelManagerRank();

  std::vector<MPI_Request> requests;
  for (size_t i = 0; i < combinedUniDSGVector_.size(); ++i) {
    auto uniDsg = combinedUniDSGVector_[i].get();
    auto dsgToUse = uniDsg;
    if (extraUniDSGVector_.size() > 0) {
      dsgToUse = extraUniDSGVector_[i].get();
    }
    assert(dsgToUse->getRawDataSize() < 2147483647 && "Dsg is larger than 2^31-1 and can not be "
                                             "sent in a single MPI Call (not "
                                             "supported yet) try a more coarse"
                                             "decomposition");
    // if we have an extra dsg for third level exchange, we use it
    if (extraUniDSGVector_.size() > 0) {
      copyFromCombinedDSGUToExtraDSGU(*uniDsg, *dsgToUse);
    }

    // send dsg data to manager
    Stats::startEvent("worker send dsg data");
    sendDsgData(dsgToUse, manager, managerComm);
    Stats::stopEvent("worker send dsg data");

    // recv combined dsgu from manager
    Stats::startEvent("worker recv dsg data");
    recvDsgData(dsgToUse, manager, managerComm);
    Stats::stopEvent("worker recv dsg data");

    if (extraUniDSGVector_.size() > 0) {
      copyFromExtraDSGUToCombinedDSGU(*dsgToUse, *uniDsg);
    }

    // distribute solution in globalReduceComm to other pgs
    auto request = asyncBcastDsgData(uniDsg, globalReduceRank, globalReduceComm);
    requests.push_back(request);
  }
  // update fgs
  integrateCombinedSolution();

  // wait for bcasts to other pgs in globalReduceComm
  Stats::startEvent("worker wait for bcasts");
  for (MPI_Request& request : requests) {
    auto returnedValue = MPI_Wait(&request, MPI_STATUS_IGNORE);
    assert(returnedValue == MPI_SUCCESS);
  }
  Stats::stopEvent("worker wait for bcasts");
}

void ProcessGroupWorker::combineThirdLevelFileBased() {
  // TODO(pollinta) steal parts from above functions
  // and write, wait, read, and combine
  throw std::runtime_error("combineThirdLevelFileBased not implemented yet");
}

/** Reduces subspace sizes with remote.
 */
void ProcessGroupWorker::reduceSubspaceSizesThirdLevel(bool thirdLevelExtraSparseGrid) {
  assert(combiParametersSet_);
  // update either old or new sparse grids
  auto* uniDSGVectorToSet = &combinedUniDSGVector_;
  if (thirdLevelExtraSparseGrid) {
    if (extraUniDSGVector_.empty()) {
      // create new vector for extra sparse grids (that will be only on this process group)
      extraUniDSGVector_.resize(combinedUniDSGVector_.size());
      for (auto& extraUniDSG : extraUniDSGVector_) {
        extraUniDSG = std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>(
            new DistributedSparseGridUniform<CombiDataType>(
                combinedUniDSGVector_[0]->getDim(), combinedUniDSGVector_[0]->getAllLevelVectors(),
                theMPISystem()->getLocalComm()));
      }
    }
    uniDSGVectorToSet = &extraUniDSGVector_;
  }

  if (uniDSGVectorToSet->size() != 1) {
    throw std::runtime_error(
        "uniDSGVectorToSet.size() > 1 -- not implemented on pg manager's side");
  }

  // prepare for MPI calls to manager
  CommunicatorType thirdLevelComm = theMPISystem()->getThirdLevelComms()[0];
  RankType thirdLevelManagerRank = theMPISystem()->getThirdLevelManagerRank();
  for (size_t i = 0; i < uniDSGVectorToSet->size(); ++i) {
    combinedUniDSGVector_[i]->sendDsgSizesWithGather(thirdLevelComm, thirdLevelManagerRank);
    // set updated sizes in dsgs
    (*uniDSGVectorToSet)[i]->receiveDsgSizesWithScatter(thirdLevelComm, thirdLevelManagerRank);
  }

  if (!thirdLevelExtraSparseGrid) {
    // distribute updated sizes to workers with same decomposition (global reduce comm)
    // cf. waitForThirdLevelSizeUpdate(), which is called in other process groups
    CommunicatorType globalReduceComm = theMPISystem()->getGlobalReduceComm();
    RankType globalReduceRank = theMPISystem()->getGlobalReduceRank();
    for (auto& dsg : combinedUniDSGVector_) {
      dsg->broadcastDsgSizes(globalReduceComm, globalReduceRank);
    }
  }
}

void ProcessGroupWorker::waitForThirdLevelSizeUpdate() {
  RankType thirdLevelPG = (RankType)combiParameters_.getThirdLevelPG();
  CommunicatorType globalReduceComm = theMPISystem()->getGlobalReduceComm();

  for (auto& dsg : combinedUniDSGVector_) {
    dsg->broadcastDsgSizes(globalReduceComm, thirdLevelPG);
  }
}

void ProcessGroupWorker::waitForThirdLevelCombiResult() {
  // receive third level combi result from third level pgroup (global reduce comm)
  RankType thirdLevelPG = (RankType)combiParameters_.getThirdLevelPG();
  CommunicatorType globalReduceComm = theMPISystem()->getGlobalReduceComm();

  for (auto& dsg : combinedUniDSGVector_) {
    auto request = asyncBcastDsgData(dsg.get(), thirdLevelPG, globalReduceComm);
    auto returnedValue = MPI_Wait(&request, MPI_STATUS_IGNORE);
    assert(returnedValue == MPI_SUCCESS);
  }

  integrateCombinedSolution();
}

void ProcessGroupWorker::zeroDsgsData() {
  for (auto& dsg : combinedUniDSGVector_)
    dsg->setZero();
}

/** free dsgus space */
void ProcessGroupWorker::deleteDsgsData() {
  for (auto& dsg : combinedUniDSGVector_)
    dsg->deleteSubspaceData();
}

void ProcessGroupWorker::writeVTKPlotFileOfTask(Task& task) {
#ifdef USE_VTK
  IndexType numGrids = combiParameters_.getNumGrids();
  for (IndexType g = 0; g < numGrids; g++) {
    DistributedFullGrid<CombiDataType>& dfg = task.getDistributedFullGrid(static_cast<int>(g));
    DFGPlotFileWriter<CombiDataType> writer {dfg, g};
    writer.writePlotFile();
  }
#else
  std::cout << "Warning: no vtk output produced as DisCoTec was compiled without VTK." << std::endl;
#endif /* USE_VTK */
}

void ProcessGroupWorker::writeVTKPlotFilesOfAllTasks() {
  for (Task* task : tasks_) writeVTKPlotFileOfTask(*task);
}

void ProcessGroupWorker::writeDSGsToDisk(std::string filenamePrefix) {
  for (size_t i = 0; i < combinedUniDSGVector_.size(); ++i) {
    auto uniDsg = combinedUniDSGVector_[i].get();
    auto dsgToUse = uniDsg;
    if (extraUniDSGVector_.size() > 0) {
      dsgToUse = extraUniDSGVector_[i].get();
      copyFromCombinedDSGUToExtraDSGU(*uniDsg, *dsgToUse);
    }
    dsgToUse->writeOneFileToDisk(filenamePrefix + "_" + std::to_string(i));
  }
}

void ProcessGroupWorker::readDSGsFromDisk(std::string filenamePrefix) {
  for (size_t i = 0; i < combinedUniDSGVector_.size(); ++i) {
    auto uniDsg = combinedUniDSGVector_[i].get();
    auto dsgToUse = uniDsg;
    if (extraUniDSGVector_.size() > 0) {
      dsgToUse = extraUniDSGVector_[i].get();
    }
    dsgToUse->readOneFileFromDisk(filenamePrefix + "_" + std::to_string(i));
    if (extraUniDSGVector_.size() > 0) {
      copyFromExtraDSGUToCombinedDSGU(*dsgToUse, *uniDsg);
    }
  }
}

void ProcessGroupWorker::readDSGsFromDiskAndReduce(std::string filenamePrefixToRead) {
  for (size_t i = 0; i < combinedUniDSGVector_.size(); ++i) {
    auto uniDsg = combinedUniDSGVector_[i].get();
    auto dsgToUse = uniDsg;
    if (extraUniDSGVector_.size() > 0) {
      dsgToUse = extraUniDSGVector_[i].get();
    }
    dsgToUse->readOneFileFromDiskAndReduce(filenamePrefixToRead + "_" + std::to_string(i));
    if (extraUniDSGVector_.size() > 0) {
      copyFromExtraDSGUToCombinedDSGU(*dsgToUse, *uniDsg);
    }
  }
}

} /* namespace combigrid */
