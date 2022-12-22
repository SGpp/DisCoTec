#include "manager/ProcessGroupWorker.hpp"

#include "boost/lexical_cast.hpp"

#include "combicom/CombiCom.hpp"
#include "fullgrid/FullGrid.hpp"
#include "hierarchization/DistributedHierarchization.hpp"
#include "manager/CombiParameters.hpp"
#include "manager/ProcessGroupSignals.hpp"
#include "mpi/MPIUtils.hpp"
#include "sparsegrid/DistributedSparseGridUniform.hpp"
#include "loadmodel/LearningLoadModel.hpp"
#include "mpi/MPISystem.hpp"


#include <algorithm>
#include <iostream>
#include <string>
#ifdef HAVE_HIGHFIVE
#include <chrono>
#include <random>
// highfive is a C++ hdf5 wrapper, available in spack (-> configure with right boost and mpi versions)
#include <highfive/H5File.hpp>
#endif
#include "mpi_fault_simulator/MPI-FT.h"

namespace combigrid {

std::string receiveStringFromManagerAndBroadcastToGroup() {
  std::string stringToReceive;
  MASTER_EXCLUSIVE_SECTION {
    MPIUtils::receiveClass(&stringToReceive, theMPISystem()->getManagerRank(),
                           theMPISystem()->getGlobalComm());
  }
  MPIUtils::broadcastClass(&stringToReceive, theMPISystem()->getMasterRank(),
                           theMPISystem()->getLocalComm());
  return stringToReceive;
}

ProcessGroupWorker::ProcessGroupWorker()
    : currentTask_(nullptr),
      status_(PROCESS_GROUP_WAIT),
      combinedUniDSGVector_(0),
      combinedFGexists_(false),
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

        if (!isGENE) {
          Stats::stopEvent("worker run");
        }
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
          tasksStream <<t->getID() << ": " << combiParameters_.getCoeff(t->getID()) << t->getLevelVector()  << "; ";
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
      Stats::startEvent("initializeCombinedUniDSGVector");
      initCombinedUniDSGVector();
      Stats::stopEvent("initializeCombinedUniDSGVector");

    } break;
    case COMBINE: {  // start combination
      Stats::startEvent("worker combine");
      combineUniform();
      currentCombi_++;
      Stats::stopEvent("worker combine");

    } break;
    case WRITE_DSGS_TO_DISK: {
      Stats::startEvent("worker write to disk");
      std::string filenamePrefix = receiveStringFromManagerAndBroadcastToGroup();
      writeDSGsToDisk(filenamePrefix);
      Stats::stopEvent("worker write to disk");
    } break;
    case READ_DSGS_FROM_DISK: {
      Stats::startEvent("worker read from disk");
      std::string filenamePrefix = receiveStringFromManagerAndBroadcastToGroup();
      readDSGsFromDisk(filenamePrefix);
      Stats::stopEvent("worker read from disk");
    } break;
    case GRID_EVAL: {  // not supported anymore

      Stats::startEvent("eval");
      gridEval();
      Stats::stopEvent("eval");

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
      Stats::startEvent("parallel eval");
      parallelEval();
      Stats::stopEvent("parallel eval");
    } break;
    case DO_DIAGNOSTICS: {  // task-specific diagnostics/post-processing
      Stats::startEvent("diagnostics");
      doDiagnostics();
      Stats::stopEvent("diagnostics");
    } break;
    case GET_L2_NORM: {  // evaluate norm on dfgs and send
      Stats::startEvent("get L2 norm");
      sendLpNorms(2);
      Stats::stopEvent("get L2 norm");
    } break;
    case GET_L1_NORM: {  // evaluate norm on dfgs and send
      Stats::startEvent("get L1 norm");
      sendLpNorms(1);
      Stats::stopEvent("get L1 norm");
    } break;
    case GET_MAX_NORM: {  // evaluate norm on dfgs and send
      Stats::startEvent("get max norm");
      sendLpNorms(0);
      Stats::stopEvent("get max norm");
    } break;
    case PARALLEL_EVAL_NORM: {  // evaluate norms on new dfg and send
      Stats::startEvent("parallel eval norm");
      parallelEvalNorm();
      Stats::stopEvent("parallel eval norm");
    } break;
    case EVAL_ANALYTICAL_NORM: {  // evaluate analytical norms on new dfg and send
      Stats::startEvent("eval analytical norm");
      evalAnalyticalOnDFG();
      Stats::stopEvent("eval analytical norm");
    } break;
    case EVAL_ERROR_NORM: {  // evaluate analytical norms on new dfg and send difference
      Stats::startEvent("eval error norm");
      evalErrorOnDFG();
      Stats::stopEvent("eval error norm");
    } break;
    case INTERPOLATE_VALUES: {  // interpolate values on given coordinates
      Stats::startEvent("worker interpolate values");
      auto values = interpolateValues();
      Stats::stopEvent("worker interpolate values");
    } break;
    case INTERPOLATE_VALUES_AND_SEND_BACK: {
      Stats::startEvent("worker interpolate values");
      auto values = interpolateValues();
      // send result
      MASTER_EXCLUSIVE_SECTION {
        MPI_Send(values.data(), values.size(),
                 abstraction::getMPIDatatype(abstraction::getabstractionDataType<CombiDataType>()),
                 theMPISystem()->getManagerRank(), TRANSFER_INTERPOLATION_TAG,
                 theMPISystem()->getGlobalComm());
      }
      Stats::stopEvent("worker interpolate values");
    } break;
    case INTERPOLATE_VALUES_AND_WRITE_SINGLE_FILE: {
      Stats::startEvent("worker interpolate values");
      auto values = interpolateValues();
      // write result
      MASTER_EXCLUSIVE_SECTION {
        std::string valuesWriteFilename =
            "interpolated_values_" + std::to_string(currentCombi_) + ".h5";
        writeInterpolatedValues(values, valuesWriteFilename);
      }
      Stats::stopEvent("worker interpolate values");
    } break;
    case WRITE_INTERPOLATED_VALUES_PER_GRID: {  // interpolate values on given coordinates and write
                                                // values to .h5
      Stats::startEvent("worker write interpolated values");
      writeInterpolatedValuesPerGrid();
      Stats::stopEvent("write interpolated values");
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
        if (!isGENE) {
          Stats::startEvent("worker run");
        }
        currentTask_->run(theMPISystem()->getLocalComm());
        Stats::Event e;
        if (!isGENE) {
          Stats::stopEvent("worker run");
        }

        // std::cout << "from ready ";
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
        //
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
      dsgu.registerDistributedFullGrid(*uniDFG);
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

  Stats::startEvent("create dsgus");
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
  Stats::stopEvent("create dsgus");

  // register dsgs in all dfgs
  Stats::startEvent("register dsgus");
  for (size_t g = 0; g < combinedUniDSGVector_.size(); ++g) {
    for (Task* t : tasks_) {
#ifdef DEBUG_OUTPUT
      MASTER_EXCLUSIVE_SECTION { std::cout << "register task " << t->getID() << std::endl; }
#endif  // def DEBUG_OUTPUT
      DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid(static_cast<int>(g));
      // set subspace sizes locally
      combinedUniDSGVector_[g]->registerDistributedFullGrid(dfg);
    }
  }
  Stats::stopEvent("register dsgus");

  // global reduce of subspace sizes
  Stats::startEvent("reduce dsgus");
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
  Stats::stopEvent("reduce dsgus");
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
      combinedUniDSGVector_[g]->addDistributedFullGrid(dfg, combiParameters_.getCoeff(t->getID()));
#ifdef DEBUG_OUTPUT
      std::cout << "Combination: added task " << t->getID() << " with coefficient "
                << combiParameters_.getCoeff(t->getID()) << "\n";
#endif
    }
  }
}

void ProcessGroupWorker::reduceUniformSG() {
  // we assume here that every task has the same number of grids, e.g. species in GENE
  auto numGrids = combiParameters_.getNumGrids();

  for (IndexType g = 0; g < numGrids; g++) {
    CombiCom::distributedGlobalReduce(*combinedUniDSGVector_[g]);
    assert(CombiCom::sumAndCheckSubspaceSizes(*combinedUniDSGVector_[g]));
  }
}


void ProcessGroupWorker::combineLocalAndGlobal() {
#ifdef DEBUG_OUTPUT
  MASTER_EXCLUSIVE_SECTION { std::cout << "start combining \n"; }
#endif
  Stats::startEvent("combine zeroDsgsData");
  zeroDsgsData();
  Stats::stopEvent("combine zeroDsgsData");

  Stats::startEvent("combine hierarchize");
  hierarchizeFullGrids();
  Stats::stopEvent("combine hierarchize");

#ifdef DEBUG_OUTPUT
  MASTER_EXCLUSIVE_SECTION { std::cout << "mid combining \n"; }
#endif

  Stats::startEvent("combine local reduce");
  addFullGridsToUniformSG();
  Stats::stopEvent("combine local reduce");

#ifdef DEBUG_OUTPUT
  MASTER_EXCLUSIVE_SECTION { std::cout << "almost done combining \n"; }
#endif

  Stats::startEvent("combine global reduce");
  reduceUniformSG();
  Stats::stopEvent("combine global reduce");

#ifdef DEBUG_OUTPUT
  MASTER_EXCLUSIVE_SECTION { std::cout << "end combining \n"; }
#endif
}

void ProcessGroupWorker::combineUniform() {
  combineLocalAndGlobal();
  integrateCombinedSolution();
}


void ProcessGroupWorker::parallelEval() {
  if (uniformDecomposition)
    parallelEvalUniform();
  else
    assert(false && "not yet implemented");
}

// cf https://stackoverflow.com/questions/874134/find-out-if-string-ends-with-another-string-in-c
static bool endsWith(const std::string& str, const std::string& suffix)
{
    return str.size() >= suffix.size() && 0 == str.compare(str.size()-suffix.size(), suffix.size(), suffix);
}

LevelVector ProcessGroupWorker::receiveLevalAndBroadcast(){
  const auto dim = combiParameters_.getDim();

  // combine must have been called before this function
  assert(combinedUniDSGVector_.size() != 0 && "you must combine before you can eval");

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
  std::string filename = receiveStringFromManagerAndBroadcastToGroup();

  for (int g = 0; g < numGrids; g++) {  // loop over all grids and plot them
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
    auto pos = filename.find(".");
    if (pos != std::string::npos){
      // if filename contains ".", insert grid number before that
      filename.insert(pos, "_" + std::to_string(g));
    }
    dfg.writePlotFile(filename.c_str());
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
  for (Task* t : tasks_) {
    auto coeff = this->combiParameters_.getCoeff(t->getID());
    for (size_t i = 0; i < numCoordinates; ++i) {
      values[i] += t->getDistributedFullGrid().evalLocal(interpolationCoords[i]) * coeff;
    }
  }
  // reduce interpolated values within process group
  MPI_Allreduce(MPI_IN_PLACE, values.data(), static_cast<int>(numCoordinates),
                abstraction::getMPIDatatype(abstraction::getabstractionDataType<CombiDataType>()),
                MPI_SUM, theMPISystem()->getLocalComm());

  // need to reduce across process groups too
  // these do not strictly need to be allreduce (could be reduce), but it is easier to maintain that
  // way (all processes end up with valid values)
  MPI_Allreduce(MPI_IN_PLACE, values.data(), static_cast<int>(numCoordinates),
                abstraction::getMPIDatatype(abstraction::getabstractionDataType<CombiDataType>()),
                MPI_SUM, theMPISystem()->getGlobalReduceComm());

  return values;
}

template <typename T>
void writeValuesToH5File(const T& values, const std::string& fileName, const std::string& groupName,
                         const std::string& dataSetName, combigrid::real simulationTime) {
#ifdef HAVE_HIGHFIVE
  // check if file already exists, if no, create
  // TODO maybe use overwrite?
  HighFive::File h5_file(fileName, HighFive::File::OpenOrCreate | HighFive::File::ReadWrite);

  // std::vector<std::string> elems = h5_file.listObjectNames();
  // std::cout << elems << std::endl;

  HighFive::Group group;
  if (h5_file.exist(groupName)) {
    group = h5_file.getGroup(groupName);
  } else {
    group = h5_file.createGroup(groupName);
  }

  HighFive::DataSet dataset =
      group.createDataSet<CombiDataType>(dataSetName, HighFive::DataSpace::From(values));
  dataset.write(values);

  HighFive::Attribute aTime = dataset.createAttribute<combigrid::real>(
      "simulation_time", HighFive::DataSpace::From(simulationTime));
  aTime.write(simulationTime);
#else  // if not compiled with hdf5
  throw std::runtime_error("requesting hdf5 write but built without hdf5 support");
#endif
}

void ProcessGroupWorker::writeInterpolatedValuesPerGrid() {
  assert(combiParameters_.getNumGrids() == 1 && "interpolate only implemented for 1 species!");
  // receive coordinates and broadcast to group members
  std::vector<std::vector<real>> interpolationCoords;
  receiveAndBroadcastInterpolationCoords(interpolationCoords, combiParameters_.getDim());
  auto numCoordinates = interpolationCoords.size();

  // call interpolation function on tasks and write out task-wise
  for (size_t i = 0; i < tasks_.size(); ++i) {
    auto taskVals = tasks_[i]->getDistributedFullGrid().getInterpolatedValues(interpolationCoords);
    // cycle through ranks to write
    if (i % (theMPISystem()->getNumProcs()) == theMPISystem()->getLocalRank()) {
      // generate a rank-local per-run random number
      // std::random_device dev;
      static std::mt19937 rng(std::chrono::high_resolution_clock::now().time_since_epoch().count());
      static std::uniform_int_distribution<std::mt19937::result_type> dist(
          1, std::numeric_limits<size_t>::max());
      static size_t rankLocalRandom = dist(rng);

      std::string saveFilePath = "interpolated_task_" + std::to_string(tasks_[i]->getID()) + ".h5";
      std::string groupName = "run_" + std::to_string(rankLocalRandom);
      std::string datasetName = "interpolated_" + std::to_string(currentCombi_);
      writeValuesToH5File(taskVals, saveFilePath, groupName, datasetName,
                          tasks_[i]->getCurrentTime());
    }
  }
}

void ProcessGroupWorker::writeInterpolatedValues(const std::vector<CombiDataType>& values,
                                                 const std::string& valuesWriteFilename) {
  assert(combiParameters_.getNumGrids() == 1 && "interpolate only implemented for 1 species!");
  std::string groupName = "all_grids";
  std::string datasetName = "interpolated_" + std::to_string(currentCombi_);
  writeValuesToH5File(values, valuesWriteFilename, groupName, datasetName,
                      tasks_[0]->getCurrentTime());
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

    MASTER_EXCLUSIVE_SECTION { fg_red.add(fg, combiParameters_.getCoeff(t->getID())); }
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
  Stats::startEvent("task init in worker");
  auto taskDecomposition = combigrid::downsampleDecomposition(
          combiParameters_.getDecomposition(),
          combiParameters_.getLMax(), currentTask_->getLevelVector(),
          combiParameters_.getBoundary());
  currentTask_->init(theMPISystem()->getLocalComm(), taskDecomposition);
  if (ENABLE_FT) {
    t_fault_ = currentTask_->initFaults(t_fault_, startTimeIteration_);
  }
  Stats::stopEvent("task init in worker");
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
  Stats::startEvent("copyDataFromDSGtoDFG");
  for (Task* taskToUpdate : tasks_) {
    for (int g = 0; g < numGrids; g++) {
      // fill dfg with hierarchical coefficients from distributed sparse grid
      taskToUpdate->getDistributedFullGrid(g).extractFromUniformSG(*combinedUniDSGVector_[g]);
    }
  }
  Stats::stopEvent("copyDataFromDSGtoDFG");

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

void ProcessGroupWorker::zeroDsgsData() {
  for (auto& dsg : combinedUniDSGVector_)
    dsg->setZero();
}

/** free dsgus space */
void ProcessGroupWorker::deleteDsgsData() {
  for (auto& dsg : combinedUniDSGVector_)
    dsg->deleteSubspaceData();
}

void ProcessGroupWorker::writeDSGsToDisk(std::string filenamePrefix) {
  for (size_t i = 0; i < combinedUniDSGVector_.size(); ++i) {
    auto uniDsg = combinedUniDSGVector_[i].get();
    auto dsgToUse = uniDsg;
    dsgToUse->writeOneFileToDisk(filenamePrefix + "_" + std::to_string(i));
  }
}

void ProcessGroupWorker::readDSGsFromDisk(std::string filenamePrefix) {
  for (size_t i = 0; i < combinedUniDSGVector_.size(); ++i) {
    auto uniDsg = combinedUniDSGVector_[i].get();
    auto dsgToUse = uniDsg;
    dsgToUse->readOneFileFromDisk(filenamePrefix + "_" + std::to_string(i));
  }
}
} /* namespace combigrid */
