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
      Stats::startEvent("combine");
      combineUniform();
      currentCombi_++;
      Stats::stopEvent("combine");

    } break;
    case COMBINE_LOCAL_AND_GLOBAL: {
      Stats::startEvent("combineLocalAndGlobal");
      combineLocalAndGlobal();
      Stats::stopEvent("combineLocalAndGlobal");

    } break;
    case COMBINE_THIRD_LEVEL: {
      Stats::startEvent("combineThirdLevel");
      combineThirdLevel();
      currentCombi_++;
      Stats::stopEvent("combineThirdLevel");

    } break;
    case WAIT_FOR_TL_COMBI_RESULT: {
      Stats::startEvent("waitForThirdLevelCombiResult");
      waitForThirdLevelCombiResult();
      currentCombi_++;
      Stats::stopEvent("waitForThirdLevelCombiResult");

    } break;
    case REDUCE_SUBSPACE_SIZES_TL_AND_ALLOCATE_EXTRA_SG: {
      Stats::startEvent("unifySubspaceSizesThirdLevelExtra");
      reduceSubspaceSizesThirdLevel(true);
      Stats::stopEvent("unifySubspaceSizesThirdLevelExtra");

    } break;
    case REDUCE_SUBSPACE_SIZES_TL: {
      Stats::startEvent("unifySubspaceSizesThirdLevel");
      reduceSubspaceSizesThirdLevel(false);
      Stats::stopEvent("unifySubspaceSizesThirdLevel");

    } break;
    case WAIT_FOR_TL_SIZE_UPDATE: {
      Stats::startEvent("waitForThirdLevelSizeUpdate");
      waitForThirdLevelSizeUpdate();
      Stats::stopEvent("waitForThirdLevelSizeUpdate");

    } break;
    case WRITE_DFGS_TO_VTK: {
      Stats::startEvent("writeVTKPlotFilesOfAllTasks");
      writeVTKPlotFilesOfAllTasks();
      Stats::stopEvent("writeVTKPlotFilesOfAllTasks");
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
        setCombinedSolutionUniform(currentTask_);
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
      Stats::startEvent("interpolate values");
      auto values = interpolateValues();
      // send result
      MASTER_EXCLUSIVE_SECTION {
        MPI_Send(values.data(), values.size(), abstraction::getMPIDatatype(
          abstraction::getabstractionDataType<CombiDataType>()), theMPISystem()->getManagerRank(),
          TRANSFER_INTERPOLATION_TAG, theMPISystem()->getGlobalComm());
      }
      Stats::stopEvent("interpolate values");
    } break;
    case WRITE_INTERPOLATED_VALUES_PER_GRID: {  // interpolate values on given coordinates and write values to .h5
      Stats::startEvent("write interpolated values");
      writeInterpolatedValuesPerGrid();
      Stats::stopEvent("write interpolated values");
    } break;
    case RESCHEDULE_ADD_TASK: {
      assert(currentTask_ == nullptr);

      receiveAndInitializeTaskAndFaults(); // receive and initalize new task
			// now the variable currentTask_ contains the newly received task
      currentTask_->setZero();
      updateTaskWithCurrentValues(*currentTask_, combiParameters_.getNumGrids());
      currentTask_->setFinished(true);
      currentTask_ = nullptr;
		} break;
    case RESCHEDULE_REMOVE_TASK: {
      assert(currentTask_ == nullptr);

      int taskID;
      MASTER_EXCLUSIVE_SECTION {
        MPI_Recv(&taskID, 1, MPI_INT, theMPISystem()->getManagerRank(), 0,
                 theMPISystem()->getGlobalComm(), MPI_STATUS_IGNORE);
      }
      MPI_Bcast(&taskID, 1, MPI_INT, theMPISystem()->getMasterRank(),
                theMPISystem()->getLocalComm());

      // search for task send to group master and remove
      for(size_t i=0; i < tasks_.size(); ++i) {
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
  auto numGrids = combiParameters_.getNumGrids();
  DimType dim = combiParameters_.getDim();
  LevelVector lmin = combiParameters_.getLMin();
  LevelVector lmax = combiParameters_.getLMax();
  const std::vector<bool>& boundary = combiParameters_.getBoundary();

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
  combinedUniDSGVector_.resize((size_t) numGrids);
  for (auto& uniDSG : combinedUniDSGVector_) {
    uniDSG = std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>(
        new DistributedSparseGridUniform<CombiDataType>(dim, lmax, lmin, boundary,
                                                        theMPISystem()->getLocalComm()));
#ifdef DEBUG_OUTPUT
    MASTER_EXCLUSIVE_SECTION {
      std::cout << "dsg size: " << uniDSG->getRawDataSize() << " * " << sizeof(CombiDataType) << std::endl;
    }
#endif // def DEBUG_OUTPUT
  }

  // set subspace sizes local and global

  // register dsgs in all dfgs
  for (int g = 0; g < numGrids; g++) {
    for (Task* t : tasks_) {
#ifdef DEBUG_OUTPUT
      MASTER_EXCLUSIVE_SECTION {
        std::cout << "register task " << t->getID() << std::endl;
      }
#endif // def DEBUG_OUTPUT
      DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid(g);
      dfg.registerUniformSG(*(combinedUniDSGVector_[(size_t) g]));
    }
    // // we may clear the levels_ member of the sparse grids here to save memory
    // // but only if we need no new full grids initialized from the sparse grids!
    // // ...such as for rescheduling or interpolation (parallelEval/ evalNorm / ...)
    // combinedUniDSGVector_[(size_t) g]->resetLevels();
  }

  // global reduce of subspace sizes
  CommunicatorType globalReduceComm = theMPISystem()->getGlobalReduceComm();
  for (int g = 0; g < numGrids; g++) {
    reduceSubspaceSizes(combinedUniDSGVector_[(size_t)g].get(), globalReduceComm);
#ifdef DEBUG_OUTPUT
    MASTER_EXCLUSIVE_SECTION {
      std::cout << "dsg size: " << combinedUniDSGVector_[(size_t)g]->getRawDataSize() << " * "
                << sizeof(CombiDataType) << std::endl;
    }
#endif  // def DEBUG_OUTPUT
  }
}

void ProcessGroupWorker::hierarchizeFullGrids() {
  auto numGrids = combiParameters_.getNumGrids();
  //real localMax(0.0);
  // std::vector<CombiDataType> beforeCombi;
  for (Task* t : tasks_) {
    for (IndexType g = 0; g < numGrids; g++) {
      DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid(static_cast<int>(g));

      // hierarchize dfg
      DistributedHierarchization::hierarchize<CombiDataType>(
          dfg, combiParameters_.getHierarchizationDims(), combiParameters_.getHierarchicalBases());
    }
  }
}

void ProcessGroupWorker::addFullGridsToUniformSG() {
  assert(combinedUniDSGVector_.size() > 0 && "Initialize dsgu first with "
                                             "initCombinedUniDSGVector()");
  auto numGrids = combiParameters_.getNumGrids();
  for (Task* t : tasks_) {
    for (IndexType g = 0; g < numGrids; g++) {
      DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid(static_cast<int>(g));

      // lokales reduce auf sg ->
      dfg.addToUniformSG(*combinedUniDSGVector_[g], combiParameters_.getCoeff(t->getID()));
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
  //Stats::startEvent("combine zeroDsgsData");
  zeroDsgsData();
  //Stats::stopEvent("combine zeroDsgsData");

  //Stats::startEvent("combine hierarchize");
  hierarchizeFullGrids();
  //Stats::stopEvent("combine hierarchize");

#ifdef DEBUG_OUTPUT
  MASTER_EXCLUSIVE_SECTION { std::cout << "mid combining \n"; }
#endif

  //Stats::startEvent("combine local reduce");
  addFullGridsToUniformSG();
  //Stats::stopEvent("combine local reduce");

#ifdef DEBUG_OUTPUT
  MASTER_EXCLUSIVE_SECTION { std::cout << "almost done combining \n"; }
#endif

  //Stats::startEvent("combine global reduce");
  reduceUniformSG();
  //Stats::stopEvent("combine global reduce");

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
  const int dim = static_cast<int>(combiParameters_.getDim());

  // receive leval and broadcast to group members
  std::vector<int> tmp(dim);
  MASTER_EXCLUSIVE_SECTION {
    MPI_Recv(&tmp[0], dim, MPI_INT, theMPISystem()->getManagerRank(), TRANSFER_LEVAL_TAG,
             theMPISystem()->getGlobalComm(), MPI_STATUS_IGNORE);
  }

  MPI_Bcast(&tmp[0], dim, MPI_INT, theMPISystem()->getMasterRank(), theMPISystem()->getLocalComm());
  LevelVector leval(tmp.begin(), tmp.end());
  return leval;
}

void ProcessGroupWorker::fillDFGFromDSGU(DistributedFullGrid<CombiDataType>& dfg, IndexType g) {
  DistributedHierarchization::fillDFGFromDSGU(dfg, *combinedUniDSGVector_[g],
                                              combiParameters_.getHierarchizationDims(),
                                              combiParameters_.getHierarchicalBases());
}

void ProcessGroupWorker::parallelEvalUniform() {
  assert(uniformDecomposition);

  assert(combiParametersSet_);
  auto numGrids = combiParameters_.getNumGrids();  // we assume here that every task has the same number of grids

  auto leval = receiveLevalAndBroadcast();
  const int dim = static_cast<int>(leval.size());

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
  const int dim = static_cast<int>(leval.size());
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
  const int dim = static_cast<int>(leval.size());
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
  const int dim = static_cast<int>(leval.size());
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
  int taskID;
  MASTER_EXCLUSIVE_SECTION {
    MPI_Recv(&taskID, 1, MPI_INT, theMPISystem()->getManagerRank(), 0,
             theMPISystem()->getGlobalComm(), MPI_STATUS_IGNORE);
  }
  MPI_Bcast(&taskID, 1, MPI_INT, theMPISystem()->getMasterRank(), theMPISystem()->getLocalComm());

  // call diagnostics on that Task
  for (auto task : tasks_) {
    if (task->getID() == taskID) {
      std::vector<DistributedSparseGridUniform<CombiDataType>*> dsgsToPassToTask;
      for (auto& dsgPtr : combinedUniDSGVector_){
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
    MPI_Recv(interpolationCoordsSerial.data(), coordsSize, realType, theMPISystem()->getManagerRank(),
             TRANSFER_INTERPOLATION_TAG, theMPISystem()->getGlobalComm(), MPI_STATUS_IGNORE);
  }
  // broadcast size of vector, and vector
  MPI_Bcast(&coordsSize, 1, MPI_INT, theMPISystem()->getMasterRank(), theMPISystem()->getLocalComm());
  interpolationCoordsSerial.resize(coordsSize);
  MPI_Bcast(interpolationCoordsSerial.data(), coordsSize, realType,
            theMPISystem()->getMasterRank(), theMPISystem()->getLocalComm());

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
    auto coeff = this->combiParameters_.getCoeff(t->getID());
    auto taskVals = t->getDistributedFullGrid().getInterpolatedValues(interpolationCoords);

    for (size_t i = 0; i < numCoordinates; ++i) {
      values[i] += taskVals[i] * coeff;
    }
  }
  return values;
}


void ProcessGroupWorker::writeInterpolatedValuesPerGrid() {
#ifdef HAVE_HIGHFIVE
  assert(combiParameters_.getNumGrids() == 1 && "interpolate only implemented for 1 species!");
  // receive coordinates and broadcast to group members
  std::vector<std::vector<real>> interpolationCoords;
  receiveAndBroadcastInterpolationCoords(interpolationCoords, combiParameters_.getDim());
  auto numCoordinates = interpolationCoords.size();

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
  const std::vector<bool>& boundary = combiParameters_.getBoundary();
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
  CombiParameters tmp;

  // local root receives task
  MASTER_EXCLUSIVE_SECTION {
    MPIUtils::receiveClass(&tmp, theMPISystem()->getManagerRank(), theMPISystem()->getGlobalComm());
    //std::cout << "master received combiparameters \n";
  }

  // broadcast task to other process of pgroup
  MPIUtils::broadcastClass(&tmp, theMPISystem()->getMasterRank(), theMPISystem()->getLocalComm());
  //std::cout << "worker received combiparameters \n";

  combiParameters_ = tmp;
  combiParametersSet_ = true;

  // overwrite local comm with cartesian communicator
  if (!isGENE  && tmp.isParallelizationSet()){
    // cf. https://www.rookiehpc.com/mpi/docs/mpi_cart_create.php
    // get decompositon from combi params
    auto par = combiParameters_.getParallelization();

    // important: note reverse ordering of dims! -- cf DistributedFullGrid //TODO(pollinta) remove reverse ordering
    std::vector<int> dims (par.rbegin(), par.rend());
    if (!reverseOrderingDFGPartitions) {
      dims.assign(par.begin(), par.end());
    }

    // Make all dimensions not periodic //TODO(pollinta) allow periodicity
    std::vector<int> periods (combiParameters_.getDim(), 0);

    // don't let MPI assign arbitrary ranks
    int reorder = false;

    // Create a communicator given the topology.
    MPI_Comm new_communicator;
    MPI_Cart_create(theMPISystem()->getLocalComm(), combiParameters_.getDim(), dims.data(), periods.data(), reorder, &new_communicator);

    theMPISystem()->storeLocalComm(new_communicator);
  }
}

void ProcessGroupWorker::setCombinedSolutionUniform(Task* t) {
  assert(combinedUniDSGVector_.size() != 0);
  assert(combiParametersSet_);

  auto numGrids = combiParameters_.getNumGrids();  // we assume here that every task has the same number of grids

  for (IndexType g = 0; g < numGrids; g++) {
    assert(combinedUniDSGVector_[g] != NULL);

    // get handle to dfg
    DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid(static_cast<int>(g));

    DistributedHierarchization::fillDFGFromDSGU(dfg, *combinedUniDSGVector_[g],
                                                combiParameters_.getHierarchizationDims(),
                                                combiParameters_.getHierarchicalBases());
  }
}

void ProcessGroupWorker::integrateCombinedSolution() {
  int numGrids = static_cast<int>(combiParameters_.getNumGrids());
  Stats::startEvent("integrateCombinedSolution");
  for (Task* t : tasks_)
    updateTaskWithCurrentValues(*t, numGrids);
  Stats::stopEvent("integrateCombinedSolution");
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
      // copy partial data from uniDSG to extraDSG
      for (size_t i = 0; i < uniDsg->getNumSubspaces(); ++i) {
        assert(dsgToUse->getDataSize(i) == 0 ||
               dsgToUse->getDataSize(i) == uniDsg->getDataSize(i));
        std::copy_n(uniDsg->getData(i), dsgToUse->getDataSize(i), dsgToUse->getData(i));
      }
    }

    // send dsg data to manager
    Stats::startEvent("combine send dsg data to manager");
    sendDsgData(dsgToUse, manager, managerComm);
    Stats::stopEvent("combine send dsg data to manager");

    // recv combined dsgu from manager
    Stats::startEvent("combine recv combined data from manager");
    recvDsgData(dsgToUse, manager, managerComm);
    Stats::stopEvent("combine recv combined data from manager");

    if (extraUniDSGVector_.size() > 0) {
      // copy partial data from extraDSG back to uniDSG
      for (size_t i = 0; i < uniDsg->getNumSubspaces(); ++i) {
        std::copy_n(dsgToUse->getData(i), dsgToUse->getDataSize(i), uniDsg->getData(i));
      }
    }

    // distribute solution in globalReduceComm to other pgs
    MPI_Request request = asyncBcastDsgData(uniDsg, globalReduceRank, globalReduceComm);
    requests.push_back(request);
  }
  // update fgs
  integrateCombinedSolution();

  // wait for bcasts to other pgs in globalReduceComm
  Stats::startEvent("combine wait for async bcasts");
  for (MPI_Request& request : requests)
    MPI_Wait(&request, MPI_STATUS_IGNORE);
  Stats::stopEvent("combine wait for async bcasts");
}

/** Reduces subspace sizes with remote.
 */
void ProcessGroupWorker::reduceSubspaceSizesThirdLevel(bool thirdLevelExtraSparseGrid) {
  assert(combiParametersSet_);

  // prepare a buffer with enough space for all subspace sizes from all dsgs
  int bufSize = 0;
  for (auto& dsg : combinedUniDSGVector_)
    bufSize += static_cast<int>(dsg->getSubspaceDataSizes().size());
  std::vector<size_t> buff;
  buff.reserve((size_t)bufSize);

  // fill the buffer with the subspace sizes
  for (auto& dsg : combinedUniDSGVector_) {
    for(size_t size : dsg->getSubspaceDataSizes()) {
      buff.push_back(size);
    }
  }

  // prepare for MPI calls to manager
  CommunicatorType thirdLevelComm = theMPISystem()->getThirdLevelComms()[0];
  RankType thirdLevelManagerRank = theMPISystem()->getThirdLevelManagerRank();

  // send size of buffer to manager
  MPI_Gather(&bufSize, 1, MPI_INT, nullptr, 0, MPI_INT, thirdLevelManagerRank,
             thirdLevelComm);

  // send subspace sizes to manager
  MPI_Datatype dtype = getMPIDatatype(
                        abstraction::getabstractionDataType<size_t>());
  MPI_Gatherv(buff.data(), bufSize, dtype, nullptr, nullptr,
              nullptr, dtype, thirdLevelManagerRank, thirdLevelComm);

  if (thirdLevelExtraSparseGrid && extraUniDSGVector_.empty()) {
    // create new vector for extra sparse grids (that will be only on this process group)
    extraUniDSGVector_.resize(combinedUniDSGVector_.size());
    for (auto& extraUniDSG : extraUniDSGVector_) {
      extraUniDSG = std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>(
          new DistributedSparseGridUniform<CombiDataType>(
              combinedUniDSGVector_[0]->getDim(), combinedUniDSGVector_[0]->getAllLevelVectors(),
              combinedUniDSGVector_[0]->getBoundaryVector(), theMPISystem()->getLocalComm()));
    }
  }

  // receive updated sizes from manager
  MPI_Scatterv(nullptr, 0, nullptr, dtype,
               buff.data(), bufSize, dtype, thirdLevelManagerRank, thirdLevelComm);

  if (!thirdLevelExtraSparseGrid) {
    // distribute updated sizes to workers with same decomposition (global reduce comm)
    // cf. waitForThirdLevelSizeUpdate(), which is called in other process groups
    CommunicatorType globalReduceComm = theMPISystem()->getGlobalReduceComm();
    RankType globalReduceRank = theMPISystem()->getGlobalReduceRank();
    MPI_Bcast(buff.data(), bufSize, dtype, globalReduceRank, globalReduceComm);
  }

  // update either old or new sparse grids
  auto* uniDSGVectorToSet = &combinedUniDSGVector_;
  if (thirdLevelExtraSparseGrid) {
    uniDSGVectorToSet = &extraUniDSGVector_;
  }
  // set updated sizes in dsgs
  std::vector<size_t>::iterator buffIt = buff.begin();
  for (auto& dsg : *uniDSGVectorToSet) {
    for(size_t i = 0; i < dsg->getNumSubspaces(); ++i) {
      dsg->setDataSize(i, *(buffIt++));
    }
    dsg->createSubspaceData();
#ifdef DEBUG_OUTPUT
    MASTER_EXCLUSIVE_SECTION {
      std::cout << "third level reduce dsg size: " << dsg->getRawDataSize() << " * "
                << sizeof(CombiDataType) << std::endl;
    }
#endif  // def DEBUG_OUTPUT
  }
}

void ProcessGroupWorker::waitForThirdLevelSizeUpdate() {
  // prepare a buffer with enough space for all subspace sizes from all dsgs
  int buffSize = 0;
  for (auto& dsg : combinedUniDSGVector_)
    buffSize += static_cast<int>(dsg->getSubspaceDataSizes().size());
  std::vector<size_t> buff((size_t) buffSize);

  // receive updated sizes from third level pgroup (global reduce comm)
  RankType thirdLevelPG = (RankType) combiParameters_.getThirdLevelPG();
  CommunicatorType globalReduceComm = theMPISystem()->getGlobalReduceComm();
  MPI_Datatype dtype = getMPIDatatype(
                        abstraction::getabstractionDataType<size_t>());
  MPI_Bcast(buff.data(), buffSize, dtype, thirdLevelPG, globalReduceComm);

  // set updated sizes in dsgs
  std::vector<size_t>::iterator buffIt = buff.begin();
  for (auto& dsg : combinedUniDSGVector_) {
    for(size_t i = 0; i < dsg->getNumSubspaces(); ++i) {
      dsg->setDataSize(i, *(buffIt++));
    }
  }
}

void ProcessGroupWorker::waitForThirdLevelCombiResult() {
  // receive third level combi result from third level pgroup (global reduce comm)
  RankType thirdLevelPG = (RankType) combiParameters_.getThirdLevelPG();
  CommunicatorType globalReduceComm = theMPISystem()->getGlobalReduceComm();

  for (auto& dsg : combinedUniDSGVector_) {
    MPI_Request request = asyncBcastDsgData(dsg.get(), thirdLevelPG, globalReduceComm);
    MPI_Wait(&request, MPI_STATUS_IGNORE);
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

void ProcessGroupWorker::updateTaskWithCurrentValues(Task& taskToUpdate, size_t numGrids) {
  for (int g = 0; g < numGrids; g++) {
    // get handle to dfg
    DistributedFullGrid<CombiDataType>& dfg = taskToUpdate.getDistributedFullGrid(g);

    DistributedHierarchization::fillDFGFromDSGU(dfg, *combinedUniDSGVector_[g],
                                                combiParameters_.getHierarchizationDims(),
                                                combiParameters_.getHierarchicalBases());

    // std::vector<CombiDataType> datavector(dfg.getElementVector());
    // afterCombi = datavector;
    // if exceeds normalization limit, normalize dfg with global max norm
    /*
    if( globalMax > 1000 ){
      dfg.mul( 1.0 / globalMax );
      std::cout << "normalized dfg with " << globalMax << std::endl;
    }
    */
  }
}

} /* namespace combigrid */
