#include "sgpp/distributedcombigrid/manager/ProcessGroupWorker.hpp"

#include "boost/lexical_cast.hpp"

#include "sgpp/distributedcombigrid/combicom/CombiCom.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/hierarchization/DistributedHierarchization.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupSignals.hpp"
#include "sgpp/distributedcombigrid/mpi/MPIUtils.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/DistributedSparseGrid.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/DistributedSparseGridUniform.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LearningLoadModel.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"


#include <algorithm>
#include <iostream>
#include <string>
#include "sgpp/distributedcombigrid/mpi_fault_simulator/MPI-FT.h"

namespace combigrid {

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

ProcessGroupWorker::~ProcessGroupWorker() {}

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
      initializeTaskAndFaults();

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
      initializeTaskAndFaults();

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
      initializeTaskAndFaults();
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
    case RESCHEDULE_ADD_TASK: {
      assert(currentTask_ == nullptr);

      initializeTaskAndFaults(); // receive and initalize new task
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
        Stats::startEvent("worker run");
        currentTask_->run(theMPISystem()->getLocalComm());
        Stats::Event e = Stats::stopEvent("worker run");

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
/* not supported anymore
void ProcessGroupWorker::combine() {
  assert( false && "not properly implemented" );

  // early exit if no tasks available
  // todo: doesnt work, each pgrouproot must call reduce function
  assert(tasks_.size() > 0);

  assert( combiParametersSet_ );
  DimType dim = combiParameters_.getDim();
  const LevelVector& lmin = combiParameters_.getLMin();
  const LevelVector& lmax = combiParameters_.getLMax();
  const std::vector<bool>& boundary = combiParameters_.getBoundary();

  // erzeug dsg
  DistributedSparseGrid<CombiDataType> dsg(dim, lmax, lmin, boundary,
theMPISystem()->getLocalComm());

  for (Task* t : tasks_) {
    DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid();

    // hierarchize dfg
    DistributedHierarchization::hierarchize<CombiDataType>(dfg);

    // lokales reduce auf sg ->
    //CombiCom::distributedLocalReduce<CombiDataType>( dfg, dsg, combiParameters_.getCoeff(
t->getID() ) );
  }

  // globales reduce
  CombiCom::distributedGlobalReduce(dsg);

  for (Task* t : tasks_) {
    // get handle to dfg
    DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid();

    // lokales scatter von dsg auf dfg
    //CombiCom::distributedLocalScatter<CombiDataType>( dfg, dsg );

    // dehierarchize dfg
    DistributedHierarchization::dehierarchize<CombiDataType>(dfg);
  }
} */




/**
 * This method reduces the lmax and lmin vectors of the sparse grid according to the reduction
 * specifications in ctparam. It is taken care of that lmin does not fall below 1 and lmax >= lmin.
 * We do not reduce the levels in the last combination as we do not want to lose any information
 * for the final checkpoint.
 */
void reduceSparseGridCoefficients(LevelVector& lmax, LevelVector& lmin,
                                  IndexType totalNumberOfCombis, IndexType currentCombi,
                                  LevelVector reduceLmin, LevelVector reduceLmax) {
  // checking for valid combi step
  assert(currentCombi >= 0 && totalNumberOfCombis >= 0 && currentCombi <= totalNumberOfCombis);

  // this if-clause is currently always true, as initCombinedUniDSGVector is called only once,
  // at the beginning of the computation.
  // Leaving it here, in case the SG subspaces need to be re-initialized at some point.
  if (currentCombi < totalNumberOfCombis - 1) {  // do not reduce in last iteration
    for (size_t i = 0; i < reduceLmin.size(); ++i) {
      assert(reduceLmax[i] >= 0 && reduceLmin[i] >= 0);  // check for valid reduce values
      if (lmin[i] > 1) {
        lmin[i] = std::max((IndexType)1, lmin[i] - reduceLmin[i]);
      }
    }
    for (size_t i = 0; i < reduceLmax.size(); ++i) {
      lmax[i] = std::max(lmin[i], lmax[i] - reduceLmax[i]);
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

  // get all subspaces in the (optimized) combischeme
  SGrid<real> sg(dim, lmax, lmin, boundary);
  std::vector<LevelVector> subspaces;
  for (size_t ssID = 0; ssID < sg.getSize(); ++ssID) {
    const LevelVector& ss = sg.getLevelVector(ssID);
    subspaces.push_back(ss);
  }

  // create dsgs
  combinedUniDSGVector_.resize((size_t) numGrids);
  for (auto& uniDSG : combinedUniDSGVector_) {
    uniDSG = std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>(
        new DistributedSparseGridUniform<CombiDataType>(dim, subspaces, boundary,
                                                        theMPISystem()->getLocalComm()));
#ifdef DEBUG_OUTPUT
    MASTER_EXCLUSIVE_SECTION {
      std::cout << "dsg size: " << uniDSG->getRawDataSize() << " * " << sizeof(CombiDataType) << std::endl;
    }
#endif // def DEBUG_OUTPUT
  }

  // set subspace sizes local and global

  // register dsgs in all dfgs
  for (Task* t : tasks_) {
    for (int g = 0; g < numGrids; g++) {
#ifdef DEBUG_OUTPUT
      MASTER_EXCLUSIVE_SECTION {
        std::cout << "register task " << t->getID() << std::endl;
      }
#endif // def DEBUG_OUTPUT
      DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid(g);
      dfg.registerUniformSG(*(combinedUniDSGVector_[(size_t) g]));
    }
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

  for (IndexType g = 0; g < numGrids; g++) {
    CombiCom::distributedGlobalReduce(*combinedUniDSGVector_[g]);
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

  Stats::startEvent("combine local reduce");
  addFullGridsToUniformSG();
  Stats::stopEvent("combine local reduce");

  Stats::startEvent("combine global reduce");
  reduceUniformSG();
  Stats::stopEvent("combine global reduce");
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
  const int dim = static_cast<int>(combiParameters_.getDim());

  // combine must have been called before this function
  assert(combinedUniDSGVector_.size() != 0 && "you must combine before you can eval");

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

  for (int g = 0; g < numGrids; g++) {  // loop over all grids and plot them
    // create dfg
    bool forwardDecomposition = combiParameters_.getForwardDecomposition();
    DistributedFullGrid<CombiDataType> dfg(
      dim, leval, theMPISystem()->getLocalComm(), combiParameters_.getBoundary(),
      combiParameters_.getParallelization(), forwardDecomposition);
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

  DistributedFullGrid<CombiDataType> dfg(
      dim, leval, theMPISystem()->getLocalComm(), combiParameters_.getBoundary(),
      combiParameters_.getParallelization(), forwardDecomposition);

  this->fillDFGFromDSGU(dfg, 0);

  sendEvalNorms(dfg);
}

void ProcessGroupWorker::evalAnalyticalOnDFG() {
  auto leval = receiveLevalAndBroadcast();
  const int dim = static_cast<int>(leval.size());
  bool forwardDecomposition = combiParameters_.getForwardDecomposition();

  DistributedFullGrid<CombiDataType> dfg(
      dim, leval, theMPISystem()->getLocalComm(), combiParameters_.getBoundary(),
      combiParameters_.getParallelization(), forwardDecomposition);

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

  DistributedFullGrid<CombiDataType> dfg(
      dim, leval, theMPISystem()->getLocalComm(), combiParameters_.getBoundary(),
      combiParameters_.getParallelization(), forwardDecomposition);

  this->fillDFGFromDSGU(dfg, 0);
  // interpolate Task's analyticalSolution
  for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
    std::vector<double> coords(leval.size());
    dfg.getCoordsLocal(li, coords);

    dfg.getData()[li] -= tasks_[0]->analyticalSolution(coords, 0);
  }

  sendEvalNorms(dfg);
}

std::vector<CombiDataType> ProcessGroupWorker::interpolateValues() {
  assert(combiParameters_.getNumGrids() == 1 && "interpolate only implemented for 1 species!");
  // receive coordinates and broadcast to group members
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
  const int dim = static_cast<int>(combiParameters_.getDim());
  auto numCoordinates = coordsSize / dim;
  assert( coordsSize % dim == 0);
  std::vector<std::vector<real>> interpolationCoords(numCoordinates);
  auto it = interpolationCoordsSerial.cbegin();
  for (auto& coord: interpolationCoords) {
    coord.insert(coord.end(), it, it+dim);
    it += dim;
  }
  assert(it == interpolationCoordsSerial.end());

  // call interpolation function on tasks and reduce with combination coefficient
  std::vector<CombiDataType> values(numCoordinates, 0.);
  for (Task* t : tasks_){
    auto coeff = combiParameters_.getCoeff(t->getID());
    auto taskVals = t->getDistributedFullGrid().getInterpolatedValues(interpolationCoords);

    for (size_t i = 0; i < numCoordinates; ++i) {
      values[i] += taskVals[i] * coeff;
    }
  }
  return values;
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

void ProcessGroupWorker::initializeTaskAndFaults(bool mayAlreadyExist /*=true*/) {
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

  // add task to task storage
  tasks_.push_back(t);

  status_ = PROCESS_GROUP_BUSY;

  // set currentTask
  currentTask_ = tasks_.back();

  // initalize task
  Stats::startEvent("task init in worker");
  currentTask_->init(theMPISystem()->getLocalComm());
  t_fault_ = currentTask_->initFaults(t_fault_, startTimeIteration_);
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

  int numGrids = static_cast<int>(combiParameters_
                     .getNumGrids());  // we assume here that every task has the same number of grids

  for (int g = 0; g < numGrids; g++) {
    assert(combinedUniDSGVector_[g] != NULL);

    // get handle to dfg
    DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid(g);

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

void ProcessGroupWorker::zeroDsgsData() {
  for (auto& dsg : combinedUniDSGVector_)
    dsg->setZero();
}

/** free dsgus space */
void ProcessGroupWorker::deleteDsgsData() {
  for (auto& dsg : combinedUniDSGVector_)
    dsg->deleteSubspaceData();
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
