#include "manager/ProcessGroupWorker.hpp"

#include "combicom/CombiCom.hpp"
#include "manager/CombiParameters.hpp"
#include "manager/ProcessGroupSignals.hpp"
#include "mpi/MPIUtils.hpp"
#include "sparsegrid/DistributedSparseGridIO.hpp"
#include "sparsegrid/DistributedSparseGridUniform.hpp"
#include "loadmodel/LearningLoadModel.hpp"
#include "mpi/MPISystem.hpp"
#include "mpi_fault_simulator/MPI-FT.h"
#include "io/H5InputOutput.hpp"
#include "utils/MonteCarlo.hpp"

#include "boost/lexical_cast.hpp"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <random>
#include <string>
#include <thread>

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

void sendNormsToManager(const std::vector<double> lpnorms) {
  // get Lp norm on every worker; reduce through dfg function
  for (int p = 0; p < lpnorms.size(); ++p) {
    // send from master to manager
    MASTER_EXCLUSIVE_SECTION {
      MPI_Send(&lpnorms[p], 1, MPI_DOUBLE, theMPISystem()->getManagerRank(), TRANSFER_NORM_TAG,
               theMPISystem()->getGlobalComm());
    }
  }
}

std::vector<std::vector<real>> receiveAndBroadcastInterpolationCoords(DimType dim) {
  std::vector<std::vector<real>> interpolationCoords;
  std::vector<real> interpolationCoordsSerial;
  auto realType = abstraction::getMPIDatatype(abstraction::getabstractionDataType<real>());
  int coordsSize = 0;
  MASTER_EXCLUSIVE_SECTION {
    MPI_Status status;
    status.MPI_ERROR = MPI_SUCCESS;
    int result = MPI_Probe(theMPISystem()->getManagerRank(), TRANSFER_INTERPOLATION_TAG,
                           theMPISystem()->getGlobalComm(), &status);
#ifndef NDEBUG
    assert(result == MPI_SUCCESS);
    if (status.MPI_ERROR != MPI_SUCCESS) {
      std::string errorString;
      errorString.resize(10000);
      int errorStringLength;
      MPI_Error_string(status.MPI_ERROR, &errorString[0], &errorStringLength);
      errorString.resize(errorStringLength);
      std::cout << "error probe: " << errorString << std::endl;
    }
#endif  // NDEBUG
    result = MPI_Get_count(&status, realType, &coordsSize);
    assert(result == MPI_SUCCESS);
    assert(coordsSize > 0);

    // resize buffer to appropriate size and receive
    interpolationCoordsSerial.resize(coordsSize);
    result = MPI_Recv(interpolationCoordsSerial.data(), coordsSize, realType,
                      theMPISystem()->getManagerRank(), TRANSFER_INTERPOLATION_TAG,
                      theMPISystem()->getGlobalComm(), &status);
#ifndef NDEBUG
    assert(result == MPI_SUCCESS);
    if (status.MPI_ERROR != MPI_SUCCESS) {
      std::string errorString;
      errorString.resize(10000);
      int errorStringLength;
      MPI_Error_string(status.MPI_ERROR, &errorString[0], &errorStringLength);
      errorString.resize(errorStringLength);
      std::cout << "error recv: " << errorString << std::endl;
    }
    assert(status.MPI_ERROR == MPI_SUCCESS);
    for (const auto& coord : interpolationCoordsSerial) {
      assert(coord >= 0.0 && coord <= 1.0);
    }
#endif  // NDEBUG
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
  interpolationCoords = deserializeInterpolationCoords(interpolationCoordsSerial, dim);
  return interpolationCoords;
}

ProcessGroupWorker::ProcessGroupWorker()
    : status_(PROCESS_GROUP_WAIT),
      combinedUniDSGVector_(0),
      combiParameters_(),
      combiParametersSet_(false),
      currentCombi_(0) {}

ProcessGroupWorker::~ProcessGroupWorker() {}

LevelVector ProcessGroupWorker::receiveLevalAndBroadcast() {
  const auto dim = combiParameters_.getDim();

  // receive leval and broadcast to group members
  std::vector<int> tmp(dim);
  MASTER_EXCLUSIVE_SECTION {
    MPI_Recv(&tmp[0], static_cast<int>(dim), MPI_INT, theMPISystem()->getManagerRank(),
             TRANSFER_LEVAL_TAG, theMPISystem()->getGlobalComm(), MPI_STATUS_IGNORE);
  }

  MPI_Bcast(&tmp[0], dim, MPI_INT, theMPISystem()->getMasterRank(), theMPISystem()->getLocalComm());
  LevelVector leval(tmp.begin(), tmp.end());
  return leval;
}

SignalType ProcessGroupWorker::wait() {
  if (status_ == PROCESS_GROUP_FAIL) {  // in this case worker got reused
    status_ = PROCESS_GROUP_WAIT;
  }
  if (status_ != PROCESS_GROUP_WAIT) {
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
  // process signal
  switch (signal) {
    case RUN_FIRST: {
      receiveAndInitializeTask();
      status_ = PROCESS_GROUP_BUSY;

      // execute task
      Stats::startEvent("run first");
      auto& currentTask = this->getTaskWorker().getTasks().back();
      currentTask->run(theMPISystem()->getLocalComm());
      Stats::Event e = Stats::stopEvent("run first");
    } break;
    case RUN_NEXT: {
      assert(this->getTaskWorker().getTasks().size() > 0);
      // // free space for computation
      // deleteDsgsData();

      this->runAllTasks();

    } break;
    case ADD_TASK: {  // add a new task to the process group
      // initalize task and set values to zero
      // the task will get the proper initial solution during the next combine
      receiveAndInitializeTask();

      auto& currentTask = this->getTaskWorker().getTasks().back();
      currentTask->setZero();
      currentTask->setFinished(true);
    } break;
    case RESET_TASKS: {  // deleta all tasks (used in process recovery)
      std::cout << "resetting tasks" << std::endl;

      this->getTaskWorker().deleteTasks();
      status_ = PROCESS_GROUP_BUSY;
    } break;
    case EXIT: {
      this->exit();
    } break;
    case INIT_DSGUS: {
      Stats::startEvent("initialize dsgu");
      initCombinedUniDSGVector();
      Stats::stopEvent("initialize dsgu");

    } break;
    case COMBINE: {  // start combination
      combineUniform();
    } break;
    case COMBINE_LOCAL_AND_GLOBAL: {
      Stats::startEvent("combine local");
      combineLocalAndGlobal();
      Stats::stopEvent("combine local");

    } break;
    case COMBINE_THIRD_LEVEL: {
      Stats::startEvent("combine third level");
      combineThirdLevel();
      Stats::stopEvent("combine third level");
    } break;
    case COMBINE_WRITE_DSGS: {
      Stats::startEvent("combine third level write");
      combineThirdLevelFileBasedWrite(receiveStringFromManagerAndBroadcastToGroup(),
                                      receiveStringFromManagerAndBroadcastToGroup());
      Stats::stopEvent("combine third level write");
    } break;
    case COMBINE_READ_DSGS_AND_REDUCE: {
      Stats::startEvent("combine third level read");
      combineThirdLevelFileBasedReadReduce(receiveStringFromManagerAndBroadcastToGroup(),
                                           receiveStringFromManagerAndBroadcastToGroup());
      Stats::stopEvent("combine third level read");
    } break;
    case COMBINE_THIRD_LEVEL_FILE: {
      Stats::startEvent("combine third level file");
      combineThirdLevelFileBased(receiveStringFromManagerAndBroadcastToGroup(),
                                 receiveStringFromManagerAndBroadcastToGroup(),
                                 receiveStringFromManagerAndBroadcastToGroup(),
                                 receiveStringFromManagerAndBroadcastToGroup());
      Stats::stopEvent("combine third level file");
    } break;
    case WAIT_FOR_TL_COMBI_RESULT: {
      Stats::startEvent("wait third level result");
      waitForThirdLevelCombiResult();
      Stats::stopEvent("wait third level result");

    } break;
    case REDUCE_SUBSPACE_SIZES_TL_AND_ALLOCATE_EXTRA_SG: {
      Stats::startEvent("unify extra third level");
      reduceSubspaceSizesThirdLevel(true);
      Stats::stopEvent("unify extra third level");

    } break;
    case REDUCE_SUBSPACE_SIZES_TL: {
      Stats::startEvent("unify sizes third level");
      reduceSubspaceSizesThirdLevel(false);
      Stats::stopEvent("unify sizes third level");

    } break;
    case WAIT_FOR_TL_SIZE_UPDATE: {
      Stats::startEvent("wait third level size");
      waitForThirdLevelSizeUpdate();
      Stats::stopEvent("wait third level size");

    } break;
    case WRITE_DFGS_TO_VTK: {
      Stats::startEvent("write vtk all tasks");
      writeVTKPlotFilesOfAllTasks();
      Stats::stopEvent("write vtk all tasks");
    } break;
    case WRITE_DSGS_TO_DISK: {
      Stats::startEvent("write to disk");
      std::string filenamePrefix = receiveStringFromManagerAndBroadcastToGroup();
      writeDSGsToDisk(filenamePrefix);
      Stats::stopEvent("write to disk");
    } break;
    case READ_DSGS_FROM_DISK: {
      Stats::startEvent("read from disk");
      std::string filenamePrefix = receiveStringFromManagerAndBroadcastToGroup();
      readDSGsFromDisk(filenamePrefix);
      Stats::stopEvent("read from disk");
    } break;
    case UPDATE_COMBI_PARAMETERS: {  // update combiparameters (e.g. in case of faults -> FTCT)

      updateCombiParameters();

    } break;
    case RECOMPUTE: {  // recompute the received task (immediately computes tasks ->
                       // difference to ADD_TASK)
      receiveAndInitializeTask();
      auto& currentTask = this->getTaskWorker().getTasks().back();

      currentTask->setZero();
      // fill task with combisolution
      fillDFGFromDSGU(*currentTask);

      // execute task
      Stats::Event e = Stats::Event();
      currentTask->run(theMPISystem()->getLocalComm());
      e.end = std::chrono::high_resolution_clock::now();
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
      auto lpnorms = this->getLpNorms(2);
      // send from master to manager
      MASTER_EXCLUSIVE_SECTION {
        MPI_Send(lpnorms.data(), static_cast<int>(lpnorms.size()), MPI_DOUBLE,
                 theMPISystem()->getManagerRank(), TRANSFER_NORM_TAG,
                 theMPISystem()->getGlobalComm());
      }
      Stats::stopEvent("get L2 norm");
    } break;
    case GET_L1_NORM: {  // evaluate norm on dfgs and send
      Stats::startEvent("get L1 norm");
      auto lpnorms = this->getLpNorms(1);
      MASTER_EXCLUSIVE_SECTION {
        MPI_Send(lpnorms.data(), static_cast<int>(lpnorms.size()), MPI_DOUBLE,
                 theMPISystem()->getManagerRank(), TRANSFER_NORM_TAG,
                 theMPISystem()->getGlobalComm());
      }
      Stats::stopEvent("get L1 norm");
    } break;
    case GET_MAX_NORM: {  // evaluate norm on dfgs and send
      Stats::startEvent("get max norm");
      auto lpnorms = this->getLpNorms(0);
      MASTER_EXCLUSIVE_SECTION {
        MPI_Send(lpnorms.data(), static_cast<int>(lpnorms.size()), MPI_DOUBLE,
                 theMPISystem()->getManagerRank(), TRANSFER_NORM_TAG,
                 theMPISystem()->getGlobalComm());
      }
      Stats::stopEvent("get max norm");
    } break;
    case INTERPOLATE_VALUES: {  // interpolate values on given coordinates
      Stats::startEvent("interpolate values");
      auto values =
          interpolateValues(receiveAndBroadcastInterpolationCoords(combiParameters_.getDim()));
      Stats::stopEvent("interpolate values");
    } break;
    case INTERPOLATE_VALUES_AND_SEND_BACK: {
      Stats::startEvent("interpolate values");
      auto values =
          interpolateValues(receiveAndBroadcastInterpolationCoords(combiParameters_.getDim()));
      // send result
      MASTER_EXCLUSIVE_SECTION {
        MPI_Send(values.data(), values.size(),
                 abstraction::getMPIDatatype(abstraction::getabstractionDataType<CombiDataType>()),
                 theMPISystem()->getManagerRank(), TRANSFER_INTERPOLATION_TAG,
                 theMPISystem()->getGlobalComm());
      }
      Stats::stopEvent("interpolate values");
    } break;
    case INTERPOLATE_VALUES_AND_WRITE_SINGLE_FILE: {
      Stats::startEvent("interpolate values");
      std::string filenamePrefix;
      MASTER_EXCLUSIVE_SECTION {
        MPIUtils::receiveClass(&filenamePrefix, theMPISystem()->getManagerRank(),
                               theMPISystem()->getGlobalComm());
      }
      auto values =
          interpolateValues(receiveAndBroadcastInterpolationCoords(combiParameters_.getDim()));
      // write result
      MASTER_EXCLUSIVE_SECTION {
        std::string valuesWriteFilename =
            filenamePrefix + "_values_" + std::to_string(currentCombi_) + ".h5";
        writeInterpolatedValues(values, valuesWriteFilename);
      }
      Stats::stopEvent("interpolate values");
    } break;
    case WRITE_INTERPOLATED_VALUES_PER_GRID: {  // interpolate values on given coordinates and write
                                                // values to .h5
      Stats::startEvent("write interpolated values");
      writeInterpolatedValuesPerGrid(
          receiveAndBroadcastInterpolationCoords(combiParameters_.getDim()),
          receiveStringFromManagerAndBroadcastToGroup());
      Stats::stopEvent("write interpolated values");
    } break;
    case RESCHEDULE_ADD_TASK: {
      receiveAndInitializeTask();  // receive and initalize new task
                                   // now the newly
                                   // received task is the last in this->getTaskWorker().getTasks()
      // now , we may need to update the kahan summation data structures
      for (auto& dsg : combinedUniDSGVector_) {
        dsg->createKahanBuffer();
      }

      auto& currentTask = this->getTaskWorker().getTasks().back();
      currentTask->setZero();
      fillDFGFromDSGU(*currentTask);
      currentTask->setFinished(true);
    } break;
    case RESCHEDULE_REMOVE_TASK: {
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
      for (size_t i = 0; i < this->getTaskWorker().getTasks().size(); ++i) {
        if (this->getTaskWorker().getTasks()[i]->getID() == taskID) {
          MASTER_EXCLUSIVE_SECTION {
            // send to group master
            Task::send(this->getTaskWorker().getTasks()[i].get(), theMPISystem()->getManagerRank(),
                       theMPISystem()->getGlobalComm());
          }
          this->getTaskWorker().removeTask(i);
          break;  // only one task has the taskID
        }
      }
    } break;
    case WRITE_DSG_MINMAX_COEFFICIENTS: {
      writeSparseGridMinMaxCoefficients(receiveStringFromManagerAndBroadcastToGroup());
    } break;
    default: {
      throw std::runtime_error("signal " + std::to_string(signal) + " not implemented");
    }
  }

  // all tasks finished -> group waiting
  if (status_ != PROCESS_GROUP_FAIL) {
    status_ = PROCESS_GROUP_WAIT;
  }

  // send ready status to manager
  MASTER_EXCLUSIVE_SECTION {
    MPI_Send(&status_, 1, MPI_INT, theMPISystem()->getManagerRank(), TRANSFER_STATUS_TAG,
             theMPISystem()->getGlobalComm());
  }

  // if failed proc in this group detected the alive procs go into recovery state
  if (ENABLE_FT) {
    if (status_ == PROCESS_GROUP_FAIL) {
      theMPISystem()->recoverCommunicators(false);
      status_ = PROCESS_GROUP_WAIT;
    }
  }
  return signal;
}

void ProcessGroupWorker::runAllTasks() {
  Stats::startEvent("run");
  status_ = PROCESS_GROUP_BUSY;  // not sure if this does anything
  this->getTaskWorker().runAllTasks();
  Stats::stopEvent("run");
}

void ProcessGroupWorker::exit() {
  // write out tasks that were in use when the computation ended
  // (i.e. after fault tolerance or rescheduling changes)
  MASTER_EXCLUSIVE_SECTION {
    // serialize tasks as string
    std::stringstream tasksStream;
    for (const auto& t : this->getTaskWorker().getTasks()) {
      tasksStream << t->getID() << ": " << t->getCoefficient() << t->getLevelVector() << "; ";
    }
    std::string tasksString = tasksStream.str();
    // Stats::setAttribute("tasks: levels", tasksString);
  }
  if (isGENE) {
    if (chdir("../ginstance")) {
    };
  }
  this->getTaskWorker().deleteTasks();
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
  if (this->getTaskWorker().getTasks().size() == 0) {
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

  // get all subspaces in the (optimized) combischeme, create dsgs
  combinedUniDSGVector_.resize(static_cast<size_t>(combiParameters_.getNumGrids()));
  for (auto& uniDSG : combinedUniDSGVector_) {
    uniDSG = std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>(
        new DistributedSparseGridUniform<CombiDataType>(combiParameters_.getDim(), lmax, lmin,
                                                        theMPISystem()->getLocalComm()));
    // // this registers all possible subspaces in the DSGU
    // // can be used to test the memory consumption of the "filled" DSGU
    // registerAllSubspacesInDSGU(*uniDSG, combiParameters_);
  }

  // register dsgs in all dfgs
  Stats::startEvent("register dsgus");
  for (size_t g = 0; g < combinedUniDSGVector_.size(); ++g) {
    for (const auto& t : this->getTaskWorker().getTasks()) {
      DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid(static_cast<int>(g));
      // set subspace sizes locally
      combinedUniDSGVector_[g]->registerDistributedFullGrid(dfg);
    }
    // // we may clear the levels_ member of the sparse grids here to save memory
    // // but only if we need no new full grids initialized from the sparse grids!
    // // ...such as for rescheduling or interpolation (parallelEval/ evalNorm / ...)
    // combinedUniDSGVector_[(size_t) g]->resetLevels();

    // create the kahan buffer now, so it has only the subspaces present on the grids in this
    // process group
    combinedUniDSGVector_[g]->createKahanBuffer();
  }
  Stats::stopEvent("register dsgus");

  // global reduce of subspace sizes
  CommunicatorType globalReduceComm = theMPISystem()->getGlobalReduceComm();
  for (auto& uniDSG : combinedUniDSGVector_) {
    CombiCom::reduceSubspaceSizes(*uniDSG, globalReduceComm);
  }
}

void ProcessGroupWorker::addFullGridsToUniformSG() {
  assert(combinedUniDSGVector_.size() > 0 &&
         "Initialize dsgu first with "
         "initCombinedUniDSGVector()");
  auto numGrids = combiParameters_.getNumGrids();
  for (const auto& t : this->getTaskWorker().getTasks()) {
    for (IndexType g = 0; g < numGrids; g++) {
      DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid(static_cast<int>(g));

      // lokales reduce auf sg ->
      combinedUniDSGVector_[g]->addDistributedFullGrid(dfg, t->getCoefficient());
    }
  }
}

void ProcessGroupWorker::reduceUniformSG(RankType globalReduceRankThatCollects) {
  // we assume here that every task has the same number of grids, e.g. species in GENE
  auto numGrids = combiParameters_.getNumGrids();

  for (IndexType g = 0; g < numGrids; g++) {
    CombiCom::distributedGlobalReduce(*combinedUniDSGVector_[g], globalReduceRankThatCollects);
    assert(CombiCom::sumAndCheckSubspaceSizes(*combinedUniDSGVector_[g]));
  }
}

void ProcessGroupWorker::combineLocalAndGlobal(RankType globalReduceRankThatCollects) {
  assert(combinedUniDSGVector_.size() > 0 &&
         "Initialize dsgu first with "
         "initCombinedUniDSGVector()");
  // assert(combinedUniDSGVector_[0]->isSubspaceDataCreated());

  zeroDsgsData();

  Stats::startEvent("hierarchize");
  this->getTaskWorker().hierarchizeFullGrids(
      combiParameters_.getBoundary(), combiParameters_.getHierarchizationDims(),
      combiParameters_.getHierarchicalBases(), combiParameters_.getLMin());
  Stats::stopEvent("hierarchize");

  Stats::startEvent("local reduce");
  addFullGridsToUniformSG();
  Stats::stopEvent("local reduce");

  Stats::startEvent("global reduce");
  reduceUniformSG(globalReduceRankThatCollects);
  Stats::stopEvent("global reduce");
}

void ProcessGroupWorker::combineUniform() {
  Stats::startEvent("combine");
  combineLocalAndGlobal();
  integrateCombinedSolution();
  Stats::stopEvent("combine");
}

void ProcessGroupWorker::parallelEval() {
  if (uniformDecomposition)
    parallelEvalUniform(receiveStringFromManagerAndBroadcastToGroup(), receiveLevalAndBroadcast());
  else
    throw std::runtime_error("parallelEval not implemented for non-uniform decomposition");
}
// cf https://stackoverflow.com/questions/874134/find-out-if-string-ends-with-another-string-in-c
static bool endsWith(const std::string& str, const std::string& suffix) {
  return str.size() >= suffix.size() &&
         0 == str.compare(str.size() - suffix.size(), suffix.size(), suffix);
}

// helper function to output bool vector
inline std::ostream& operator<<(std::ostream& os, const std::vector<bool>& l) {
  os << "[";

  for (size_t i = 0; i < l.size(); ++i) os << l[i] << " ";

  os << "]";

  return os;
}

void ProcessGroupWorker::fillDFGFromDSGU(DistributedFullGrid<CombiDataType>& dfg,
                                         IndexType g) const {
  // fill dfg with hierarchical coefficients from distributed sparse grid
  dfg.extractFromUniformSG(*combinedUniDSGVector_[g]);

  bool anyNotBoundary =
      std::any_of(combiParameters_.getBoundary().begin(), combiParameters_.getBoundary().end(),
                  [](BoundaryType b) { return b == 0; });

  if (anyNotBoundary) {
    LevelVector zeroLMin = LevelVector(combiParameters_.getDim(), 0);
    DistributedHierarchization::dehierarchizeDFG(dfg, combiParameters_.getHierarchizationDims(),
                                                 combiParameters_.getHierarchicalBases(), zeroLMin);
  } else {
    DistributedHierarchization::dehierarchizeDFG(dfg, combiParameters_.getHierarchizationDims(),
                                                 combiParameters_.getHierarchicalBases(),
                                                 combiParameters_.getLMin());
  }
}

void ProcessGroupWorker::fillDFGFromDSGU(Task& t) const {
  auto numGrids = static_cast<int>(
      combiParameters_
          .getNumGrids());  // we assume here that every task has the same number of grids
  for (int g = 0; g < numGrids; g++) {
    assert(combinedUniDSGVector_[g] != nullptr);
    this->fillDFGFromDSGU(t.getDistributedFullGrid(g), g);
  }
}

void ProcessGroupWorker::parallelEvalUniform(std::string filename, LevelVector leval) {
  assert(uniformDecomposition);

  assert(combiParametersSet_);
  auto numGrids =
      combiParameters_
          .getNumGrids();  // we assume here that every task has the same number of grids

  const auto dim = static_cast<DimType>(leval.size());

  for (IndexType g = 0; g < numGrids; g++) {  // loop over all grids and plot them
    // create dfg
    bool forwardDecomposition = combiParameters_.getForwardDecomposition();
    auto levalDecomposition = combigrid::downsampleDecomposition(
        combiParameters_.getDecomposition(), combiParameters_.getLMax(), leval,
        combiParameters_.getBoundary());

    DistributedFullGrid<CombiDataType> dfg(
        dim, leval, theMPISystem()->getLocalComm(), combiParameters_.getBoundary(),
        combiParameters_.getParallelization(), forwardDecomposition, levalDecomposition);
    this->fillDFGFromDSGU(dfg, g);
    // save dfg to file with MPI-IO
    if (endsWith(filename, ".vtk")) {
      dfg.writePlotFileVTK(filename.c_str());
    } else {
      std::string fn = filename;
      auto pos = fn.find(".");
      if (pos != std::string::npos) {
        // if filename contains ".", insert grid number before that
        fn.insert(pos, "_" + std::to_string(g));
      }
      dfg.writePlotFile(fn.c_str());
    }
  }
}

std::vector<double> ProcessGroupWorker::getLpNorms(int p) const {
  return this->getTaskWorker().getLpNorms(p);
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
  for (const auto& task : this->getTaskWorker().getTasks()) {
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

std::vector<CombiDataType> ProcessGroupWorker::interpolateValues(
    const std::vector<std::vector<real>>& interpolationCoords) const {
  assert(combiParameters_.getNumGrids() == 1 && "interpolate only implemented for 1 species!");
  auto numCoordinates = interpolationCoords.size();

  // call interpolation function on tasks and reduce with combination coefficient
  std::vector<CombiDataType> values(numCoordinates, 0.);
  std::vector<CombiDataType> kahanTrailingTerm(numCoordinates, 0.);

  for (const auto& t : this->getTaskWorker().getTasks()) {
    const auto coeff = t->getCoefficient();
    for (size_t i = 0; i < numCoordinates; ++i) {
      auto localValue = t->getDistributedFullGrid().evalLocal(interpolationCoords[i]);
      auto summand = localValue * coeff;
      // cf. https://en.wikipedia.org/wiki/Kahan_summation_algorithm
      volatile auto y = summand - kahanTrailingTerm[i];
      volatile auto t = values[i] + y;
      kahanTrailingTerm[i] = (t - values[i]) - y;
      values[i] = t;
    }
  }
  // reduce interpolated values within process group
  MPI_Allreduce(MPI_IN_PLACE, values.data(), static_cast<int>(numCoordinates),
                abstraction::getMPIDatatype(abstraction::getabstractionDataType<CombiDataType>()),
                MPI_SUM, theMPISystem()->getLocalComm());
  //TODO is it necessary to correct for the kahan terms across process groups too?
  // need to reduce across process groups too
  // these do not strictly need to be allreduce (could be reduce), but it is easier to maintain that
  // way (all processes end up with valid values)
  MPI_Allreduce(MPI_IN_PLACE, values.data(), static_cast<int>(numCoordinates),
                abstraction::getMPIDatatype(abstraction::getabstractionDataType<CombiDataType>()),
                MPI_SUM, theMPISystem()->getGlobalReduceComm());

  // hope for RVO or change
  return values;
}

void ProcessGroupWorker::writeInterpolatedValuesPerGrid(
    const std::vector<std::vector<real>>& interpolationCoords, std::string fileNamePrefix) const {
  assert(combiParameters_.getNumGrids() == 1 && "interpolate only implemented for 1 species!");
  // call interpolation function on tasks and write out task-wise
  for (size_t i = 0; i < this->getTaskWorker().getTasks().size(); ++i) {
    auto taskVals =
        this->getTaskWorker().getTasks()[i]->getDistributedFullGrid().getInterpolatedValues(
            interpolationCoords);
    // cycle through ranks to write
    if (i % (theMPISystem()->getNumProcs()) == theMPISystem()->getLocalRank()) {
      std::string saveFilePath = fileNamePrefix + "_task_" +
                                 std::to_string(this->getTaskWorker().getTasks()[i]->getID()) +
                                 ".h5";
      std::string groupName = "run_";
      std::string datasetName = "interpolated_" + std::to_string(currentCombi_);
      h5io::writeValuesToH5File(taskVals, saveFilePath, groupName, datasetName,
                                this->getTaskWorker().getTasks()[i]->getCurrentTime());
    }
  }
}

void ProcessGroupWorker::writeInterpolatedValues(const std::vector<CombiDataType>& values,
                                                 const std::string& valuesWriteFilename) const {
  assert(combiParameters_.getNumGrids() == 1 && "interpolate only implemented for 1 species!");
  std::string groupName = "all_grids";
  std::string datasetName = "interpolated_" + std::to_string(currentCombi_);
  assert(this->getTaskWorker().getTasks().size() > 0);
  assert(currentCombi_ >= 0);
  assert(this->getTaskWorker().getTasks()[0]->getCurrentTime() >= 0.0);
  assert(values.size() > 0);
  assert(valuesWriteFilename.size() > 0);
  h5io::writeValuesToH5File(values, valuesWriteFilename, groupName, datasetName,
                            this->getTaskWorker().getTasks()[0]->getCurrentTime());
}

void ProcessGroupWorker::writeInterpolatedValuesSingleFile(
    const std::vector<std::vector<real>>& interpolationCoords,
    const std::string& filenamePrefix) const {
  // all processes interpolate
  assert(combiParameters_.getNumGrids() == 1 && "interpolate only implemented for 1 species!");
  auto values = interpolateValues(interpolationCoords);
  // one process writes
  OTHER_OUTPUT_GROUP_EXCLUSIVE_SECTION {
    MASTER_EXCLUSIVE_SECTION {
      writeInterpolatedValues(values,
                              filenamePrefix + +"_values_" + std::to_string(currentCombi_) + ".h5");
    }
  }
}

void ProcessGroupWorker::writeSparseGridMinMaxCoefficients(std::string fileNamePrefix) const {
  for (size_t i = 0; i < combinedUniDSGVector_.size(); ++i) {
    DistributedSparseGridIO::writeMinMaxCoefficents(*combinedUniDSGVector_[i], fileNamePrefix, i);
  }
}

void ProcessGroupWorker::receiveAndInitializeTask() {
  Task* t;

  // local root receives task
  MASTER_EXCLUSIVE_SECTION {
    Task::receive(&t, theMPISystem()->getManagerRank(), theMPISystem()->getGlobalComm());
  }

  // broadcast task to other process of pgroup
  Task::broadcast(&t, theMPISystem()->getMasterRank(), theMPISystem()->getLocalComm());

  this->initializeTask(std::unique_ptr<Task>(t));
}

void ProcessGroupWorker::initializeTask(std::unique_ptr<Task> t) {
  auto taskDecomposition = combigrid::downsampleDecomposition(
      combiParameters_.getDecomposition(), combiParameters_.getLMax(), t->getLevelVector(),
      combiParameters_.getBoundary());

  this->getTaskWorker().initializeTask(std::move(t), taskDecomposition,
                                       theMPISystem()->getLocalComm());
}

void ProcessGroupWorker::setCombiParameters(const CombiParameters& combiParameters) {
  combiParameters_ = combiParameters;
  combiParametersSet_ = true;

  // overwrite local comm with cartesian communicator
  if (!isGENE && combiParameters_.isParallelizationSet()) {
    // cf. https://www.rookiehpc.com/mpi/docs/mpi_cart_create.php
    // get decompositon from combi params
    auto par = combiParameters_.getParallelization();

    // important: note reverse ordering of dims! -- cf DistributedFullGrid
    //TODO(pollinta) remove reverse ordering
    std::vector<int> dims(par.size());
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

void ProcessGroupWorker::updateCombiParameters() {
  // local root receives combi parameters
  CombiParameters combiParametersReceived;
  MASTER_EXCLUSIVE_SECTION {
    MPIUtils::receiveClass(&combiParametersReceived, theMPISystem()->getManagerRank(), theMPISystem()->getGlobalComm());
  }
  // broadcast parameters to other processes of pgroup
  MPIUtils::broadcastClass(&combiParametersReceived, theMPISystem()->getMasterRank(), theMPISystem()->getLocalComm());

  this->setCombiParameters(combiParametersReceived);
}

void ProcessGroupWorker::integrateCombinedSolution() {
  auto numGrids = static_cast<int>(combiParameters_.getNumGrids());
  for (const auto& taskToUpdate : this->getTaskWorker().getTasks()) {
    for (int g = 0; g < numGrids; g++) {
      // fill dfg with hierarchical coefficients from distributed sparse grid
      taskToUpdate->getDistributedFullGrid(g).extractFromUniformSG(*combinedUniDSGVector_[g]);
    }
  }
  Stats::startEvent("dehierarchize");
  this->getTaskWorker().dehierarchizeFullGrids(
      combiParameters_.getBoundary(), combiParameters_.getHierarchizationDims(),
      combiParameters_.getHierarchicalBases(), combiParameters_.getLMin());
  Stats::stopEvent("dehierarchize");
  currentCombi_++;
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
    assert(dsgToUse->getRawDataSize() < 2147483647 &&
           "Dsg is larger than 2^31-1 and can not be "
           "sent in a single MPI Call (not "
           "supported yet) try a more coarse"
           "decomposition");
    // if we have an extra dsg for third level exchange, we use it
    if (extraUniDSGVector_.size() > 0) {
      dsgToUse->copyDataFrom(*uniDsg);
    }

    // send dsg data to manager
    Stats::startEvent("send dsg data");
    CombiCom::sendDsgData(*dsgToUse, manager, managerComm);
    Stats::stopEvent("send dsg data");

    // recv combined dsgu from manager
    Stats::startEvent("recv dsg data");
    CombiCom::recvDsgData(*dsgToUse, manager, managerComm);
    Stats::stopEvent("recv dsg data");

    if (extraUniDSGVector_.size() > 0) {
      // copy partial data from extraDSG back to uniDSG
      uniDsg->copyDataFrom(*dsgToUse);
    }

    // distribute solution in globalReduceComm to other pgs
    auto request = CombiCom::asyncBcastDsgData(*uniDsg, globalReduceRank, globalReduceComm);
    requests.push_back(request);
  }
  // update fgs
  integrateCombinedSolution();

  // wait for bcasts to other pgs in globalReduceComm
  Stats::startEvent("wait for bcasts");
  for (MPI_Request& request : requests) {
    auto returnedValue = MPI_Wait(&request, MPI_STATUS_IGNORE);
    assert(returnedValue == MPI_SUCCESS);
  }
  Stats::stopEvent("wait for bcasts");
}

void ProcessGroupWorker::combineThirdLevelFileBasedWrite(std::string filenamePrefixToWrite,
                                                         std::string writeCompleteTokenFileName) {
  assert(combinedUniDSGVector_.size() != 0);
  assert(combiParametersSet_);

  // write sparse grid and corresponding token file
  Stats::startEvent("write SG");
  this->writeDSGsToDisk(filenamePrefixToWrite);
  MASTER_EXCLUSIVE_SECTION { std::ofstream tokenFile(writeCompleteTokenFileName); }
  Stats::stopEvent("write SG");
}

void ProcessGroupWorker::combineThirdLevelFileBasedReadReduce(std::string filenamePrefixToRead,
                                                              std::string startReadingTokenFileName,
                                                              bool overwrite,
                                                              bool keepSparseGridFiles) {
  // wait until we can start to read
  Stats::startEvent("wait SG");
  MASTER_EXCLUSIVE_SECTION {
    std::cout << "Waiting for token file " << startReadingTokenFileName << std::endl;
    while (!std::filesystem::exists(startReadingTokenFileName)) {
      // wait for token file to appear
      std::this_thread::sleep_for(std::chrono::seconds(1));
    }
  }
  MPI_Barrier(theMPISystem()->getOutputGroupComm());
  Stats::stopEvent("wait SG");

  if (overwrite) {
    Stats::startEvent("read SG");
    this->readDSGsFromDisk(filenamePrefixToRead);
    Stats::stopEvent("read SG");
  } else {
    Stats::startEvent("read/reduce SG");
    this->readDSGsFromDiskAndReduce(filenamePrefixToRead);
    Stats::stopEvent("read/reduce SG");
  }
  if (combinedUniDSGVector_.size() != 1) {
    throw std::runtime_error("Combining more than one DSG is not implemented yet");
  }
  // distribute solution in globalReduceComm to other pgs
  auto request =
      CombiCom::asyncBcastDsgData(*combinedUniDSGVector_[0], theMPISystem()->getGlobalReduceRank(),
                        theMPISystem()->getGlobalReduceComm());

  // update fgs
  integrateCombinedSolution();

  // remove reading token
  MASTER_EXCLUSIVE_SECTION {
    std::filesystem::remove(startReadingTokenFileName);
    // remove sparse grid file
    if (!keepSparseGridFiles) {
      std::filesystem::remove(filenamePrefixToRead + "_" + std::to_string(0));
    }
  }

  auto returnedValue = MPI_Wait(&request, MPI_STATUS_IGNORE);
  assert(returnedValue == MPI_SUCCESS);
}

void ProcessGroupWorker::combineThirdLevelFileBased(std::string filenamePrefixToWrite,
                                                    std::string writeCompleteTokenFileName,
                                                    std::string filenamePrefixToRead,
                                                    std::string startReadingTokenFileName) {
  this->combineThirdLevelFileBasedWrite(filenamePrefixToWrite, writeCompleteTokenFileName);
  this->combineThirdLevelFileBasedReadReduce(filenamePrefixToRead, startReadingTokenFileName);
}

void ProcessGroupWorker::setExtraSparseGrid(bool initializeSizes) {
  if (combinedUniDSGVector_.size() != 1) {
    throw std::runtime_error("combinedUniDSGVector_ is empty");
  }

  if (extraUniDSGVector_.empty()) {
    // create new vector for extra sparse grids (that will be only on this process group)
    extraUniDSGVector_.resize(combinedUniDSGVector_.size());
    for (auto& extraUniDSG : extraUniDSGVector_) {
      extraUniDSG = std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>(
          new DistributedSparseGridUniform<CombiDataType>(
              combinedUniDSGVector_[0]->getDim(), combinedUniDSGVector_[0]->getAllLevelVectors(),
              theMPISystem()->getOutputGroupComm()));
      // create Kahan buffer now (at zero size), because summation is not needed on this sparse grid
      extraUniDSG->createKahanBuffer();
      if (initializeSizes) {
        for (size_t i = 0; i < extraUniDSG->getNumSubspaces(); ++i) {
          extraUniDSG->setDataSize(i, combinedUniDSGVector_[0]->getDataSize(i));
        }
      }
      // level vectors are not required; read from the initial sparse grid if needed
      extraUniDSG->resetLevels();
    }
  } else {
    throw std::runtime_error(
        "extraUniDSGVector_ is not empty -- if you think this is ok, try to remove the if-else "
        "here");
  }
}

/** Reduces subspace sizes with remote.
 */
void ProcessGroupWorker::reduceSubspaceSizesThirdLevel(bool thirdLevelExtraSparseGrid) {
  assert(combiParametersSet_);
  // update either old or new sparse grids
  auto* uniDSGVectorToSet = &combinedUniDSGVector_;
  if (thirdLevelExtraSparseGrid) {
    this->setExtraSparseGrid(false);  // don't initialize, would be overwritten
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
    CombiCom::sendSubspaceSizesWithGather(*combinedUniDSGVector_[i], thirdLevelComm, thirdLevelManagerRank);
    // set updated sizes in dsgs
    CombiCom::receiveSubspaceSizesWithScatter(*(*uniDSGVectorToSet)[i], thirdLevelComm, thirdLevelManagerRank);
  }

  if (!thirdLevelExtraSparseGrid) {
    // distribute updated sizes to workers with same decomposition (global reduce comm)
    // cf. waitForThirdLevelSizeUpdate(), which is called in other process groups
    CommunicatorType globalReduceComm = theMPISystem()->getGlobalReduceComm();
    RankType globalReduceRank = theMPISystem()->getGlobalReduceRank();
    for (auto& dsg : combinedUniDSGVector_) {
      CombiCom::broadcastSubspaceSizes(*dsg, globalReduceComm, globalReduceRank);
    }
  }
}

void ProcessGroupWorker::waitForThirdLevelSizeUpdate() {
  RankType thirdLevelPG = (RankType)combiParameters_.getThirdLevelPG();
  CommunicatorType globalReduceComm = theMPISystem()->getGlobalReduceComm();

  for (auto& dsg : combinedUniDSGVector_) {
    CombiCom::broadcastSubspaceSizes(*dsg, globalReduceComm, thirdLevelPG);
  }
}

void ProcessGroupWorker::reduceSubspaceSizes(const std::string& filenameToRead,
                                             bool extraSparseGrid, bool overwrite) {
  if (extraSparseGrid) {
    OUTPUT_GROUP_EXCLUSIVE_SECTION {
      this->setExtraSparseGrid(true);
#ifndef NDEBUG
      // duplicate subspace sizes to validate later
      std::vector<SubspaceSizeType> subspaceSizesToValidate =
          extraUniDSGVector_[0]->getSubspaceDataSizes();
#endif
      // use extra sparse grid
      if (overwrite) {
        DistributedSparseGridIO::readSubspaceSizesFromFile(*extraUniDSGVector_[0], filenameToRead,
                                                           false);
      } else {
        auto minFunctionInstantiation = [](SubspaceSizeType a, SubspaceSizeType b) {
          return std::min(a, b);
        };
        DistributedSparseGridIO::readReduceSubspaceSizesFromFile(
            *extraUniDSGVector_[0], filenameToRead, minFunctionInstantiation, 0, false);
      }
#ifndef NDEBUG
      assert(subspaceSizesToValidate.size() ==
             extraUniDSGVector_[0]->getSubspaceDataSizes().size());
      for (size_t i = 0; i < subspaceSizesToValidate.size(); ++i) {
        assert(extraUniDSGVector_[0]->getSubspaceDataSizes()[i] == 0 ||
               extraUniDSGVector_[0]->getSubspaceDataSizes()[i] == subspaceSizesToValidate[i]);
      }
      auto numDOFtoValidate =
          std::accumulate(subspaceSizesToValidate.begin(), subspaceSizesToValidate.end(), 0);
      auto numDOFnow = std::accumulate(extraUniDSGVector_[0]->getSubspaceDataSizes().begin(),
                                       extraUniDSGVector_[0]->getSubspaceDataSizes().end(), 0);
      assert(numDOFtoValidate >= numDOFnow);
#endif
    }
  } else {
    if (!extraUniDSGVector_.empty()) {
      throw std::runtime_error("extraUniDSGVector_ not empty, but extraSparseGrid is false");
    }
#ifndef NDEBUG
    std::vector<SubspaceSizeType> subspaceSizesToValidate =
        combinedUniDSGVector_[0]->getSubspaceDataSizes();
#endif
    FIRST_GROUP_EXCLUSIVE_SECTION {
      if (overwrite) {
        DistributedSparseGridIO::readSubspaceSizesFromFile(*combinedUniDSGVector_[0],
                                                           filenameToRead, true);
      } else {
        // if no extra sparse grid, max-reduce the normal one
        auto maxFunctionInstantiation = [](SubspaceSizeType a, SubspaceSizeType b) {
          return std::max(a, b);
        };
        DistributedSparseGridIO::readReduceSubspaceSizesFromFile(
            *combinedUniDSGVector_[0], filenameToRead, maxFunctionInstantiation, 0, true);
      }
      if (theMPISystem()->getGlobalReduceRank() != 0) {
        throw std::runtime_error("read rank is not the global reduce rank");
      }
    }
    else {
      if (theMPISystem()->getGlobalReduceRank() == 0) {
        throw std::runtime_error("read rank IS the global reduce rank");
      }
    }
    // reduce to all other process groups
    CommunicatorType globalReduceComm = theMPISystem()->getGlobalReduceComm();
    RankType senderRank = 0;
    CombiCom::broadcastSubspaceSizes(*combinedUniDSGVector_[0], globalReduceComm, senderRank);
#ifndef NDEBUG
    assert(subspaceSizesToValidate.size() ==
           combinedUniDSGVector_[0]->getSubspaceDataSizes().size());
    for (size_t i = 0; i < subspaceSizesToValidate.size(); ++i) {
      assert(subspaceSizesToValidate[i] == 0 ||
             subspaceSizesToValidate[i] == combinedUniDSGVector_[0]->getSubspaceDataSizes()[i]);
    }
    auto numDOFtoValidate =
        std::accumulate(subspaceSizesToValidate.begin(), subspaceSizesToValidate.end(), 0);
    auto numDOFnow = std::accumulate(combinedUniDSGVector_[0]->getSubspaceDataSizes().begin(),
                                     combinedUniDSGVector_[0]->getSubspaceDataSizes().end(), 0);
    assert(numDOFtoValidate <= numDOFnow);
#endif
  }
}

void ProcessGroupWorker::reduceSubspaceSizesFileBased(std::string filenamePrefixToWrite,
                                                      std::string writeCompleteTokenFileName,
                                                      std::string filenamePrefixToRead,
                                                      std::string startReadingTokenFileName,
                                                      bool extraSparseGrid) {
  assert(combinedUniDSGVector_.size() == 1);
  FIRST_GROUP_EXCLUSIVE_SECTION {
    DistributedSparseGridIO::writeSubspaceSizesToFile(*combinedUniDSGVector_[0],
                                                      filenamePrefixToWrite);
    MASTER_EXCLUSIVE_SECTION { std::ofstream tokenFile(writeCompleteTokenFileName); }
  }

  // if extra sparse grid, we only need to do something if we are the I/O group
  OUTPUT_GROUP_EXCLUSIVE_SECTION {}
  else if (extraSparseGrid) {
    // otherwise, return
    return;
  }

  // wait until we can start to read
  MASTER_EXCLUSIVE_SECTION {
    while (!std::filesystem::exists(startReadingTokenFileName)) {
      // wait for token file to appear
      std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }
  }
  if (extraSparseGrid) {
    MPI_Barrier(theMPISystem()->getOutputGroupComm());
  } else {
    MPI_Barrier(theMPISystem()->getLocalComm());
  }
  this->reduceSubspaceSizes(filenamePrefixToRead, extraSparseGrid);
}

void ProcessGroupWorker::waitForThirdLevelCombiResult(bool fromOutputGroup) {
  assert(extraUniDSGVector_.empty());
  RankType broadcastSender;
  if (fromOutputGroup) {
    broadcastSender = theMPISystem()->getOutputRankInGlobalReduceComm();
  } else {
  // receive third level combi result from third level pgroup (global reduce comm)
    broadcastSender = (RankType)combiParameters_.getThirdLevelPG();
  }
  CommunicatorType globalReduceComm = theMPISystem()->getGlobalReduceComm();

  Stats::startEvent("wait for bcasts");
  for (auto& dsg : combinedUniDSGVector_) {
    auto request = CombiCom::asyncBcastDsgData(*dsg, broadcastSender, globalReduceComm);
    auto returnedValue = MPI_Wait(&request, MPI_STATUS_IGNORE);
    assert(returnedValue == MPI_SUCCESS);
  }
  Stats::stopEvent("wait for bcasts");

  integrateCombinedSolution();
}

void ProcessGroupWorker::zeroDsgsData() {
  for (auto& dsg : combinedUniDSGVector_)
    dsg->setZero();
  for (auto& dsg : extraUniDSGVector_)
    dsg->setZero();
}

/** free dsgus space */
void ProcessGroupWorker::deleteDsgsData() {
  for (auto& dsg : combinedUniDSGVector_)
    dsg->deleteSubspaceData();
  for (auto& dsg : extraUniDSGVector_)
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
  for (const auto& task : this->getTaskWorker().getTasks()) {
    writeVTKPlotFileOfTask(*task);
  }
}

void ProcessGroupWorker::writeDSGsToDisk(std::string filenamePrefix) {
  for (size_t i = 0; i < combinedUniDSGVector_.size(); ++i) {
    auto filename = filenamePrefix + "_" + std::to_string(i);
    auto uniDsg = combinedUniDSGVector_[i].get();
    auto dsgToUse = uniDsg;
    if (extraUniDSGVector_.size() > 0) {
      dsgToUse = extraUniDSGVector_[i].get();
      dsgToUse->copyDataFrom(*uniDsg);
    }
    DistributedSparseGridIO::writeOneFile(*dsgToUse, filename);
  }
}

void ProcessGroupWorker::readDSGsFromDisk(std::string filenamePrefix, bool alwaysReadFullDSG) {
  for (size_t i = 0; i < combinedUniDSGVector_.size(); ++i) {
    auto uniDsg = combinedUniDSGVector_[i].get();
    auto dsgToUse = uniDsg;
    if (extraUniDSGVector_.size() > 0 && !alwaysReadFullDSG) {
      dsgToUse = extraUniDSGVector_[i].get();
    }
    DistributedSparseGridIO::readOneFile(*dsgToUse, filenamePrefix + "_" + std::to_string(i));
    if (extraUniDSGVector_.size() > 0) {
      // copy partial data from extraDSG back to uniDSG
      uniDsg->copyDataFrom(*dsgToUse);
    }
  }
}

void ProcessGroupWorker::readDSGsFromDiskAndReduce(std::string filenamePrefixToRead,
                                                   bool alwaysReadFullDSG) {
  for (size_t i = 0; i < combinedUniDSGVector_.size(); ++i) {
    auto uniDsg = combinedUniDSGVector_[i].get();
    auto dsgToUse = uniDsg;
    if (extraUniDSGVector_.size() > 0 && !alwaysReadFullDSG) {
      dsgToUse = extraUniDSGVector_[i].get();
    }
    // assume that at least for four process groups, we should have enough spare RAM
    // to read all of the sparse grid at once
    // if fewer, chunk the read/reduce
    int numberReduceChunks = 1;
    if (theMPISystem()->getNumGroups() == 1) {
      numberReduceChunks = 4;
    } else if (theMPISystem()->getNumGroups() < 4) {
      numberReduceChunks = 2;
    }
    DistributedSparseGridIO::readOneFileAndReduce(
        *dsgToUse, filenamePrefixToRead + "_" + std::to_string(i), numberReduceChunks);
    if (extraUniDSGVector_.size() > 0) {
      // copy partial data from extraDSG back to uniDSG
      uniDsg->copyDataFrom(*dsgToUse);
    }
  }
}

} /* namespace combigrid */
