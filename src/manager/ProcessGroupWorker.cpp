#include "manager/ProcessGroupWorker.hpp"

#include "combicom/CombiCom.hpp"
#include "manager/InterpolationWorker.hpp"
#include "manager/ProcessGroupSignals.hpp"
#include "mpi/MPIUtils.hpp"
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
  for (int p = 0; p < lpnorms.size(); ++p) {
    // send from master to manager
    MASTER_EXCLUSIVE_SECTION {
      MPI_Send(&lpnorms[p], 1, MPI_DOUBLE, theMPISystem()->getManagerRank(), TRANSFER_NORM_TAG,
               theMPISystem()->getGlobalComm());
    }
  }
}

SignalType receiveSignalFromManagerAndBroadcast() {
  SignalType signal = -1;

  MASTER_EXCLUSIVE_SECTION {
    // receive signal from manager
    MPI_Recv(&signal, 1, MPI_INT, theMPISystem()->getManagerRank(), TRANSFER_SIGNAL_TAG,
             theMPISystem()->getGlobalComm(), MPI_STATUS_IGNORE);
  }
  // distribute signal to other processes of pgroup
  MPI_Bcast(&signal, 1, MPI_INT, theMPISystem()->getMasterRank(), theMPISystem()->getLocalComm());
  return signal;
}

LevelVector receiveLevalAndBroadcast(DimType dim) {
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
      combiParameters_(),
      combiParametersSet_(false),
      currentCombi_(0) {}

ProcessGroupWorker::~ProcessGroupWorker() {}

SignalType ProcessGroupWorker::wait() {
  if (status_ == PROCESS_GROUP_FAIL) {  // in this case worker got reused
    status_ = PROCESS_GROUP_WAIT;
  }

  SignalType signal = receiveSignalFromManagerAndBroadcast();
  // process signal
  switch (signal) {
    case RUN_FIRST: {
      receiveAndInitializeTask();
      status_ = PROCESS_GROUP_BUSY;

      // execute task
      Stats::startEvent("run first");
      auto& currentTask = this->getTaskWorker().getLastTask();
      currentTask->run(theMPISystem()->getLocalComm());
      Stats::Event e = Stats::stopEvent("run first");
    } break;
    case RUN_NEXT: {
      assert(this->getTaskWorker().getTasks().size() > 0);
      // // free space for computation
      // this->getSparseGridWorker().deleteDsgsData();

      this->runAllTasks();

    } break;
    case ADD_TASK: {  // add a new task to the process group
      // initalize task and set values to zero
      // the task will get the proper initial solution during the next combine
      receiveAndInitializeTask();

      auto& currentTask = this->getTaskWorker().getLastTask();
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
      this->initCombinedDSGVector();
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
      combigrid::writeVTKPlotFilesOfAllTasks(this->getTaskWorker().getTasks(),
                                             combiParameters_.getNumGrids());
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
      auto& currentTask = this->getTaskWorker().getLastTask();

      currentTask->setZero();
      // fill task with combisolution
      this->getSparseGridWorker().fillDFGFromDSGU(
          *currentTask, combiParameters_.getHierarchizationDims(),
          combiParameters_.getHierarchicalBases(), combiParameters_.getLMin());

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
      writeInterpolatedValuesSingleFile(
          receiveAndBroadcastInterpolationCoords(combiParameters_.getDim()),
          receiveStringFromManagerAndBroadcastToGroup());
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
      for (auto& dsg : this->getSparseGridWorker().getCombinedUniDSGVector()) {
        dsg->createKahanBuffer();
      }

      auto& currentTask = this->getTaskWorker().getLastTask();
      currentTask->setZero();
      this->getSparseGridWorker().fillDFGFromDSGU(
          *currentTask, combiParameters_.getHierarchizationDims(),
          combiParameters_.getHierarchicalBases(), combiParameters_.getLMin());
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
      this->getSparseGridWorker().writeMinMaxCoefficients(
          receiveStringFromManagerAndBroadcastToGroup());
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

void ProcessGroupWorker::initCombinedDSGVector() {
  assert(combiParametersSet_);
  Stats::startEvent("init dsgus");
  this->getSparseGridWorker().initCombinedUniDSGVector(
      combiParameters_.getLMin(), combiParameters_.getLMax(),
      combiParameters_.getLMaxReductionVector(), combiParameters_.getNumGrids(),
      combiParameters_.getCombinationVariant());
  Stats::stopEvent("init dsgus");
}

void ProcessGroupWorker::combineLocalAndGlobal(RankType globalReduceRankThatCollects) {
  assert(this->getSparseGridWorker().getNumberOfGrids() > 0 &&
         "Initialize dsgu first with "
         "initCombinedUniDSGVector()");
  // assert(this->getSparseGridWorker().getCombinedUniDSGVector()[0]->isSubspaceDataCreated());

  this->getSparseGridWorker().zeroDsgsData();

  Stats::startEvent("hierarchize");
  this->getTaskWorker().hierarchizeFullGrids(
      combiParameters_.getBoundary(), combiParameters_.getHierarchizationDims(),
      combiParameters_.getHierarchicalBases(), combiParameters_.getLMin());
  Stats::stopEvent("hierarchize");

  Stats::startEvent("local reduce");
  this->getSparseGridWorker().addFullGridsToUniformSG();
  Stats::stopEvent("local reduce");

  Stats::startEvent("global reduce");
  this->getSparseGridWorker().reduceUniformSG(combiParameters_.getCombinationVariant(),
                                              globalReduceRankThatCollects);
  Stats::stopEvent("global reduce");
}

void ProcessGroupWorker::combineUniform() {
  Stats::startEvent("combine");
  combineLocalAndGlobal();
  updateFullFromCombinedSparseGrids();
  Stats::stopEvent("combine");
}

void ProcessGroupWorker::parallelEval() {
  if (uniformDecomposition) {
    parallelEvalUniform(receiveStringFromManagerAndBroadcastToGroup(),
                        receiveLevalAndBroadcast(combiParameters_.getDim()));
  } else {
    throw std::runtime_error("parallelEval not implemented for non-uniform decomposition");
  }
}

void ProcessGroupWorker::parallelEvalUniform(const std::string& filename,
                                             const LevelVector& leval) const {
  assert(uniformDecomposition);
  bool forwardDecomposition = combiParameters_.getForwardDecomposition();
  auto levalDecomposition = combigrid::downsampleDecomposition(combiParameters_.getDecomposition(),
                                                               combiParameters_.getLMax(), leval,
                                                               combiParameters_.getBoundary());
  this->getSparseGridWorker().interpolateAndPlotOnLevel(
      filename, leval, combiParameters_.getBoundary(), combiParameters_.getHierarchizationDims(),
      combiParameters_.getHierarchicalBases(), combiParameters_.getLMin(),
      combiParameters_.getParallelization(), levalDecomposition);
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
      task->doDiagnostics();
      return;
    }
  }
  assert(false && "this taskID is not here");
}

std::vector<CombiDataType> ProcessGroupWorker::interpolateValues(
    const std::vector<std::vector<real>>& interpolationCoords) const {
  assert(combiParameters_.getNumGrids() == 1 && "interpolate only implemented for 1 species!");
  return combigrid::interpolateValues<CombiDataType>(this->getTaskWorker().getTasks(),
                                                     interpolationCoords);
}

void ProcessGroupWorker::writeInterpolatedValuesPerGrid(
    const std::vector<std::vector<real>>& interpolationCoords,
    const std::string& fileNamePrefix) const {
  assert(combiParameters_.getNumGrids() == 1 && "interpolate only implemented for 1 species!");
  combigrid::writeInterpolatedValuesPerGrid(this->getTaskWorker().getTasks(), interpolationCoords,
                                            fileNamePrefix, currentCombi_);
}

void ProcessGroupWorker::writeInterpolatedValuesSingleFile(
    const std::vector<std::vector<real>>& interpolationCoords,
    const std::string& filenamePrefix) const {
  // all processes interpolate
  assert(combiParameters_.getNumGrids() == 1 && "interpolate only implemented for 1 species!");
  combigrid::writeInterpolatedValuesSingleFile<CombiDataType>(
      this->getTaskWorker().getTasks(), interpolationCoords, filenamePrefix, currentCombi_);
}

void ProcessGroupWorker::writeSparseGridMinMaxCoefficients(
    const std::string& fileNamePrefix) const {
  this->getSparseGridWorker().writeMinMaxCoefficients(fileNamePrefix);
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

void ProcessGroupWorker::setCombiParameters(CombiParameters&& combiParameters) {
  combiParameters_ = std::move(combiParameters);
  combiParametersSet_ = true;

  // overwrite local comm with cartesian communicator
  if (!isGENE && combiParameters_.isParallelizationSet()) {
    // cf. https://www.rookiehpc.com/mpi/docs/mpi_cart_create.php
    // get decompositon from combi params
    auto& dims = combiParameters_.getParallelization();

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

  this->setCombiParameters(std::move(combiParametersReceived));
}

void ProcessGroupWorker::updateFullFromCombinedSparseGrids() {
  Stats::startEvent("dehierarchize");
  this->getSparseGridWorker().integrateCombinedSolutionToTasks(
      combiParameters_.getBoundary(), combiParameters_.getHierarchizationDims(),
      combiParameters_.getHierarchicalBases(), combiParameters_.getLMin());
  Stats::stopEvent("dehierarchize");
  currentCombi_++;
}

void ProcessGroupWorker::combineThirdLevel() {
  assert(this->getSparseGridWorker().getNumberOfGrids() != 0);
  assert(combiParametersSet_);

  assert(theMPISystem()->getThirdLevelComms().size() == 1 && "init thirdLevel communicator failed");
  const CommunicatorType& managerComm = theMPISystem()->getThirdLevelComms()[0];
  const CommunicatorType& globalReduceComm = theMPISystem()->getGlobalReduceComm();
  const RankType& globalReduceRank = theMPISystem()->getGlobalReduceRank();
  const RankType& manager = theMPISystem()->getThirdLevelManagerRank();

  std::vector<MPI_Request> requests;
  for (size_t i = 0; i < this->getSparseGridWorker().getNumberOfGrids(); ++i) {
    auto uniDsg = this->getSparseGridWorker().getCombinedUniDSGVector()[i].get();
    auto dsgToUse = uniDsg;
    if (this->getSparseGridWorker().getExtraUniDSGVector().size() > 0) {
      dsgToUse = this->getSparseGridWorker().getExtraUniDSGVector()[i].get();
    }
    assert(dsgToUse->getRawDataSize() < 2147483647 &&
           "Dsg is larger than 2^31-1 and can not be "
           "sent in a single MPI Call (not "
           "supported yet) try a more coarse"
           "decomposition");
    // if we have an extra dsg for third level exchange, we use it
    if (this->getSparseGridWorker().getExtraUniDSGVector().size() > 0) {
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

    if (this->getSparseGridWorker().getExtraUniDSGVector().size() > 0) {
      // copy partial data from extraDSG back to uniDSG
      uniDsg->copyDataFrom(*dsgToUse);
    }

    // distribute solution in globalReduceComm to other pgs
    auto request =
        this->getSparseGridWorker().startBroadcastDSGs(combiParameters_.getCombinationVariant());
    requests.push_back(request);
  }
  // update fgs
  updateFullFromCombinedSparseGrids();

  // wait for bcasts to other pgs in globalReduceComm
  Stats::startEvent("wait for bcasts");
  for (MPI_Request& request : requests) {
    auto returnedValue = MPI_Wait(&request, MPI_STATUS_IGNORE);
    assert(returnedValue == MPI_SUCCESS);
  }
  Stats::stopEvent("wait for bcasts");
}

void ProcessGroupWorker::combineThirdLevelFileBasedWrite(
    const std::string& filenamePrefixToWrite, const std::string& writeCompleteTokenFileName) {
  assert(this->getSparseGridWorker().getNumberOfGrids() > 0);
  assert(combiParametersSet_);

  // write sparse grid and corresponding token file
  Stats::startEvent("write SG");
  this->writeDSGsToDisk(filenamePrefixToWrite);
  MASTER_EXCLUSIVE_SECTION { std::ofstream tokenFile(writeCompleteTokenFileName); }
  Stats::stopEvent("write SG");
}

void ProcessGroupWorker::combineThirdLevelFileBasedReadReduce(
    const std::string& filenamePrefixToRead, const std::string& startReadingTokenFileName,
    bool overwrite, bool keepSparseGridFiles) {
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

  overwrite ? Stats::startEvent("read SG") : Stats::startEvent("read/reduce SG");
  auto request = this->getSparseGridWorker().readReduceStartBroadcastDSGs(
      filenamePrefixToRead, this->combiParameters_.getCombinationVariant(), overwrite);
  overwrite ? Stats::stopEvent("read SG") : Stats::stopEvent("read/reduce SG");

  // update fgs
  updateFullFromCombinedSparseGrids();

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

void ProcessGroupWorker::combineThirdLevelFileBased(const std::string& filenamePrefixToWrite,
                                                    const std::string& writeCompleteTokenFileName,
                                                    const std::string& filenamePrefixToRead,
                                                    const std::string& startReadingTokenFileName) {
  this->combineThirdLevelFileBasedWrite(filenamePrefixToWrite, writeCompleteTokenFileName);
  this->combineThirdLevelFileBasedReadReduce(filenamePrefixToRead, startReadingTokenFileName);
}

void ProcessGroupWorker::setExtraSparseGrid(bool initializeSizes) {
  return this->getSparseGridWorker().setExtraSparseGrid(initializeSizes);
}

/** Reduces subspace sizes with remote.
 */
void ProcessGroupWorker::reduceSubspaceSizesThirdLevel(bool thirdLevelExtraSparseGrid) {
  assert(combiParametersSet_);
  // update either old or new sparse grid
  auto& uniDSGToSet =
      this->getSparseGridWorker().getSparseGridToUseForThirdLevel(thirdLevelExtraSparseGrid);

  // prepare for MPI calls to manager
  CommunicatorType thirdLevelComm = theMPISystem()->getThirdLevelComms()[0];
  RankType thirdLevelManagerRank = theMPISystem()->getThirdLevelManagerRank();
  CombiCom::sendSubspaceSizesWithGather(*this->getSparseGridWorker().getCombinedUniDSGVector()[0],
                                        thirdLevelComm, thirdLevelManagerRank);
  // set updated sizes in dsgs
  CombiCom::receiveSubspaceSizesWithScatter(*uniDSGToSet, thirdLevelComm, thirdLevelManagerRank);

  if (!thirdLevelExtraSparseGrid) {
    // distribute updated sizes to workers with same decomposition (global reduce comm)
    // cf. waitForThirdLevelSizeUpdate(), which is called in other process groups
    CommunicatorType globalReduceComm = theMPISystem()->getGlobalReduceComm();
    RankType globalReduceRank = theMPISystem()->getGlobalReduceRank();
    for (auto& dsg : this->getSparseGridWorker().getCombinedUniDSGVector()) {
      CombiCom::broadcastSubspaceSizes(*dsg, globalReduceComm, globalReduceRank);
    }
  }
}

void ProcessGroupWorker::waitForThirdLevelSizeUpdate() {
  RankType thirdLevelPG = (RankType)combiParameters_.getThirdLevelPG();
  CommunicatorType globalReduceComm = theMPISystem()->getGlobalReduceComm();

  for (auto& dsg : this->getSparseGridWorker().getCombinedUniDSGVector()) {
    CombiCom::broadcastSubspaceSizes(*dsg, globalReduceComm, thirdLevelPG);
  }
}

void ProcessGroupWorker::reduceSubspaceSizes(const std::string& filenameToRead,
                                             bool extraSparseGrid, bool overwrite) {
  this->getSparseGridWorker().reduceSubspaceSizes(
      filenameToRead, this->combiParameters_.getCombinationVariant(), extraSparseGrid, overwrite);
}

void ProcessGroupWorker::reduceSubspaceSizesFileBased(const std::string& filenamePrefixToWrite,
                                                      const std::string& writeCompleteTokenFileName,
                                                      const std::string& filenamePrefixToRead,
                                                      const std::string& startReadingTokenFileName,
                                                      bool extraSparseGrid) {
  FIRST_GROUP_EXCLUSIVE_SECTION {
    this->getSparseGridWorker().writeSubspaceSizesToFile(filenamePrefixToWrite);
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
  assert(this->getSparseGridWorker().getExtraUniDSGVector().empty());
  RankType broadcastSender;
  if (fromOutputGroup) {
    broadcastSender = theMPISystem()->getOutputRankInGlobalReduceComm();
  } else {
    // receive third level combi result from third level pgroup (global reduce comm)
    broadcastSender = (RankType)combiParameters_.getThirdLevelPG();
  }
  CommunicatorType globalReduceComm = theMPISystem()->getGlobalReduceComm();

  Stats::startEvent("wait for bcasts");
  auto request =
      this->getSparseGridWorker().startBroadcastDSGs(combiParameters_.getCombinationVariant());
  auto returnedValue = MPI_Wait(&request, MPI_STATUS_IGNORE);
  assert(returnedValue == MPI_SUCCESS);
  Stats::stopEvent("wait for bcasts");

  updateFullFromCombinedSparseGrids();
}

void ProcessGroupWorker::zeroDsgsData() { this->getSparseGridWorker().zeroDsgsData(); }

void ProcessGroupWorker::writeDSGsToDisk(const std::string& filenamePrefix) {
  this->getSparseGridWorker().writeDSGsToDisk(filenamePrefix);
}

void ProcessGroupWorker::readDSGsFromDisk(const std::string& filenamePrefix,
                                          bool alwaysReadFullDSG) {
  this->getSparseGridWorker().readDSGsFromDisk(filenamePrefix, alwaysReadFullDSG);
}

} /* namespace combigrid */
