#include "discotec/manager/ProcessManager.hpp"
#include <algorithm>
#include <iostream>

#include <boost/asio.hpp>

#include "combicom/CombiCom.hpp"
#include "discotec/io/H5InputOutput.hpp"
#include "discotec/utils/Types.hpp"
#include "discotec/mpi/MPIUtils.hpp"

namespace combigrid {

void ProcessManager::sortTasks(){
  LoadModel* lm = loadModel_.get();
  assert(lm);
  std::sort(tasks_.begin(), tasks_.end(),
            [lm](const Task* instance1, const Task* instance2){
                assert(instance1);
                return (lm->eval(instance1->getLevelVector()) > lm->eval(instance2->getLevelVector()));
            }
  );
}

bool ProcessManager::runfirst(bool doInitDSGUs) {
  // sort instances in decreasing order
  sortTasks();

  assert(tasks_.size() >= pgroups_.size());

  for (size_t i = 0; i < tasks_.size(); ++i) {
    // wait for available process group
    ProcessGroupManagerID g = wait();

    // assign instance to group
    g->runfirst(tasks_[i]);
  }

  bool group_failed = waitAllFinished();
  //size_t numDurationsToReceive = tasks_.size(); //TODO make work for failure
  //receiveDurationsOfTasksFromGroupMasters(0);

  if (doInitDSGUs) {
    // initialize dsgus
    initDsgus();
  }

  // return true if no group failed
  return !group_failed;
}

void ProcessManager::receiveDurationsOfTasksFromGroupMasters(size_t numDurationsToReceive = 0){
  if (numDurationsToReceive == 0){
    numDurationsToReceive = pgroups_.size();
  }
  for (size_t i = 0; i < numDurationsToReceive; ++i) {
    DurationInformation recvbuf;
    for([[maybe_unused]] const auto& t : pgroups_[i]->getTaskContainer()){
      // this assumes that the manager rank is the highest in globalComm
      MPIUtils::receiveClass(&recvbuf, i, theMPISystem()->getGlobalComm());

      const auto& levelVector = getLevelVectorFromTaskID(tasks_, recvbuf.task_id);
      if(LearningLoadModel* llm = dynamic_cast<LearningLoadModel*>(loadModel_.get())){
        llm->addDurationInformation(recvbuf, levelVector);
      }
      levelVectorToLastTaskDuration_[levelVector] = recvbuf.duration;
    }
  }
}

bool ProcessManager::runnext() {
  bool group_failed = waitAllFinished();

  assert(!group_failed && "runnext must not be called when there are failed groups");

  for (size_t i = 0; i < pgroups_.size(); ++i) {
    pgroups_[i]->runnext();
  }

  group_failed = waitAllFinished();
  //size_t numDurationsToReceive = tasks_.size(); //TODO make work for failure
  // if(!group_failed) {
    // receiveDurationsOfTasksFromGroupMasters(0);
  // }
  // return true if no group failed
  return !group_failed;
}

void ProcessManager::exit() {
  // wait until all process groups are in wait state
  // after sending the exit signal checking the status might not be possible
  size_t numWaiting = 0;

  while (numWaiting != pgroups_.size()) {
    numWaiting = 0;

    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (pgroups_[i]->getStatus() == PROCESS_GROUP_WAIT)
        ++numWaiting;
    }
  }

  // send exit signal to each group
  for (size_t i = 0; i < pgroups_.size(); ++i) {
    [[maybe_unused]] bool success = pgroups_[i]->exit();
    assert(success);
  }
}

void ProcessManager::initDsgus() {
  Stats::startEvent("manager init dsgus");
  // wait until all process groups are in wait state
  // after sending the exit signal checking the status might not be possible
  size_t numWaiting = 0;

  while (numWaiting != pgroups_.size()) {
    numWaiting = 0;

    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (pgroups_[i]->getStatus() == PROCESS_GROUP_WAIT) ++numWaiting;
    }
  }

  // tell groups to init Dsgus
  for (size_t i = 0; i < pgroups_.size(); ++i) {
    [[maybe_unused]] bool success = pgroups_[i]->initDsgus();
    assert(success);
  }

  waitAllFinished();
  Stats::stopEvent("manager init dsgus");
}

void ProcessManager::updateCombiParameters() {
  Stats::startEvent("manager update parameters");
  {
    [[maybe_unused]] bool fail = waitAllFinished();
    assert(!fail && "should not fail here");
  }

  for (auto g : pgroups_) g->updateCombiParameters(params_);
  {
    [[maybe_unused]] bool fail = waitAllFinished();
    assert(!fail && "should not fail here");
  }
  Stats::stopEvent("manager update parameters");
}

/*
 * Compute the group faults that occured at this combination step using the
 * fault simulator
 */
void ProcessManager::getGroupFaultIDs(std::vector<size_t>& faultsID,
                                      std::vector<ProcessGroupManagerID>& groupFaults) {
  for (auto p : pgroups_) {
    StatusType status = p->waitStatus();

    if (status == PROCESS_GROUP_FAIL) {
      TaskContainer failedTasks = p->getTaskContainer();
      groupFaults.push_back(p);
      for (auto task : failedTasks) faultsID.push_back(task->getID());
    }
  }
}

void ProcessManager::redistribute(std::vector<size_t>& taskID) {
  for (size_t i = 0; i < taskID.size(); ++i) {
    // find id in list of tasks
    Task* t = NULL;

    for (Task* tmp : tasks_) {
      if (tmp->getID() == taskID[i]) {
        t = tmp;
        break;
      }
    }

    assert(t != NULL);

    // wait for available process group
    ProcessGroupManagerID g = wait();

    // assign instance to group
    g->addTask(t);
  }

  size_t numWaiting = 0;

  while (numWaiting != pgroups_.size()) {
    numWaiting = 0;

    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (pgroups_[i]->getStatus() == PROCESS_GROUP_WAIT) ++numWaiting;
    }
  }

  std::cout << "Redistribute finished" << std::endl;
}

void ProcessManager::reInitializeGroup(std::vector<ProcessGroupManagerID>& recoveredGroups,
                                       std::vector<size_t>& tasksToIgnore) {
  std::vector<Task*> removeTasks;
  for (auto g : recoveredGroups) {
    // erase existing tasks in group members to avoid doubled tasks
    g->resetTasksWorker();
    for (Task* t : g->getTaskContainer()) {
      assert(t != NULL);
      if (std::find(tasksToIgnore.begin(), tasksToIgnore.end(), t->getID()) ==
          tasksToIgnore.end()) {  // ignore tasks that are recomputed
        StatusType status = g->waitStatus();
        std::cout << "status of g: " << status << "\n";
        // assign instance to group
        g->refreshTask(t);
      } else {
        if (std::find(removeTasks.begin(), removeTasks.end(), t) != removeTasks.end()) {
          std::cout << "Error task " << t->getID() << "twice in container! Processor"
                    << theMPISystem()->getWorldRank() << " \n";
        }
        removeTasks.push_back(t);
      }
    }
    for (Task* t : removeTasks) {
      g->removeTask(t);
    }
    removeTasks.clear();
  }

  size_t numWaiting = 0;

  while (numWaiting != pgroups_.size()) {
    numWaiting = 0;
    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (pgroups_[i]->getStatus() == PROCESS_GROUP_WAIT) ++numWaiting;
    }
  }

  std::cout << "Reinitialization finished" << std::endl;
}

void ProcessManager::recompute(std::vector<size_t>& taskID, bool failedRecovery,
                               std::vector<ProcessGroupManagerID>& recoveredGroups) {
  for (size_t i = 0; i < taskID.size(); ++i) {
    // find id in list of tasks
    Task* t = NULL;

    for (Task* tmp : tasks_) {
      if (tmp->getID() == taskID[i]) {
        t = tmp;
        break;
      }
    }

    assert(t != NULL);

    // wait for available process group
    if (failedRecovery) {
      ProcessGroupManagerID g = wait();
      // assign instance to group
      g->recompute(t);
    } else {
      ProcessGroupManagerID g = waitAvoid(recoveredGroups);
      // assign instance to group
      g->recompute(t);
    }
  }

  size_t numWaiting = 0;

  while (numWaiting != pgroups_.size()) {
    numWaiting = 0;

    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (pgroups_[i]->getStatus() == PROCESS_GROUP_WAIT) ++numWaiting;
    }
  }

  std::cout << "Recompute finished" << std::endl;
}

bool ProcessManager::recoverCommunicators(std::vector<ProcessGroupManagerID> failedGroups) {
  if (pgroups_.size() == failedGroups.size()) {
    std::cout << "last process groups failed! Aborting! \n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  waitAllFinished();

  // send recover communicators signal to alive groups
  for (ProcessGroupManagerID g : pgroups_) {
    if (g->getStatus() == PROCESS_GROUP_WAIT) {
      g->recoverCommunicators();
    }
  }

  bool failedRecovery = theMPISystem()->recoverCommunicators(true, failedGroups);

  // remove failed groups from group list and set new
  // todo: this is rather error prone. this relies on the previous functions
  // to have removed all processes of failed groups and that the order of
  // processes has not changed
  if (!failedRecovery) {
    for (auto pg : failedGroups) {
      pg->setStatus(PROCESS_GROUP_WAIT);
    }
  }
  if (failedRecovery) {
    pgroups_.erase(std::remove_if(pgroups_.begin(), pgroups_.end(),
                                  [](const ProcessGroupManagerID& p) {
                                    return (p->getStatus() == PROCESS_GROUP_FAIL);
                                  }),
                   pgroups_.end());

    for (size_t i = 0; i < pgroups_.size(); ++i) {
      pgroups_[i]->setMasterRank(int(i));
    }
  }
  return failedRecovery;
}


void ProcessManager::restoreCombischeme() {
  LevelVector lmin = params_.getLMin();
  LevelVector lmax = params_.getLMax();
  CombiMinMaxScheme combischeme(params_.getDim(), lmin, lmax);
  combischeme.createAdaptiveCombischeme();
  combischeme.makeFaultTolerant();
  std::vector<LevelVector> levels = combischeme.getCombiSpaces();
  std::vector<combigrid::real> coeffs = combischeme.getCoeffs();

  for (size_t i = 0; i < levels.size(); ++i) {
    params_.setCoeff(params_.getID(levels[i]), coeffs[i]);
  }

  this->updateCombiParameters();
}

bool ProcessManager::waitAllFinished() {
  bool group_failed = false;
  for (auto p : pgroups_) {
    if (waitForPG(p))
      group_failed = true;
  }
  return group_failed;
}

bool ProcessManager::waitForPG(ProcessGroupManagerID pg) {
  StatusType status = pg->waitStatus();
  if (status == PROCESS_GROUP_FAIL)
    return true;
  return false;
}

void ProcessManager::waitForAllGroupsToWait() const {
// wait until all process groups are in wait state
  size_t numWaiting = 0;

  while (numWaiting != pgroups_.size()) {
    numWaiting = 0;

    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (pgroups_[i]->getStatus() == PROCESS_GROUP_WAIT) ++numWaiting;
    }
  }
}

void ProcessManager::parallelEval(const LevelVector& leval, std::string& filename, size_t groupID) {
  // actually it would be enough to wait for the group which does the eval
  {
    [[maybe_unused]] bool fail = waitAllFinished();

    assert(!fail && "should not fail here");
  }

  assert(groupID < pgroups_.size());

  auto g = pgroups_[groupID];
  g->parallelEval(leval, filename);

  {
    [[maybe_unused]] bool fail = waitAllFinished();

    assert(!fail && "should not fail here");
  }
}

void ProcessManager::doDiagnostics(size_t taskID) {
  auto g = getProcessGroupWithTaskID(taskID);
  g->doDiagnostics(taskID);
  // call manager-side diagnostics on that Task
  for (auto task : tasks_) {
    if (task->getID() == taskID) {
      task->receiveDiagnostics();
      return;
    }
  }
}

std::map<size_t, double> ProcessManager::getLpNorms(int p) {
  std::map<size_t, double> norms;

  for (const auto& pg : pgroups_) {
    pg->getLpNorms(p, norms);
  }
  return norms;
}

double ProcessManager::getLpNorm(int p) {
  std::map<size_t, double> norms = this->getLpNorms(p);

  double norm = 0.0;
  double kahanTrailingTerm = 0.0;
  for (const auto& pair : norms) {
    auto summand = pair.second * this->params_.getCoeff(pair.first);
    // cf. https://en.wikipedia.org/wiki/Kahan_summation_algorithm
    auto y = summand - kahanTrailingTerm;
    auto t = norm + y;
    kahanTrailingTerm = (t - norm) - y;
    norm = t;
  }
  return norm;
}

void ProcessManager::setupThirdLevel() {
  Stats::startEvent("manager connect third level");
  std::string hostnameInfo = "manager = " + boost::asio::ip::host_name();
  std::cout << hostnameInfo << std::endl;
  thirdLevel_.connectToThirdLevelManager(10.);
  Stats::stopEvent("manager connect third level");
}

std::vector<CombiDataType> ProcessManager::interpolateValues(
    const std::vector<std::vector<real>>& interpolationCoords) {
  auto numValues = interpolationCoords.size();
  auto values = std::vector<CombiDataType>(numValues, std::numeric_limits<double>::quiet_NaN());
  MPI_Request request = MPI_REQUEST_NULL;

  // send interpolation coords as a single array
  std::vector<real> interpolationCoordsSerial = serializeInterpolationCoords(interpolationCoords);
  // have the last process group return the all-reduced values
  pgroups_[pgroups_.size() - 1]->interpolateValues(interpolationCoordsSerial, values, &request);
  for (size_t i = 0; i < pgroups_.size() - 1; ++i) {
    // all other groups only communicate the interpolation values to last group
    auto dummyValues = std::vector<CombiDataType>(0);
    pgroups_[i]->interpolateValues(interpolationCoordsSerial, dummyValues);
  }

  MPI_Wait(&request, MPI_STATUS_IGNORE);
  return values;
}

void ProcessManager::writeInterpolatedValuesPerGrid(
    const std::vector<std::vector<real>>& interpolationCoords, const std::string& filenamePrefix) {
  // send interpolation coords as a single array
  std::vector<real> interpolationCoordsSerial = serializeInterpolationCoords(interpolationCoords);

  for (size_t i = 0; i < pgroups_.size(); ++i) {
    pgroups_[i]->writeInterpolatedValuesPerGrid(interpolationCoordsSerial, filenamePrefix);
  }
}

void ProcessManager::writeInterpolatedValuesSingleFile(
    const std::vector<std::vector<real>>& interpolationCoords, const std::string& filenamePrefix) {
  Stats::startEvent("manager write interpolated");
  // send interpolation coords as a single array
  std::vector<real> interpolationCoordsSerial = serializeInterpolationCoords(interpolationCoords);

  // have the last process group write the all-reduced values
  auto dummyValuesEmpty = std::vector<CombiDataType>(0);
  pgroups_[pgroups_.size() - 1]->interpolateValues(interpolationCoordsSerial, dummyValuesEmpty,
                                                   nullptr, filenamePrefix);
  for (size_t i = 0; i < pgroups_.size() - 1; ++i) {
    // all other groups only communicate the interpolation values to last group
    pgroups_[i]->interpolateValues(interpolationCoordsSerial, dummyValuesEmpty);
  }
  Stats::stopEvent("manager write interpolated");
}

void ProcessManager::writeInterpolationCoordinates(
    const std::vector<std::vector<real>>& interpolationCoords, const std::string& filenamePrefix) const {
  std::string saveFilePath = filenamePrefix + "_coords.h5";
  h5io::writeValuesToH5File(interpolationCoords, saveFilePath, "manager", "coordinates");
}

void ProcessManager::monteCarloThirdLevel(size_t numPoints,
                                          std::vector<std::vector<real>>& coordinates,
                                          std::vector<CombiDataType>& values) {
  Stats::startEvent("manager MC third level");
  coordinates = montecarlo::getRandomCoordinates(numPoints, params_.getDim());
  auto ourCoordinatesSerial = serializeInterpolationCoords(coordinates);

  // obtain instructions from third level manager
  thirdLevel_.signalReadyToExchangeData();
  std::string instruction = thirdLevel_.fetchInstruction();

  // exchange coordinates with remote
  if (instruction == "send_first") {
    thirdLevel_.sendData(ourCoordinatesSerial.data(), ourCoordinatesSerial.size());
#ifndef NDEBUG
    // this part is redundant but also doesn't hurt (?)
    auto theirCoordinates = ourCoordinatesSerial; // to reserve the size
    thirdLevel_.recvData(theirCoordinates.data(), theirCoordinates.size());
    for (size_t i = 0; i < theirCoordinates.size(); ++i) {
      assert(ourCoordinatesSerial[i] == theirCoordinates[i]);
    }
#endif // !NDEBUG
  } else if (instruction == "recv_first") {
    thirdLevel_.recvData(ourCoordinatesSerial.data(), ourCoordinatesSerial.size());
#ifndef NDEBUG
    thirdLevel_.sendData(ourCoordinatesSerial.data(), ourCoordinatesSerial.size());
#endif // !NDEBUG
    coordinates = deserializeInterpolationCoords(ourCoordinatesSerial, params_.getDim());
  } else {
    throw std::runtime_error("unknown instruction: " + instruction);
  }
  thirdLevel_.signalReady();

  // interpolate locally
  values = this->interpolateValues(coordinates);

  // obtain instructions from third level manager
  thirdLevel_.signalReadyToExchangeData();
  instruction = thirdLevel_.fetchInstruction();

  // exchange values with remote
  auto buffSize = numPoints;
  std::vector<CombiDataType> remoteValues(buffSize);
  if (instruction == "send_first") {
    thirdLevel_.sendData(values.data(), buffSize);
    thirdLevel_.recvData(remoteValues.data(), buffSize);
  } else if (instruction == "recv_first") {
    thirdLevel_.recvData(remoteValues.data(), buffSize);
    thirdLevel_.sendData(values.data(), buffSize);
  }
  thirdLevel_.signalReady();

  // add them up
  for (size_t i = 0; i < numPoints; ++i){
    values[i] += remoteValues[i];
  }
  Stats::stopEvent("manager MC third level");
}

void ProcessManager::writeSparseGridMinMaxCoefficients(const std::string& filename) {
  pgroups_.back()->writeSparseGridMinMaxCoefficients(filename);
}

void ProcessManager::reschedule() {
  std::map<LevelVector, int> levelVectorToProcessGroupIndex;
  for (size_t i = 0; i < pgroups_.size(); ++i) {
    for (const auto& t : pgroups_[i]->getTaskContainer()) {
      levelVectorToProcessGroupIndex.insert({t->getLevelVector(), i});
    }
  }
  auto tasksToMigrate = rescheduler_->eval(levelVectorToProcessGroupIndex,
                                           levelVectorToLastTaskDuration_,
                                           loadModel_.get());
  for (const auto& t : tasksToMigrate) {
    auto levelvectorToMigrate = t.first;
    auto processGroupIndexToAddTaskTo = t.second;
    auto processGroupIndexToRemoveTaskFrom =
      levelVectorToProcessGroupIndex.at(levelvectorToMigrate);

    Task *removedTask =
      pgroups_[processGroupIndexToRemoveTaskFrom]->rescheduleRemoveTask(
          levelvectorToMigrate);
    waitAllFinished();
    assert(removedTask != nullptr);
    pgroups_[processGroupIndexToAddTaskTo]->rescheduleAddTask(removedTask);
    waitAllFinished();
  }

  // update local tasks_ vector!
  tasks_.clear();
  for (auto& pg : pgroups_) {
    for (auto t : pg->getTaskContainer()) {
      tasks_.push_back(t);
    }
  }
}

void ProcessManager::writeCombigridsToVTKPlotFile(ProcessGroupManagerID pg) {
#if defined(USE_VTK)
  pg->writeCombigridsToVTKPlotFile();
  waitForPG(pg);
#else
  std::cout << "Warning: no vtk output produced as DisCoTec was compiled without VTK." << std::endl;
#endif /* defined(USE_VTK) */
}

void ProcessManager::writeDSGsToDisk(const std::string& filenamePrefix) {
  thirdLevelPGroup_->writeDSGsToDisk(filenamePrefix);
  waitForPG(thirdLevelPGroup_);
}

void ProcessManager::readDSGsFromDisk(const std::string& filenamePrefix) {
  thirdLevelPGroup_->readDSGsFromDisk(filenamePrefix);
  waitForPG(thirdLevelPGroup_);
}

} /* namespace combigrid */
