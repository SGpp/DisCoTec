#include "sgpp/distributedcombigrid/manager/ProcessManager.hpp"
#include <algorithm>
#include <iostream>
#include "sgpp/distributedcombigrid/combicom/CombiCom.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/mpi/MPIUtils.hpp"

namespace combigrid {

ProcessManager::~ProcessManager() {}

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

bool ProcessManager::runfirst() {
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
  receiveDurationsOfTasksFromGroupMasters(0);

  // return true if no group failed
  return !group_failed;
}

void ProcessManager::receiveDurationsOfTasksFromGroupMasters(size_t numDurationsToReceive = 0){
  if (numDurationsToReceive == 0){
    numDurationsToReceive = pgroups_.size();
  }
  for (size_t i = 0; i < numDurationsToReceive; ++i) {
    DurationInformation recvbuf;
    for(const auto& t : pgroups_[i]->getTaskContainer()){
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
  if(!group_failed)
    receiveDurationsOfTasksFromGroupMasters(0);
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
      if (pgroups_[i]->getStatus() == PROCESS_GROUP_WAIT) ++numWaiting;
    }
  }

  // send exit signal to each group
  for (size_t i = 0; i < pgroups_.size(); ++i) {
    bool success = pgroups_[i]->exit();
    assert(success);
  }
}

void ProcessManager::updateCombiParameters() {
  {
    bool fail = waitAllFinished();
    assert(!fail && "should not fail here");
  }

  for (auto g : pgroups_) g->updateCombiParameters(params_);
  {
    bool fail = waitAllFinished();
    assert(!fail && "should not fail here");
  }
}

/*
 * Compute the group faults that occured at this combination step using the fault simulator
 */
void ProcessManager::getGroupFaultIDs(std::vector<int>& faultsID,
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

void ProcessManager::redistribute(std::vector<int>& taskID) {
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
                                       std::vector<int>& tasksToIgnore) {
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

void ProcessManager::recompute(std::vector<int>& taskID, bool failedRecovery,
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

void ProcessManager::recover(int i, int nsteps) {  // outdated

  std::vector<int> faultsID;
  std::vector<ProcessGroupManagerID> groupFaults;
  getGroupFaultIDs(faultsID, groupFaults);

  /* call optimization code to find new coefficients */
  const std::string prob_name = "interpolation based optimization";
  std::vector<int> redistributeFaultsID, recomputeFaultsID;
  recomputeOptimumCoefficients(prob_name, faultsID, redistributeFaultsID, recomputeFaultsID);
  // time does not need to be updated in gene but maybe in other applications
  /*  for ( auto id : redistributeFaultsID ) {
      GeneTask* tmp = static_cast<GeneTask*>(getTask(id));
      tmp->setStepsTotal((i+1)*nsteps);
      tmp->setCombiStep(i+1);
    }

    for ( auto id : recomputeFaultsID ) {
      GeneTask* tmp = static_cast<GeneTask*>(getTask(id));
      tmp->setStepsTotal((i)*nsteps);
      tmp->setCombiStep(i);
    }*/
  /* recover communicators*/
  bool failedRecovery = recoverCommunicators(groupFaults);
  /* communicate new combination scheme*/
  if (failedRecovery) {
    std::cout << "redistribute \n";
    redistribute(redistributeFaultsID);
  } else {
    std::cout << "reinitializing group \n";
    reInitializeGroup(groupFaults, recomputeFaultsID);
  }

  /* if some tasks have to be recomputed, do so*/
  if (!recomputeFaultsID.empty()) {
    recompute(recomputeFaultsID, failedRecovery, groupFaults);
  }
  std::cout << "updateing Combination Parameters \n";
  // needs to be after reInitialization!
  updateCombiParameters();
  /* redistribute failed tasks to living groups */
  // redistribute(faultsID);
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
}

bool ProcessManager::waitAllFinished() {
  bool group_failed = false;
  for (auto p : pgroups_) {
    StatusType status = p->waitStatus();
    if (status == PROCESS_GROUP_FAIL) {
      group_failed = true;
    }
  }

  return group_failed;
}

void ProcessManager::parallelEval(const LevelVector& leval, std::string& filename, size_t groupID) {
  // actually it would be enough to wait for the group which does the eval
  {
    bool fail = waitAllFinished();

    assert(!fail && "should not fail here");
  }

  assert(groupID < pgroups_.size());

  auto g = pgroups_[groupID];
  g->parallelEval(leval, filename);

  {
    bool fail = waitAllFinished();

    assert(!fail && "should not fail here");
  }
}

std::map<int, double> ProcessManager::getLpNorms(int p) {
  SignalType signal;
  if (p == 2) {
    signal = GET_L2_NORM;
  } else if (p == 1) {
    signal = GET_L1_NORM;
  } else if (p == 0) {
    signal = GET_MAX_NORM;
  } else {
    assert(false && "please implement signal for this norm");
  }

  for (const auto& pg : pgroups_) {
    pg->sendSignalToProcessGroup(signal);
  }

  std::map<int, double> norms;
  for (RankType i = 0; i < pgroups_.size(); ++i) {
    std::vector<double> recvbuf;
    auto numTasks = pgroups_[i]->getTaskContainer().size();
    recvbuf.resize(numTasks);

    // this assumes that the manager rank is the highest in globalComm
    MPI_Recv(recvbuf.data(), static_cast<int>(numTasks), MPI_DOUBLE, i, TRANSFER_NORM_TAG,
             theMPISystem()->getGlobalComm(), MPI_STATUS_IGNORE);

    for (size_t j = 0; j < numTasks; ++j) {
      norms[pgroups_[i]->getTaskContainer()[j]->getID()] = recvbuf[j];
    }
  }

  for (const auto& pg : pgroups_) {
    pg->setProcessGroupBusyAndReceive();
  }
  return norms;
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

} /* namespace combigrid */
