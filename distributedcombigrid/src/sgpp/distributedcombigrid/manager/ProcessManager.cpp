/*
 * ProcessManager.cpp
 *
 *  Created on: Oct 8, 2013
 *      Author: heenemo
 */

#include <algorithm>
#include <iostream>
#include "sgpp/distributedcombigrid/manager/ProcessManager.hpp"
#include "sgpp/distributedcombigrid/combicom/CombiCom.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupManager.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"

namespace combigrid {

bool compareInstances(const Task* instance1, const Task* instance2) {
  return (instance1->estimateRuntime() > instance2->estimateRuntime());
}

ProcessManager::ProcessManager( ProcessGroupManagerContainer& pgroups,
                                TaskContainer& tasks,
                                CombiParameters& params ) :
  pgroups_(pgroups),
  tasks_(tasks),
  params_(params)
{
}

ProcessManager::~ProcessManager() {
}

bool ProcessManager::runfirst() {
  // sort instances in decreasing order
  std::sort(tasks_.begin(), tasks_.end(), compareInstances);

  for (size_t i = 0; i < tasks_.size(); ++i) {
    // wait for available process group
    ProcessGroupManagerID g = wait();

    // assign instance to group
    g->runfirst(tasks_[i]);
  }

  bool group_failed = waitAllFinished();
  // return true if no group failed
  return !group_failed;
}

bool ProcessManager::runnext() {
  bool group_failed = waitAllFinished();

  assert( !group_failed
      && "runnext must not be called when there are failed groups" );

  for (size_t i = 0; i < pgroups_.size(); ++i) {
    pgroups_[i]->runnext();
  }

  group_failed = waitAllFinished();
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
    bool success = pgroups_[i]->exit();
    assert(success);
  }
}

void ProcessManager::updateCombiParameters() {
  {
    bool fail = waitAllFinished();

    assert( !fail && "should not fail here" );
  }

  for( auto g : pgroups_ )
    g->updateCombiParameters(params_);

  {
    bool fail = waitAllFinished();

    assert( !fail && "should not fail here" );
  }
}

/*
 * Compute the group faults that occured at this combination step using the fault simulator
 */
void ProcessManager::getGroupFaultIDs( std::vector< int>& faultsID, std::vector< ProcessGroupManagerID>& groupFaults ) {
  for( auto p : pgroups_ ){
    StatusType status = p->waitStatus();

    if( status == PROCESS_GROUP_FAIL ){
      TaskContainer failedTasks = p->getTaskContainer();
      groupFaults.push_back(p);
      for( auto task : failedTasks )
        faultsID.push_back(task->getID());
    }
  }
}


void ProcessManager::redistribute( std::vector<int>& taskID ) {
  for (size_t i = 0; i < taskID.size(); ++i) {
    // find id in list of tasks
    Task* t = NULL;

    for ( Task* tmp : tasks_ ) {
      if ( tmp->getID() == taskID[i] ) {
        t = tmp;
        break;
      }
    }

    assert( t != NULL );

    // wait for available process group
    ProcessGroupManagerID g = wait();

    // assign instance to group
    g->addTask( t );
  }

  size_t numWaiting = 0;

  while (numWaiting != pgroups_.size()) {
    numWaiting = 0;

    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (pgroups_[i]->getStatus() == PROCESS_GROUP_WAIT)
        ++numWaiting;
    }

  }

  std::cout << "Redistribute finished" << std::endl;
}

void ProcessManager::reInitializeGroup(std::vector< ProcessGroupManagerID>& recoveredGroups ) {
  for (auto g : recoveredGroups) {
    //erase existing tasks in group members to avoid doubled tasks
    g->resetTasksWorker();
    for ( Task* t : g->getTaskContainer()) {
      assert( t != NULL );
      // assign instance to group
      g->addTask( t );
    }
  }

  size_t numWaiting = 0;

  while (numWaiting != pgroups_.size()) {
    numWaiting = 0;
    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (pgroups_[i]->getStatus() == PROCESS_GROUP_WAIT)
        ++numWaiting;
    }
  }

  std::cout << "Reinitialization finished" << std::endl;
}


void ProcessManager::recompute( std::vector<int>& taskID ) {
  for (size_t i = 0; i < taskID.size(); ++i) {
    // find id in list of tasks
    Task* t = NULL;

    for ( Task* tmp : tasks_ ) {
      if ( tmp->getID() == taskID[i] ) {
        t = tmp;
        break;
      }
    }

    assert( t != NULL );

    // wait for available process group
    ProcessGroupManagerID g = wait();

    // assign instance to group
    g->recompute( t );
  }

  size_t numWaiting = 0;

  while (numWaiting != pgroups_.size()) {
    numWaiting = 0;

    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (pgroups_[i]->getStatus() == PROCESS_GROUP_WAIT)
        ++numWaiting;
    }

  }

  std::cout << "Recompute finished" << std::endl;
}


bool ProcessManager::recoverCommunicators(std::vector< ProcessGroupManagerID> failedGroups){
  waitAllFinished();

  // send recover communicators signal to alive groups
  for( ProcessGroupManagerID g : pgroups_ ){
    if( g->getStatus() == PROCESS_GROUP_WAIT ){
      g->recoverCommunicators();
    }
  }

  bool failedRecovery = theMPISystem()->recoverCommunicators( true,failedGroups );

  // remove failed groups from group list and set new
  // todo: this is rather error prone. this relies on the previous functions
  // to have removed all processes of failed groups and that the order of
  // processes has not changed
  pgroups_.erase(
    std::remove_if( pgroups_.begin(),
                    pgroups_.end(),
                    [] (const ProcessGroupManagerID& p) {
                       return (p->getStatus() == PROCESS_GROUP_FAIL); } ),
    pgroups_.end() );

  for( size_t i=0; i<pgroups_.size(); ++i ){
    pgroups_[i]->setMasterRank( int(i) );
  }
  return failedRecovery;
}

void ProcessManager::recover(){ //outdated

  std::vector<int> faultsID;
  std::vector< ProcessGroupManagerID> groupFaults;
  getGroupFaultIDs(faultsID, groupFaults);

  /* call optimization code to find new coefficients */
  const std::string prob_name = "interpolation based optimization";
  std::vector<int> redistributeFaultsID, recomputeFaultsID;
  recomputeOptimumCoefficients(prob_name, faultsID, redistributeFaultsID, recomputeFaultsID);
  //time does not need to be updated in gene but maybe in other applications
  /*for ( auto id : redistributeFaultsID ) {
    GeneTask* tmp = static_cast<GeneTask*>(manager.getTask(id));
    tmp->setStepsTotal((i+1)*nsteps);
    tmp->setCombiStep(i+1);
  }

  for ( auto id : recomputeFaultsID ) {
    GeneTask* tmp = static_cast<GeneTask*>(manager.getTask(id));
    tmp->setStepsTotal((i)*nsteps);
    tmp->setCombiStep(i);
  }*/
  /* recover communicators*/
  bool failedRecovery = recoverCommunicators(groupFaults);

  /* redistribute failed tasks to living groups */
  //redistribute(faultsID);
  if(failedRecovery){
    redistribute(redistributeFaultsID);
  }
  else{
    reInitializeGroup(groupFaults);
  }
  /* communicate new combination scheme*/
  updateCombiParameters();

  /* if some tasks have to be recomputed, do so*/
  if(!recomputeFaultsID.empty())
    recompute(recomputeFaultsID);
}

void ProcessManager::restoreCombischeme() {

  LevelVector lmin = params_.getLMin();
  LevelVector lmax = params_.getLMax();
  CombiMinMaxScheme combischeme(params_.getDim(), lmin, lmax);
  combischeme.createAdaptiveCombischeme();
  combischeme.makeFaultTolerant();
  std::vector<LevelVector> levels = combischeme.getCombiSpaces();
  std::vector<combigrid::real> coeffs = combischeme.getCoeffs();

  for (size_t i = 0; i < levels.size(); ++i){
    params_.setCoeff( params_.getID(levels[i]), coeffs[i] );
  }
}

bool ProcessManager::waitAllFinished(){
  bool group_failed = false;
  int i = 0;
  for( auto p : pgroups_ ){
    StatusType status = p->waitStatus();
    if( status == PROCESS_GROUP_FAIL ){
      group_failed = true;
    }
    ++i;
  }

  return group_failed;
}


void ProcessManager::parallelEval( const LevelVector& leval,
                                   std::string& filename,
                                   size_t groupID ){
  // actually it would be enough to wait for the group which does the eval
  {
    bool fail = waitAllFinished();

    assert( !fail && "should not fail here" );
  }

  assert( groupID < pgroups_.size() );

  auto g = pgroups_[ groupID ];
  g->parallelEval( leval, filename );

  {
    bool fail = waitAllFinished();

    assert( !fail && "should not fail here" );
  }
}

} /* namespace combigrid */

