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
    // Commented out since not relevant when SDC occurs
//    assert( !fail && "should not fail here" );
  }

  for( auto g : pgroups_ )
    g->updateCombiParameters(params_);
  std::cout<<"Updated all.";
  {
    bool fail = waitAllFinished();
    // Commented out since not relevant when SDC occurs
//    assert( !fail && "should not fail here" );
  }
  std::cout<<"Waiting done";
}

/*
 * Compute the group faults that occurred at this combination step using the fault simulator
 */
void
ProcessManager::getGroupFaultIDs( std::vector<int>& faultsID ) {
  for( auto p : pgroups_ ){
    StatusType status = p->waitStatus();

    if( status == PROCESS_GROUP_FAIL ){
      TaskContainer failedTasks = p->getTaskContainer();
      for( auto task : failedTasks )
        faultsID.push_back(task->getID());
    }
  }
}

void
ProcessManager::getSDCFaultIDs( std::vector<int>& faultsID ) {
  for( auto p : pgroups_ ){
    StatusType status = p->waitStatus();
    if( status == PROCESS_GROUP_FAIL ){
      size_t numTasks = p->getTaskContainer().size();
      faultsID.resize( numTasks );
      MPI_Status statusSDC;
      MPI_Recv( &faultsID[0], numTasks , MPI_INT, MPI_ANY_SOURCE, infoTag, theMPISystem()->getGlobalComm(), &statusSDC );
      int numTasksSDC;
      MPI_Get_count( &statusSDC, MPI_INT, &numTasksSDC );
      faultsID.resize( numTasksSDC );
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
  std::cout<<"Waiting...";
  while (numWaiting != pgroups_.size()) {
    numWaiting = 0;

    for (size_t i = 0; i < pgroups_.size(); ++i) {
      if (pgroups_[i]->getStatus() == PROCESS_GROUP_WAIT)
        ++numWaiting;
      std::cout<<"Status = "<<pgroups_[i]->getStatus();
    }

  }

  std::cout << "Recompute finished" << std::endl;
}


void ProcessManager::recoverCommunicators(){
  waitAllFinished();

  // send recover communicators signal to alive groups
  for( ProcessGroupManagerID g : pgroups_ ){
    if( g->getStatus() == PROCESS_GROUP_WAIT ){
      g->recoverCommunicators();
    }
  }

//  theMPISystem()->recoverCommunicators( true );

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
}

void ProcessManager::recover(){

  std::vector<int> faultsID;
  getGroupFaultIDs(faultsID);

  /* call optimization code to find new coefficients */
  const std::string prob_name = "interpolation based optimization";
  std::vector<int> redistributeFaultsID, recomputeFaultsID;
  recomputeOptimumCoefficients(prob_name, faultsID, redistributeFaultsID, recomputeFaultsID);

  /* recover communicators*/
  recoverCommunicators();

  /* redistribute failed tasks to living groups */
  redistribute(faultsID);

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
  for( auto p : pgroups_ ){
    StatusType status = p->waitStatus();
    if( status == PROCESS_GROUP_FAIL ){
      group_failed = true;
    }
  }

  return group_failed;
}

} /* namespace combigrid */

