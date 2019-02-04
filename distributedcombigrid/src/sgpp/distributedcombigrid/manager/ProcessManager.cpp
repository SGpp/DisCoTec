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
  thirdLevelPGroup_(pgroups[0]), // TODO: changing PG needs adjustment in MPISystem
  tasks_(tasks),
  params_(params),
  thirdLevel_(params.getThirdLevelHost(), params.getThirdLevelPort(), params.getSystemName())
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

void ProcessManager::setupThirdLevel()
{
  thirdLevel_.connectToThirdLevelManager();
}


} /* namespace combigrid */

