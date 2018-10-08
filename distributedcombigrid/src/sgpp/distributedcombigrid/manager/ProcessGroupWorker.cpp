/*
 * ProcessGroupWorker.cpp
 *
 *  Created on: Jun 24, 2014
 *      Author: heenemo
 */

#include "sgpp/distributedcombigrid/manager/ProcessGroupWorker.hpp"

#include "boost/lexical_cast.hpp"

#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupSignals.hpp"
#include "sgpp/distributedcombigrid/combicom/CombiCom.hpp"
#include "sgpp/distributedcombigrid/hierarchization/DistributedHierarchization.hpp"
#include "sgpp/distributedcombigrid/mpi/MPIUtils.hpp"
#include "sgpp/distributedcombigrid/utils/Config.hpp"

#include "sgpp/distributedcombigrid/mpi_fault_simulator/MPI-FT.h"
#include <string>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <thread>

namespace combigrid {

ProcessGroupWorker::ProcessGroupWorker() :
    currentTask_( NULL),
    status_(PROCESS_GROUP_WAIT),
    combinedFG_( NULL),
    combinedUniDSGVector_(0),
    combinedFGexists_(false),
    combiParameters_(),
	combiScheme_ (LevelVector {1}, LevelVector {1}),
    combiParametersSet_(false),
    currentCombi_(0)
{
    t_fault_ = -1;
    startTimeIteration_ = (std::chrono::high_resolution_clock::now());
  MASTER_EXCLUSIVE_SECTION {
    std::string fname ="out/all-betas-"+std::to_string(theMPISystem()->getGlobalRank())+".txt";
    //betasFile_ = std::ofstream( fname, std::ofstream::out );
  }
}

ProcessGroupWorker::~ProcessGroupWorker() {
  delete combinedFG_;
}

void ProcessGroupWorker::runFirst() {
	Task* t;
	// local root receives task
	MASTER_EXCLUSIVE_SECTION{
	Task::receive( &t,
			theMPISystem()->getManagerRank(),
			theMPISystem()->getGlobalComm() );
}
	std::cout << "running task: " <<  t->getID() << " ";
// broadcast task to other process of pgroup
	Task::broadcast(&t, theMPISystem()->getMasterRank(),
			theMPISystem()->getLocalComm());
	MPI_Barrier(theMPISystem()->getLocalComm());
	// add task to task storage
	tasks_.push_back(t);
	status_ = PROCESS_GROUP_BUSY;
	// set currentTask
	currentTask_ = tasks_.back();
	std::cout << "current task1 " <<  currentTask_->getID() << " ";
	// initalize task
	Stats::startEvent("worker init");
	currentTask_->init(theMPISystem()->getLocalComm());
	t_fault_ = currentTask_->initFaults(t_fault_, startTimeIteration_);
	std::cout << "current task2 " <<  currentTask_->getID() << " ";
	Stats::stopEvent("worker init");
	// execute task
	Stats::startEvent("worker run first");
	currentTask_->run(theMPISystem()->getLocalComm());
	std::cout << "current task3 " <<  currentTask_->getID() << " ";
	Stats::stopEvent("worker run first");
}

void ProcessGroupWorker::runNewTask(){
	std::this_thread::sleep_for(std::chrono::seconds {1});
	std::cout << "dummytest\n";
	Task* t;
	// local root receives task
	MASTER_EXCLUSIVE_SECTION{
	Task::receive( &t,
			theMPISystem()->getManagerRank(),
			theMPISystem()->getGlobalComm() );
}
// broadcast task to other process of pgroup
	Task::broadcast(&t, theMPISystem()->getMasterRank(),
			theMPISystem()->getLocalComm());
	MPI_Barrier(theMPISystem()->getLocalComm());
	// add task to task storage
	tasks_.push_back(t);
	status_ = PROCESS_GROUP_BUSY;
	// set currentTask
	currentTask_ = tasks_.back();
	// initalize task
	Stats::startEvent("worker init");
	currentTask_->init(theMPISystem()->getLocalComm());
	t_fault_ = currentTask_->initFaults(t_fault_, startTimeIteration_);
	Stats::stopEvent("worker init");

	if(!isGENE){
		setCombinedSolutionUniform(currentTask_);
	}

	// execute task
	Stats::startEvent("worker run first");
	currentTask_->run(theMPISystem()->getLocalComm());
	Stats::stopEvent("worker run first");
}

void ProcessGroupWorker::runNext() {
	// this should not happen
	//assert(tasks_.size() > 0);
	// reset finished status of all tasks
	if (tasks_.size() != 0) {
		for (size_t i = 0; i < tasks_.size(); ++i)
			tasks_[i]->setFinished(false);
		status_ = PROCESS_GROUP_BUSY;
		// set currentTask
		currentTask_ = tasks_[0];
		// run first task
		if (!isGENE) {
			Stats::startEvent("worker run");
		}
		currentTask_->run(theMPISystem()->getLocalComm());
		if (!isGENE) {
			Stats::stopEvent("worker run");
		}
	} else {
		std::cout << "Possible error: No tasks! \n";
	}
}

void ProcessGroupWorker::addTask() {
	//add a new task to the process group
	std::cout << "adding a single task" << std::endl;
	Task* t;
	// local root receives task
	MASTER_EXCLUSIVE_SECTION{
	Task::receive(&t, theMPISystem()->getManagerRank(), theMPISystem()->getGlobalComm());

	std::cout << "received task" << std::endl;
	}
// broadcast task to other process of pgroup
	Task::broadcast(&t, 0,
			theMPISystem()->getLocalComm());
	std::cout << "added task id: " << t->getID();
	MPI_Barrier(theMPISystem()->getLocalComm());
	// check if task already exists on this group
	for (auto tmp : tasks_)
		assert(tmp->getID() != t->getID());
	// initalize task and set values to zero
	// the task will get the proper initial solution during the next combine
	t->init(theMPISystem()->getLocalComm());
	t_fault_ = t->initFaults(t_fault_, startTimeIteration_);
	t->setZero();
	t->setFinished(true);
	// add task to task storage
	tasks_.push_back(t);
	currentTask_ = tasks_.back(); //important for updating values
	currentTask_->changeDir(theMPISystem()->getLocalComm());
	status_ = PROCESS_GROUP_BUSY;
}

void ProcessGroupWorker::recompute() {
	//recompute the received task (immediately computes tasks -> difference to ADD_TASK)
	Task* t;
	// local root receives task
	MASTER_EXCLUSIVE_SECTION{
	Task::receive(&t, theMPISystem()->getManagerRank(), theMPISystem()->getGlobalComm());
}
// broadcast task to other process of pgroup
	Task::broadcast(&t, theMPISystem()->getMasterRank(),
			theMPISystem()->getLocalComm());
	MPI_Barrier(theMPISystem()->getLocalComm());
	// add task to task storage
	tasks_.push_back(t);
	status_ = PROCESS_GROUP_BUSY;
	// set currentTask
	currentTask_ = tasks_.back();
	// initalize task
	currentTask_->init(theMPISystem()->getLocalComm());
	t_fault_ = currentTask_->initFaults(t_fault_, startTimeIteration_);
	currentTask_->setZero();
	// fill task with combisolution
	if (!isGENE) {
		setCombinedSolutionUniform(currentTask_);
	}
	// execute task
	currentTask_->run(theMPISystem()->getLocalComm());
}

SignalType ProcessGroupWorker::wait() {
  if(status_ == PROCESS_GROUP_FAIL){ //in this case worker got reused
    status_ = PROCESS_GROUP_WAIT;
  }
  if (status_ != PROCESS_GROUP_WAIT){
    int myRank;
    MPI_Comm_rank(theMPISystem()->getWorldComm(), &myRank);
#ifdef DEBUG_OUTPUT
    std::cout << "status is " << status_ << "of rank " << myRank << "\n";
    std::cout << "executing next task\n";
#endif
    return RUN_NEXT;
  }
  SignalType signal = -1;

  MASTER_EXCLUSIVE_SECTION {
    // receive signal from manager
    MPI_Recv( &signal, 1, MPI_INT,
              theMPISystem()->getManagerRank(),
              signalTag,
              theMPISystem()->getGlobalComm(),
              MPI_STATUS_IGNORE);
  }
  // distribute signal to other processes of pgroup
  MPI_Bcast( &signal, 1, MPI_INT,
             theMPISystem()->getMasterRank(),
             theMPISystem()->getLocalComm() );
#ifdef DEBUG_OUTPUT
  std::cout << theMPISystem()->getWorldRank() << " waits for signal " << signal << " \n";
#endif
  // process signal
  if (signal == RUN_FIRST) {
		runFirst();
  } else if (signal == RUN_NEXT) {
    // this should not happen
    //assert(tasks_.size() > 0);
    // reset finished status of all tasks
		runNext();

  } else if (signal == RUN_NEWTASK){
	  	runNewTask();
  } else if (signal == ADD_TASK) { //add a new task to the process group
		addTask();
  } else if (signal == RESET_TASKS) { //deleta all tasks (used in process recovery)
    std::cout << "resetting tasks" << std::endl;


    // freeing tasks
    for ( auto tmp : tasks_ )
      delete(tmp);

    tasks_.clear();
    status_ = PROCESS_GROUP_BUSY;

  } else if (signal == EVAL) {
    // receive x

    // loop over all tasks
    // t.eval(x)
  } else if (signal == EXIT) {
    if(isGENE){
      chdir( "../ginstance" );
    }

  } else if (signal == SYNC_TASKS) {
    MASTER_EXCLUSIVE_SECTION {
      for (size_t i = 0; i < tasks_.size(); ++i) {
        Task::send(&tasks_[i], theMPISystem()->getManagerRank(), theMPISystem()->getGlobalComm());
      }
    }
  } else if (signal == COMBINE) { //start combination

    Stats::startEvent("combine");
    combineUniform();
    currentCombi_++;
    Stats::stopEvent("combine");

  } else if (signal == GRID_EVAL) { // not supported anymore

    Stats::startEvent("eval");
    gridEval();
    Stats::stopEvent("eval");

    return signal;

  } else if (signal == COMBINE_FG) {

    combineFG();

  } else if (signal == UPDATE_COMBI_PARAMETERS) { //update combiparameters (e.g. in case of faults -> FTCT)

    updateCombiParameters();

  } else if (signal == RECOMPUTE) { //recompute the received task (immediately computes tasks -> difference to ADD_TASK)
	recompute();
  } else if ( signal ==  RECOVER_COMM ){ //start recovery in case of faults
    theMPISystem()->recoverCommunicators( true );
    return signal;
  } else if( signal == PARALLEL_EVAL ){ //output final grid

    Stats::startEvent("parallel eval");
    parallelEval();
    Stats::stopEvent("parallel eval");

  } else if (signal == TASK_TO_PROC){
	  MASTER_EXCLUSIVE_SECTION {
	  MPIUtils::receiveClass(&taskToProc_, theMPISystem()->getManagerRank(), theMPISystem()->getGlobalComm());

	  }
	  //TODO rank
	  MPIUtils::broadcastClass(&taskToProc_, 0, theMPISystem()->getLocalComm());

	  MASTER_EXCLUSIVE_SECTION {
		  std::cout << "finished taskToProc bcast\n";
	  }
  } else if (signal == ADD_EXPANSION){
	  constexpr int addExpansionTag = 1235;
	  LevelVector vec (combiScheme_.dim());

	  MASTER_EXCLUSIVE_SECTION {
	  MPI_Recv(vec.data(), vec.size(), MPI_LONG, theMPISystem()->getManagerRank(),
			  addExpansionTag, theMPISystem()->getGlobalComm(), MPI_STATUS_IGNORE);

	  std::cout << "master received: " << vec << "\n";
	  }


	  MPI_Bcast(vec.data(), vec.size(), MPI_LONG, 0, theMPISystem()->getLocalComm());
	  std::cout << "added Expansion: " << theMPISystem()->getWorldRank()<< "\n";
	  combiScheme_.addExpansion(vec);
  } else if (signal == BEST_EXPANSION){
	  std::cout << "Worker Start expansion\n";
	  findBestExpansion();
  }
  if(isGENE){
    // special solution for GENE
    // todo: find better solution and remove this
	if(signal == RUN_FIRST  || signal == RUN_NEXT || signal == RECOMPUTE || signal == RUN_NEWTASK)
	  std::cout << "isFinished: " << currentTask_->isFinished() << "\n";
    if( ( signal == RUN_FIRST  || signal == RUN_NEXT || signal == RECOMPUTE || signal == RUN_NEWTASK) && /*!currentTask_->isFinished() &&*/ omitReadySignal )
      return signal;
  }
  // in the general case: send ready signal.
  //if(!omitReadySignal)
  std::cout << "readCalled with signal: " << signal << "\n";
  ready();
  if(isGENE){
    if(signal == ADD_TASK){ //ready resets currentTask but needs to be set for GENE
      currentTask_ = tasks_.back();
    }
  }
  return signal;
}
void ProcessGroupWorker::decideToKill(){
  //decide if processor was killed during this iteration
  currentTask_->decideToKill();
}

void ProcessGroupWorker::ready() {
  if( ENABLE_FT ){
    // with this barrier the local root but also each other process can detect
    // whether a process in the group has failed
    int globalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);
    //std::cout << "rank " << globalRank << " is ready \n";
    int err = simft::Sim_FT_MPI_Barrier( theMPISystem()->getLocalCommFT() );

    if( err == MPI_ERR_PROC_FAILED ){
      status_ = PROCESS_GROUP_FAIL;


      std::cout << "rank " << globalRank << " fault detected" << std::endl;
    }
  }

  if( status_ != PROCESS_GROUP_FAIL ){
    // check if there are unfinished tasks
    for (size_t i = 0; i < tasks_.size(); ++i) {
      if (!tasks_[i]->isFinished()) {
        status_ = PROCESS_GROUP_BUSY;

        // set currentTask
        currentTask_ = tasks_[i];
        Stats::startEvent("worker run");
        currentTask_->run(theMPISystem()->getLocalComm());
        Stats::stopEvent("worker run");
 	if( ENABLE_FT ){
          // with this barrier the local root but also each other process can detect
          // whether a process in the group has failed
          int err = simft::Sim_FT_MPI_Barrier( theMPISystem()->getLocalCommFT() );

          if( err == MPI_ERR_PROC_FAILED ){
            status_ = PROCESS_GROUP_FAIL;
            break;
          }
        }
	//merge problem?
	// todo: gene specific voodoo 
 	      if(isGENE && !currentTask_->isFinished()){
 	        return;
 	      }
	//
      }
    }

    // all tasks finished -> group waiting
    if(status_ != PROCESS_GROUP_FAIL){
      status_ = PROCESS_GROUP_WAIT;
    }
  }

  // send ready status to manager
  MASTER_EXCLUSIVE_SECTION{
    StatusType status = status_;

    if( ENABLE_FT ){
      simft::Sim_FT_MPI_Send( &status, 1, MPI_INT,  theMPISystem()->getManagerRank(), statusTag,
          theMPISystem()->getGlobalCommFT() );
    } else{
      MPI_Send(&status, 1, MPI_INT, theMPISystem()->getManagerRank(), statusTag, theMPISystem()->getGlobalComm());
    }
  }

  // reset current task
  currentTask_ = NULL;

  // if failed proc in this group detected the alive procs go into recovery state
  if( ENABLE_FT ){
    if( status_ == PROCESS_GROUP_FAIL ){
      theMPISystem()->recoverCommunicators( false );
      status_ = PROCESS_GROUP_WAIT;
    }
  }
}

void reduceSparseGridCoefficients(LevelVector& lmax,LevelVector& lmin, IndexType totalNumberOfCombis,
    IndexType currentCombi, LevelVector reduceLmin, LevelVector reduceLmax){
  for (size_t i = 0; i < reduceLmin.size(); ++i){
     if (lmin[i] > 1){
       lmin[i] -= reduceLmin[i];
     }
  }
  for (size_t i = 0; i < reduceLmax.size(); ++i){
     lmax[i] = std::max(lmin[i],lmax[i] - reduceLmax[i]);
  }
}
void ProcessGroupWorker::combineUniform() {
#ifdef DEBUG_OUTPUT

  MASTER_EXCLUSIVE_SECTION{
    std::cout << "start combining \n";
  }
#endif
  Stats::startEvent("combine init");

  // each pgrouproot must call reduce function
  //assert(tasks_.size() > 0);
  if(tasks_.size() == 0){
    std::cout << "Possible error: task size is 0! \n";
  }
  assert( combiParametersSet_ );
  int numGrids = combiParameters_.getNumGrids(); //we assume here that every task has the same number of grids

  DimType dim = combiParameters_.getDim();
  LevelVector lmin = combiParameters_.getLMin();
  LevelVector lmax = combiParameters_.getLMax();
  const std::vector<bool>& boundary = combiParameters_.getBoundary();

  // the  can be smaller than lmax because the highest subspaces do not have
  // to be exchanged
  // todo: use a flag to switch on/off optimized combination

  reduceSparseGridCoefficients(lmax,lmin,combiParameters_.getNumberOfCombinations(),currentCombi_,
      combiParameters_.getLMinReductionVector(), combiParameters_.getLMaxReductionVector());

  /*for (size_t i = 0; i < lmax.size(); ++i)
        if (lmin[i] > 1)
          lmin[i] -= 01;
  for (size_t i = 0; i < lmax.size(); ++i)
      lmax[i] = std::max(lmin[i],lmax[i] - 2);
  */
#ifdef DEBUG_OUTPUT
  MASTER_EXCLUSIVE_SECTION{
    std::cout << "lmin: "<< lmin << std::endl;
    std::cout << "lmax: "<< lmax << std::endl;
  }
#endif

  //delete old dsgs
  for(int g=0; g<combinedUniDSGVector_.size(); g++){

    if (combinedUniDSGVector_[g] != NULL)
      delete combinedUniDSGVector_[g];
  }
  combinedUniDSGVector_.clear();
  // erzeug dsgs
  combinedUniDSGVector_.resize(numGrids);
  for(int g=0; g<numGrids; g++){
    combinedUniDSGVector_[g] = new DistributedSparseGridUniform<CombiDataType>(combiScheme_, boundary,
      theMPISystem()->getLocalComm());
  }
  // todo: move to init function to avoid reregistering
  // register dsgs in all dfgs
  for (Task* t : tasks_) {
    for(int g=0; g<numGrids; g++){

      DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid(g);

      dfg.registerUniformSG(*(combinedUniDSGVector_[g]));
    }
  }
  Stats::stopEvent("combine init");
  Stats::startEvent("combine hierarchize");

  real localMax(0.0);
  //std::vector<CombiDataType> beforeCombi;
  for (Task* t : tasks_) {
    for(int g=0; g<numGrids; g++){

      DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid(g);

      // hierarchize dfg
      DistributedHierarchization::hierarchize<CombiDataType>(
          dfg, combiParameters_.getHierarchizationDims() );

      // lokales reduce auf sg ->
      dfg.addToUniformSG( *combinedUniDSGVector_[g], combiParameters_.getCoeff( t->getID() ) );
#ifdef DEBUG_OUTPUT
      std::cout << "Combination: added task " << t->getID() << " with coefficient " << combiParameters_.getCoeff( t->getID() ) <<"\n";
#endif
    }
  }
  Stats::stopEvent("combine hierarchize");

  Stats::startEvent("combine global reduce");

  for(int g=0; g<numGrids; g++){
    CombiCom::distributedGlobalReduce( *combinedUniDSGVector_[g] );
  }
  Stats::stopEvent("combine global reduce");

  //std::vector<CombiDataType> afterCombi;
  Stats::startEvent("combine dehierarchize");

  for (Task* t : tasks_) {
    for(int g=0; g<numGrids; g++){

      // get handle to dfg
      DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid(g);

      // extract dfg vom dsg
      dfg.extractFromUniformSG( *combinedUniDSGVector_[g] );

      // dehierarchize dfg
      DistributedHierarchization::dehierarchize<CombiDataType>(
          dfg, combiParameters_.getHierarchizationDims() );
    }
  }
  Stats::stopEvent("combine dehierarchize");
}


void ProcessGroupWorker::parallelEval(){
  if(uniformDecomposition)
    parallelEvalUniform();
  else
    assert( false && "not yet implemented" );
}

void ProcessGroupWorker::parallelEvalUniform(){
  assert(uniformDecomposition);

  assert(combiParametersSet_);
  int numGrids = combiParameters_.getNumGrids(); //we assume here that every task has the same number of grids

  const int dim = static_cast<int>( combiParameters_.getDim() );

  // combine must have been called before this function
  assert( combinedUniDSGVector_.size() != 0 && "you must combine before you can eval" );

  // receive leval and broadcast to group members
  std::vector<int> tmp(dim);
  MASTER_EXCLUSIVE_SECTION{
    MPI_Recv( &tmp[0], dim, MPI_INT,
              theMPISystem()->getManagerRank(), 0,
              theMPISystem()->getGlobalComm(), MPI_STATUS_IGNORE);
  }

  MPI_Bcast( &tmp[0], dim, MPI_INT,
             theMPISystem()->getMasterRank(),
             theMPISystem()->getLocalComm() );
  LevelVector leval( tmp.begin(), tmp.end() );

  // receive filename and broadcast to group members
  std::string filename;
  MASTER_EXCLUSIVE_SECTION{
    MPIUtils::receiveClass( &filename,
                            theMPISystem()->getManagerRank(),
                            theMPISystem()->getGlobalComm() );
  }

  MPIUtils::broadcastClass( &filename,
                            theMPISystem()->getMasterRank(),
                            theMPISystem()->getLocalComm() );

  for(int g=0; g < numGrids; g++){//loop over all grids and plot them
    // create dfg
    bool forwardDecomposition = !isGENE;
    DistributedFullGrid<CombiDataType> dfg( dim, leval,
                                            combiParameters_.getApplicationComm(),
                                            combiParameters_.getBoundary(),
                                            combiParameters_.getParallelization(),
                                            forwardDecomposition
                                            );

    // register dsg
    dfg.registerUniformSG(*combinedUniDSGVector_[g]);

    // fill dfg with hierarchical coefficients from distributed sparse grid
    dfg.extractFromUniformSG( *combinedUniDSGVector_[g] );

    // dehierarchize dfg
    DistributedHierarchization::dehierarchize<CombiDataType>(
        dfg, combiParameters_.getHierarchizationDims() );
    std::string fn = filename;
    fn = fn + std::to_string(g);
    // save dfg to file with MPI-IO
    dfg.writePlotFile( fn.c_str() );
  }
}


void ProcessGroupWorker::gridEval() { //not supported anymore
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
  MASTER_EXCLUSIVE_SECTION{
    // receive size of levelvector = dimensionality
    MPI_Status status;
    int bsize;
    MPI_Probe( theMPISystem()->getManagerRank(), 0, theMPISystem()->getGlobalComm(), &status);
    MPI_Get_count(&status, MPI_INT, &bsize);

    assert(bsize == static_cast<int>(dim));

    std::vector<int> tmp(dim);
    MPI_Recv( &tmp[0], bsize, MPI_INT,
              theMPISystem()->getManagerRank(), 0,
              theMPISystem()->getGlobalComm(), MPI_STATUS_IGNORE);
    leval = LevelVector(tmp.begin(), tmp.end());
  }

  assert( combiParametersSet_ );
  const std::vector<bool>& boundary = combiParameters_.getBoundary();
  FullGrid<CombiDataType> fg_red(dim, leval, boundary);

  // create the empty grid on only on localroot
  MASTER_EXCLUSIVE_SECTION {
    fg_red.createFullGrid();
  }

  // collect fg on pgrouproot and reduce
  for (size_t i = 0; i < tasks_.size(); ++i) {
    Task* t = tasks_[i];

    FullGrid<CombiDataType> fg(t->getDim(), t->getLevelVector(), boundary );

    MASTER_EXCLUSIVE_SECTION {
      fg.createFullGrid();
    }

    t->getFullGrid( fg,
                    theMPISystem()->getMasterRank(),
                    theMPISystem()->getLocalComm() );

    MASTER_EXCLUSIVE_SECTION{
      fg_red.add(fg, combiParameters_.getCoeff( t->getID() ) );
    }
  }
  // global reduce of f_red
  MASTER_EXCLUSIVE_SECTION {
    CombiCom::FGReduce( fg_red,
                        theMPISystem()->getManagerRank(),
                        theMPISystem()->getGlobalComm() );
  }
}

//todo: this is just a temporary function which will drop out some day
// also this function requires a modified fgreduce method which uses allreduce
// instead reduce in manger
void ProcessGroupWorker::combineFG() {
  //gridEval();

  // TODO: Sync back to fullgrids
}

void ProcessGroupWorker::updateCombiParameters() {
  CombiParameters tmp;

  // local root receives task
  MASTER_EXCLUSIVE_SECTION {
    MPIUtils::receiveClass(
        &tmp,
        theMPISystem()->getManagerRank(),
        theMPISystem()->getGlobalComm() );
    std::cout << "master received combiparameters \n";
	std::cout << "master size: " << tmp.getLevelsToIDs().size() << "\n";
  }

  // broadcast task to other process of pgroup
  MPIUtils::broadcastClass(
      &tmp,
      theMPISystem()->getMasterRank(),
      theMPISystem()->getLocalComm() );
  std::cout << "worker received combiparameters \n";
  if(combiParameters_.isApplicationCommSet()){
    CommunicatorType free = combiParameters_.getApplicationComm();
    if(free != NULL && free != MPI_COMM_NULL){
      MPI_Comm_free(&free);
    }
  }

  combiParameters_ = tmp;
  if(!combiParametersSet_){
	  combiScheme_ = DimAdaptiveCombiScheme {combiParameters_.getLMin(), combiParameters_.getLMax()};
  }

  combiParametersSet_ = true;

}


void ProcessGroupWorker::setCombinedSolutionUniform( Task* t ) {
  assert( combinedUniDSGVector_.size() != 0 );
  assert(combiParametersSet_);

  int numGrids = combiParameters_.getNumGrids(); //we assume here that every task has the same number of grids

  for(int g=0; g < numGrids;g++){
    assert( combinedUniDSGVector_[g] != NULL );

    // get handle to dfg
    DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid(g);

    // extract dfg vom dsg
    dfg.extractFromUniformSG( *combinedUniDSGVector_[g] );

    // dehierarchize dfg
    DistributedHierarchization::dehierarchize<CombiDataType>(
        dfg, combiParameters_.getHierarchizationDims() );
  }
}

int ProcessGroupWorker::getProcTask(int taskID){
	return taskToProc_.at(taskID);
}

void ProcessGroupWorker::findBestExpansion(){
	int messageTag = 1;
	double bestError = -1;
	LevelVector bestExpansion {};

	for(const auto& activeNode : combiScheme_.getActiveNodes()){
		for(DimType i = 0; i < combiScheme_.dim(); ++i){
			LevelVector potExpansion = activeNode;
			++potExpansion.at(i);
			if(!combiScheme_.isExpansion(potExpansion) || combiScheme_.isBorder(activeNode, i)){
				continue;
			}

			const auto cmpPair = combiScheme_.getPosNegPair(activeNode, i);
			std::cout << "active: " << cmpPair.first << "\n";
			std::cout << "bwd: " << cmpPair.second << "\n";
			std::cout << "size: " << combiParameters_.getLevelsToIDs().size() << "\n";

			for(auto test : combiParameters_.getLevelsToIDs()){
				std::cout << test.first << " " << test.second << "\n";
			}
			const int activeTaskID = combiParameters_.getID(cmpPair.first);
			const int bwdTaskID = combiParameters_.getID(cmpPair.second);


			//This calculation depends on the current implementation of then
			//rank calculations in MPI System so it might break easily
			const int rankOffset = theMPISystem()->getLocalRank();

			const int activeNodeRank = theMPISystem()->getWorldRankFromGlobalRank(getProcTask(activeTaskID)) + rankOffset;
			const int bwdNeighRank = theMPISystem()->getWorldRankFromGlobalRank(getProcTask(bwdTaskID)) + rankOffset;
			const bool activeNodeOwned = activeNodeRank == theMPISystem()->getWorldRank();
			const bool bwdNeighOwned = bwdNeighRank == theMPISystem()->getWorldRank();



			if(activeNodeOwned){
				std::vector<CombiDataType> activeSubGrid {};
				std::vector<CombiDataType> bwdSubGrid {};
				if(bwdNeighOwned){ //The proc itself owns both tasks
					Task *activeNodeTask = getTask(activeTaskID);
					Task *bwdNeighTask = getTask(bwdTaskID);
					activeSubGrid = activeNodeTask->getDistributedFullGrid().getSubGrid(0);
					bwdSubGrid= bwdNeighTask->getDistributedFullGrid().getSubGrid(0);
					assert(activeSubGrid.size() == bwdSubGrid.size());
				} else {
					Task *activeNodeTask = getTask(activeTaskID);
					activeSubGrid = activeNodeTask->getDistributedFullGrid().getSubGrid(0);
					bwdSubGrid.resize(activeSubGrid.size());
					MPI_Recv(bwdSubGrid.data(), bwdSubGrid.size(), MPI_DOUBLE, bwdNeighRank, messageTag, theMPISystem()->getWorldComm(), MPI_STATUS_IGNORE);
				}

				double error = 0;
				if(!activeSubGrid.empty()){
					const double expansionGridPoints = std::accumulate(std::begin(activeNode), std::end(activeNode), 1, [](double acc, LevelType val) {
						return acc * (1 << (val - 1));
					});

					error = maxRelativeError(activeSubGrid, bwdSubGrid);
					//error /= expansionGridPoints;
				}
				MASTER_EXCLUSIVE_SECTION{
					MPI_Reduce(MPI_IN_PLACE, &error, 1, MPI_DOUBLE, MPI_SUM, 0, theMPISystem()->getLocalComm());
					if(error > bestError){
						bestError = error;
						bestExpansion = activeNode;
					}
				} else {
					MPI_Reduce(&error, &error, 1, MPI_DOUBLE, MPI_SUM, 0, theMPISystem()->getLocalComm());
				}

			} else if (bwdNeighOwned){ //the proc owns the bwdNeighbour
				//since the bwd neighbour is smaller it is sent to the proc
				//with the active node
				Task *bwdNeighTask = getTask(bwdTaskID);
				auto bwdSubGrid= bwdNeighTask->getDistributedFullGrid().getSubGrid(0);
				MPI_Send(bwdSubGrid.data(), static_cast<int>(bwdSubGrid.size()), MPI_DOUBLE, activeNodeRank, messageTag, theMPISystem()->getWorldComm());
			}
			//if no if-branch was executed the proc owns nothing so we can simply continue
			++messageTag;
		}
	}

	MASTER_EXCLUSIVE_SECTION{
		MPI_Send(&bestError, 1, MPI_DOUBLE, theMPISystem()->getManagerRank(), 1234, theMPISystem()->getGlobalComm());
		MPI_Send(bestExpansion.data(), bestExpansion.size(), MPI_LONG, theMPISystem()->getManagerRank(), 1235, theMPISystem()->getGlobalComm());
		std::cout << "Sent best expansion\n";
	}
}

} /* namespace combigrid */
