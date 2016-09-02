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
#include "sgpp/distributedcombigrid/sparsegrid/DistributedSparseGrid.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/DistributedSparseGridUniform.hpp"
#include "sgpp/distributedcombigrid/combicom/CombiCom.hpp"
#include "sgpp/distributedcombigrid/hierarchization/DistributedHierarchization.hpp"
#include "sgpp/distributedcombigrid/mpi/MPIUtils.hpp"
#include <set>

namespace combigrid {

ProcessGroupWorker::ProcessGroupWorker() :
        currentTask_( NULL),
        status_(PROCESS_GROUP_WAIT),
        combinedFG_( NULL),
        combinedUniDSG_(NULL),
        combinedFGexists_(false),
        combiParameters_(),
        combiParametersSet_(false)
{
}

ProcessGroupWorker::~ProcessGroupWorker() {
  delete combinedFG_;
}

SignalType ProcessGroupWorker::wait() {
  if (status_ == PROCESS_GROUP_BUSY)
    return RUN_NEXT;

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
  // process signal
  if (signal == RUN_FIRST) {

    Task* t;

    // local root receives task
    MASTER_EXCLUSIVE_SECTION {
      Task::receive(  &t,
          theMPISystem()->getManagerRank(),
          theMPISystem()->getGlobalComm() );
    }

    // broadcast task to other process of pgroup
    Task::broadcast(&t, theMPISystem()->getMasterRank(), theMPISystem()->getLocalComm());

    MPI_Barrier(theMPISystem()->getLocalComm());

    // add task to task storage
    tasks_.push_back(t);

    status_ = PROCESS_GROUP_BUSY;

    // set currentTask
    currentTask_ = tasks_.back();

    // initalize task
    currentTask_->init(theMPISystem()->getLocalComm());

    // execute task
    currentTask_->run(theMPISystem()->getLocalComm());

  } else if (signal == RUN_NEXT) {
    // this should not happen
    assert(tasks_.size() > 0);

    // reset finished status of all tasks
    for (size_t i = 0; i < tasks_.size(); ++i)
      tasks_[i]->setFinished(false);

    status_ = PROCESS_GROUP_BUSY;

    // set currentTask
    currentTask_ = tasks_[0];

    // run first task
    currentTask_->run(theMPISystem()->getLocalComm());

  } else if (signal == ADD_TASK) {
    std::cout << "adding a single task" << std::endl;

    Task* t;

    // local root receives task
    MASTER_EXCLUSIVE_SECTION {
      Task::receive(&t, theMPISystem()->getManagerRank(), theMPISystem()->getGlobalComm());
    }

    // broadcast task to other process of pgroup
    Task::broadcast(&t, theMPISystem()->getMasterRank(), theMPISystem()->getLocalComm());

    MPI_Barrier(theMPISystem()->getLocalComm());

    // check if task already exists on this group
    for ( auto tmp : tasks_ )
      assert( tmp->getID() != t->getID() );

    // initalize task and set values to zero
    // the task will get the proper initial solution during the next combine
    t->init( theMPISystem()->getLocalComm() );

    t->setZero();

    t->setFinished( true );

    // add task to task storage
    tasks_.push_back(t);

    status_ = PROCESS_GROUP_BUSY;

  } else if (signal == EVAL) {
    // receive x

    // loop over all tasks
    // t.eval(x)
  } else if (signal == EXIT) {

  } else if (signal == SYNC_TASKS) {
    MASTER_EXCLUSIVE_SECTION {
      for (size_t i = 0; i < tasks_.size(); ++i) {
        Task::send(&tasks_[i], theMPISystem()->getManagerRank(), theMPISystem()->getGlobalComm());
      }
    }
  } else if (signal == COMBINE) {

    combineUniform();

  } else if (signal == GRID_EVAL) {

    gridEval();
    return signal;

  } else if (signal == COMBINE_FG) {

    combineFG();

  } else if (signal == UPDATE_COMBI_PARAMETERS) {

    updateCombiParameters();

  } else if (signal == RECOMPUTE) {
    Task* t;

    // local root receives task
    MASTER_EXCLUSIVE_SECTION {
      Task::receive(&t, theMPISystem()->getManagerRank(), theMPISystem()->getGlobalComm());
    }

    // broadcast task to other process of pgroup
    Task::broadcast(&t, theMPISystem()->getMasterRank(), theMPISystem()->getLocalComm());

    MPI_Barrier(theMPISystem()->getLocalComm());

    // add task to task storage
    tasks_.push_back(t);

    status_ = PROCESS_GROUP_BUSY;

    // set currentTask
    currentTask_ = tasks_.back();

    // initalize task
    currentTask_->init(theMPISystem()->getLocalComm());

    currentTask_->setZero();

    // fill task with combisolution
    setCombinedSolutionUniform( currentTask_ );

    // execute task
    currentTask_->run(theMPISystem()->getLocalComm());
    //  } else if ( signal ==  RECOVER_COMM ){
    //    theMPISystem()->recoverCommunicators( true );
    //    return signal;
  } else if (signal == SEARCH_SDC) {
      searchSDC();
  }

  // in the general case: send ready signal.
  if(!omitReadySignal)
    ready();

  return signal;
}

void ProcessGroupWorker::ready() {

  // check if there are unfinished tasks
  for (size_t i = 0; i < tasks_.size(); ++i) {
    if (!tasks_[i]->isFinished()) {
      status_ = PROCESS_GROUP_BUSY;

      // set currentTask
      currentTask_ = tasks_[i];
      currentTask_->run(theMPISystem()->getLocalComm());
    }
  }

  // all tasks finished -> group waiting
  if( status_ != PROCESS_GROUP_FAIL )
    status_ = PROCESS_GROUP_WAIT;

  MASTER_EXCLUSIVE_SECTION{
    StatusType status = status_;
    MPI_Send(&status, 1, MPI_INT, theMPISystem()->getManagerRank(), statusTag, theMPISystem()->getGlobalComm());
  }

  // reset current task
  currentTask_ = NULL;
}

void ProcessGroupWorker::combine() {
  // early exit if no tasks available
  // todo: doesnt work, each pgrouproot must call reduce function
  assert(tasks_.size() > 0);

  assert( combiParametersSet_ );
  DimType dim = combiParameters_.getDim();
  const LevelVector& lmin = combiParameters_.getLMin();
  const LevelVector& lmax = combiParameters_.getLMax();
  const std::vector<bool>& boundary = combiParameters_.getBoundary();

  // erzeug dsg
  DistributedSparseGrid<CombiDataType> dsg(dim, lmax, lmin, boundary, theMPISystem()->getLocalComm());

  for (Task* t : tasks_) {
    DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid();

    // hierarchize dfg
    DistributedHierarchization::hierarchize<CombiDataType>(dfg);

    // lokales reduce auf sg ->
    //CombiCom::distributedLocalReduce<CombiDataType>( dfg, dsg, combiParameters_.getCoeff( t->getID() ) );
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
}

void ProcessGroupWorker::combineUniform() {
  // early exit if no tasks available
  // todo: doesnt work, each pgrouproot must call reduce function
  assert(tasks_.size() > 0);

  assert( combiParametersSet_ );
  DimType dim = combiParameters_.getDim();
  const LevelVector& lmin = combiParameters_.getLMin();
  LevelVector lmax = combiParameters_.getLMax();
  const std::vector<bool>& boundary = combiParameters_.getBoundary();;

  for (size_t i = 0; i < lmax.size(); ++i)
    if (lmax[i] > 1)
      lmax[i] -= 1;

  // todo: delete old dsg
  if (combinedUniDSG_ != NULL)
    delete combinedUniDSG_;

  // erzeug dsg
  combinedUniDSG_ = new DistributedSparseGridUniform<CombiDataType>(dim, lmax,
      lmin, boundary,
      theMPISystem()->getLocalComm());

  // todo: move to init function to avoid reregistering
  // register dsg in all dfgs
  for (Task* t : tasks_) {
    DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid();

    dfg.registerUniformSG(*combinedUniDSG_);
  }

  for (Task* t : tasks_) {

    DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid();

    // hierarchize dfg
    DistributedHierarchization::hierarchize<CombiDataType>(dfg);

    // lokales reduce auf sg ->
    dfg.addToUniformSG( *combinedUniDSG_, combiParameters_.getCoeff( t->getID() ) );
  }


  CombiCom::distributedGlobalReduce( *combinedUniDSG_ );

  for (Task* t : tasks_) {
    // get handle to dfg
    DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid();

    // extract dfg vom dsg
    dfg.extractFromUniformSG( *combinedUniDSG_ );

    // dehierarchize dfg
    DistributedHierarchization::dehierarchize<CombiDataType>( dfg );
  }

}

void ProcessGroupWorker::gridEval() {
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
  }

  // broadcast task to other process of pgroup
  MPIUtils::broadcastClass(
      &tmp,
      theMPISystem()->getMasterRank(),
      theMPISystem()->getLocalComm() );

  combiParameters_ = tmp;

  combiParametersSet_ = true;

  status_ = PROCESS_GROUP_BUSY;

}


void ProcessGroupWorker::setCombinedSolutionUniform( Task* t ) {
  assert( combinedUniDSG_ != NULL );

  // get handle to dfg
  DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid();

  // extract dfg vom dsg
  dfg.extractFromUniformSG( *combinedUniDSG_ );

  // dehierarchize dfg
  DistributedHierarchization::dehierarchize<CombiDataType>( dfg );
}

void ProcessGroupWorker::searchSDC(){
  SDCMethodType method = -1;
  MASTER_EXCLUSIVE_SECTION {
    // receive signal from manager
    MPI_Recv( &method, 1, MPI_INT,
        theMPISystem()->getManagerRank(),
        infoTag,
        theMPISystem()->getGlobalComm(),
        MPI_STATUS_IGNORE);
  }
  // distribute signal to other processes of pgroup
  MPI_Bcast( &method, 1, MPI_INT,
      theMPISystem()->getMasterRank(),
      theMPISystem()->getLocalComm() );

  int success;
  if( method == COMPARE_PAIRS )
    success = comparePairs( 5 );
  if( method == COMPARE_VALUES )
    success = compareValues();

  std::vector<int> successes;
  successes.resize(theMPISystem()->getNumProcs());

  MPI_Gather( &success, 1, MPI_INT,
      &successes[0], 1, MPI_INT,
      theMPISystem()->getMasterRank(), theMPISystem()->getLocalComm() );

  MASTER_EXCLUSIVE_SECTION{
    int failedWorker = std::find(successes.begin(), successes.end(), 0)-successes.begin();
    if ( failedWorker < successes.size() ){
      status_ = PROCESS_GROUP_FAIL;
      std::vector<int> levelsSDC;
      levelsSDC.resize(tasks_.size());
      MPI_Status statusSDC;
      MPI_Recv( &levelsSDC[0], tasks_.size(), MPI_INT, failedWorker, infoTag, theMPISystem()->getLocalComm(), &statusSDC );
      int numTasksSDC;
      MPI_Get_count( &statusSDC, MPI_INT, &numTasksSDC );
      levelsSDC.resize( numTasksSDC );
      MPI_Send( &levelsSDC[0], numTasksSDC, MPI_INT, theMPISystem()->getManagerRank(), infoTag, theMPISystem()->getGlobalComm());
    }
  }
}

int ProcessGroupWorker::comparePairs( int numNearestNeighbors ){

  /* Generate all pairs of grids */
  std::vector<std::vector<Task*>> allPairs;

  generatePairs( numNearestNeighbors, allPairs );

  std::vector<CombiDataType> allBetas;
  std::vector<CombiDataType> allBetasSum;
  std::vector<LevelVector> allSubs;
  std::vector<size_t> allJs;

  MPI_File_open(theMPISystem()->getLocalComm(), "out/all-betas-0.txt", MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &betasFile_ );

  for (auto pair : allPairs){

    DistributedFullGrid<CombiDataType>& dfg_t = pair[0]->getDistributedFullGrid();
    DistributedFullGrid<CombiDataType>& dfg_s = pair[1]->getDistributedFullGrid();

    LevelVector t_level = pair[0]->getLevelVector();
    LevelVector s_level = pair[1]->getLevelVector();

    DistributedSparseGridUniform<CombiDataType>* SDCUniDSG = new DistributedSparseGridUniform<CombiDataType>(
        combiParameters_.getDim(), combiParameters_.getLMax(), combiParameters_.getLMin(),
        combiParameters_.getBoundary(), theMPISystem()->getLocalComm());

    dfg_t.registerUniformSG( *SDCUniDSG );
    dfg_s.registerUniformSG( *SDCUniDSG );

    dfg_t.addToUniformSG( *SDCUniDSG, 1.0 );
    dfg_s.addToUniformSG( *SDCUniDSG, -1.0, false );

    LevelVector common_level;
    for (size_t i = 0; i < t_level.size(); ++i)
      common_level.push_back( (t_level[i] <= s_level[i]) ? t_level[i] : s_level[i] );

    CombiDataType localBetaMax(0.0);

    LevelVector subMax;
    size_t jMax = 0;

    for (size_t i = 0; i < SDCUniDSG->getNumSubspaces(); ++i){
      if (SDCUniDSG->getLevelVector(i) <= common_level){
        auto subData = SDCUniDSG->getData(i);
        auto subSize = SDCUniDSG->getDataSize(i);
        for (size_t j = 0; j < subSize; ++j){
          if (std::abs(subData[j]) > std::abs(localBetaMax)){
            localBetaMax = subData[j];
            subMax = SDCUniDSG->getLevelVector(i);
            jMax = j;
          }
        }
      }
    }
    allBetas.push_back( localBetaMax );
    allSubs.push_back( subMax );
    allJs.push_back( jMax );
  }

  std::vector<CombiDataType> allBetasReduced;
  allBetasReduced.resize(allBetas.size());

  CombiCom::BetasReduce( allBetas, allBetasReduced, theMPISystem()->getLocalComm() );

  auto globalBetaMax = std::max_element(allBetasReduced.begin(), allBetasReduced.end(),
      [](CombiDataType a, CombiDataType b){ return std::abs(a) < std::abs(b); } );

  auto b = std::find( allBetas.begin(), allBetas.end(), *globalBetaMax );

  betas_.clear();

  if(b != allBetas.end()) {

    size_t indMax = std::distance(allBetas.begin(), b);

    LevelVector subMax = allSubs[indMax];
    size_t jMax = allJs[indMax];

    for (auto pair : allPairs){

      DistributedFullGrid<CombiDataType>& dfg_t = pair[0]->getDistributedFullGrid();
      DistributedFullGrid<CombiDataType>& dfg_s = pair[1]->getDistributedFullGrid();

      LevelVector t_level = pair[0]->getLevelVector();
      LevelVector s_level = pair[1]->getLevelVector();

      DistributedSparseGridUniform<CombiDataType>* SDCUniDSG = new DistributedSparseGridUniform<CombiDataType>(
          combiParameters_.getDim(), combiParameters_.getLMax(), combiParameters_.getLMin(),
          combiParameters_.getBoundary(), theMPISystem()->getLocalComm());

      dfg_t.registerUniformSG( *SDCUniDSG );
      dfg_s.registerUniformSG( *SDCUniDSG );

      dfg_t.addToUniformSG( *SDCUniDSG, 1.0 );
      dfg_s.addToUniformSG( *SDCUniDSG, -1.0, false );

      LevelVector common_level;
      for (size_t i = 0; i < t_level.size(); ++i)
        common_level.push_back( (t_level[i] <= s_level[i]) ? t_level[i] : s_level[i] );

      if (subMax <= common_level){
        auto subData = SDCUniDSG->getData(subMax);
        CombiDataType localBetaMax = subData[jMax];
        betas_[std::make_pair(t_level, s_level)] = localBetaMax;
      }
      else
        betas_[std::make_pair(t_level, s_level)] = 0;
    }

    std::vector<int> levelsSDC;
    filterSDC( levelsSDC );

    MPI_File_close( &betasFile_ );

    if( levelsSDC.size() != 0 ){
      //      status_ = PROCESS_GROUP_FAIL;
      // send SDC grids to master
      MPI_Request req;
      MPI_Isend( &levelsSDC[0], levelsSDC.size(), MPI_INT, theMPISystem()->getMasterRank(), infoTag, theMPISystem()->getLocalComm(), &req );
      return 0;
    }
  }
  return 1;
}

int ProcessGroupWorker::compareValues(){
}

void ProcessGroupWorker::generatePairs( int numNearestNeighbors, std::vector<std::vector<Task*>> &allPairs ){

  std::vector<LevelVector> levels;
  std::map<LevelVector, int> numPairs;

  for ( auto tt: tasks_ ){
    levels.push_back(tt->getLevelVector());
    numPairs[tt->getLevelVector()] = 0;
  }

  for (Task* s : tasks_ ){

    std::sort(levels.begin(), levels.end(), [s](LevelVector const& a, LevelVector const& b) {
      return l1(a - s->getLevelVector()) < l1(b - s->getLevelVector());
    });

    int k = 0;

    for( size_t t_i = 1; t_i < levels.size(); ++t_i ){
      std::vector<Task*> currentPair;

      Task* t = *std::find_if(tasks_.begin(), tasks_.end(),
          [levels,t_i](Task* const &tt) -> bool { return tt->getLevelVector() == levels[t_i]; });

      currentPair.push_back(t);
      currentPair.push_back(s);

      if(std::find(allPairs.begin(), allPairs.end(), currentPair) == allPairs.end()
          && numPairs[s->getLevelVector()] < numNearestNeighbors
          && numPairs[t->getLevelVector()] < numNearestNeighbors ){
        allPairs.push_back({currentPair[1],currentPair[0]});
        numPairs[s->getLevelVector()]++;
        numPairs[t->getLevelVector()]++;
        k++;
      }

      if (k == numNearestNeighbors)
        break;
    }
  }
}

void ProcessGroupWorker::filterSDC( std::vector<int> &levelsSDC ){

  // Number of measurements (beta values)
  size_t n = betas_.size();

  auto lmin = combiParameters_.getLMin();
  auto lmax = combiParameters_.getLMax();
  size_t diff = lmax[0] - lmin[0] + 1;

  // Number of unknowns (functions D1, D2, and D12)
  //  size_t p = 2*diff + diff*(diff+1)/2;
  size_t p = 2*diff;

  std::cout<<"n = "<< n <<", p = "<<p<<std::endl;
  if ( n < p )
    return;

  gsl_multifit_robust_workspace *regressionWsp = gsl_multifit_robust_alloc (gsl_multifit_robust_cauchy, n , p );

  gsl_multifit_robust_tune( .1, regressionWsp );

  gsl_multifit_robust_maxiter( 100, regressionWsp );

  gsl_matrix *X = gsl_matrix_alloc( n, p );
  gsl_vector *y = gsl_vector_alloc( n );
  gsl_vector *c = gsl_vector_alloc( p );
  gsl_vector *r = gsl_vector_alloc( n );
  gsl_matrix *cov = gsl_matrix_alloc( p, p );

  // Initialize matrix with zeros
  gsl_matrix_set_zero( X );

  IndexType idx_D1, idx_D2, idx_D12;

  // Helper function to compute indices of D12
  auto idx = [](IndexType d, IndexType i){
    std::vector<IndexType> numbers(i);
    std::iota(numbers.begin(), numbers.end(),d-i+1);
    return std::accumulate(numbers.begin(), numbers.end(), 0);
  };

  int row  = 0;
  CombiDataType val;
  for( auto const &entry : betas_ ){
    LevelVector key_t = entry.first.first;
    LevelVector key_s = entry.first.second;
    CombiDataType beta = entry.second;

    std::vector<CombiDataType> ht = {1.0/pow(2.0,key_t[0]), 1.0/pow(2.0,key_t[1])};
    std::vector<CombiDataType> hs = {1.0/pow(2.0,key_s[0]), 1.0/pow(2.0,key_s[1])};
    for( size_t i = 0; i < key_t.size(); ++i ){
      key_t[i] -= lmin[0];
      key_s[i] -= lmin[0];
    }

    idx_D1 =  key_t[0];
    gsl_matrix_set( X, row, idx_D1, pow(ht[0],2) );

    idx_D1 =  key_s[0];
    val = gsl_matrix_get( X, row, idx_D1 );
    gsl_matrix_set( X, row, idx_D1, val - pow(hs[0],2) );

    idx_D2 = diff + key_t[1];
    gsl_matrix_set( X, row, idx_D2, pow(ht[1],2) );

    idx_D2 = diff + key_s[1];
    val = gsl_matrix_get( X, row, idx_D2 );
    gsl_matrix_set( X, row, idx_D2, val - pow(hs[1],2) );

    //    idx_D12 = 2*diff + idx(diff,key_t[1]) + key_t[0];
    //    gsl_matrix_set( X, row, idx_D12, pow(ht[0]*ht[1],2) );
    //
    //    idx_D12 = 2*diff + idx(diff,key_s[1]) + key_s[0];
    //    val = gsl_matrix_get( X, row, idx_D12 );
    //    gsl_matrix_set( X, row, idx_D12, val - pow(hs[0]*hs[1],2) );

    gsl_vector_set( y, row, beta );

    row++;
  }

  gsl_set_error_handler_off();
  gsl_multifit_robust( X, y, c, cov, regressionWsp );
  gsl_multifit_robust_residuals( X, y, c, r, regressionWsp );

  // Write pairs and their beta values to file
  std::stringstream buf;
  buf<<n<< std::endl;
  row = 0;
  for ( auto const &entry : betas_ ){
    LevelVector key_t = entry.first.first;
    LevelVector key_s = entry.first.second;
    CombiDataType beta = entry.second;
    CombiDataType res = r->data[row];
    buf << key_t <<","<< key_s <<","<< beta << ", "<< res <<std::endl;
    MPI_File_seek(betasFile_, 0, MPI_SEEK_END);
    MPI_File_write(betasFile_, buf.str().c_str(), buf.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    std::cout<<buf.str();
    buf.str("");
    row++;
  }
  // Write regression coefficients to file
  buf << c->size << std::endl;
  for( size_t i = 0; i < c->size; ++i){
    buf << c->data[i]<<std::endl;
    MPI_File_seek(betasFile_, 0, MPI_SEEK_END);
    MPI_File_write(betasFile_, buf.str().c_str(), buf.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    buf.str("");
  }

  std::map<LevelVector,int> mapSDC;
  for ( auto const &entry : betas_ ){
    LevelVector key_t = entry.first.first;
    LevelVector key_s = entry.first.second;
    mapSDC[key_t] = 0;
    mapSDC[key_s] = 0;
  }

  row = 0;
  for ( auto const &entry : betas_ ){
    LevelVector key_t = entry.first.first;
    LevelVector key_s = entry.first.second;
    CombiDataType beta = entry.second;
    CombiDataType res = r->data[row];
    if ( std::abs(res) > 2.5 && beta != 0 ){
      mapSDC[key_t]++;
      mapSDC[key_s]++;
    }
    row++;
  }
  std::cout<< "SDC grid: " << std::endl;
  for (auto s : mapSDC){
    if (s.second >= 2){
      int id = combiParameters_.getID(s.first);
      levelsSDC.push_back(id);
    }
  }

  gsl_matrix_free(X);
  gsl_matrix_free(cov);
  gsl_vector_free(y);
  gsl_vector_free(c);
  gsl_vector_free(r);
  gsl_multifit_robust_free(regressionWsp);
}

//  void addToUniformSG(DistributedSparseGridUniform<FG_ELEMENT>& dsg,
//                      real coeff) {
//    // test if dsg has already been registered
//    if (&dsg != dsg_)
//      registerUniformSG(dsg);
//
//    // create iterator for each subspace in dfg
//    typedef typename std::vector<FG_ELEMENT>::iterator SubspaceIterator;
//    typename std::vector<SubspaceIterator> it_sub(
//      subspaceAssigmentList_.size());
//
//    for (size_t subFgId = 0; subFgId < it_sub.size(); ++subFgId) {
//      if (subspaceAssigmentList_[subFgId] < 0)
//        continue;
//
//      IndexType subSgId = subspaceAssigmentList_[subFgId];
//
//      it_sub[subFgId] = dsg.getDataVector(subSgId).begin();
//    }
//
//    // loop over all grid points
//    for (size_t i = 0; i < fullgridVector_.size(); ++i) {
//      // get subspace_fg id
//      size_t subFgId(assigmentList_[i]);
//
//      if (subspaceAssigmentList_[subFgId] < 0)
//        continue;
//
//      IndexType subSgId = subspaceAssigmentList_[subFgId];
//
//      assert(it_sub[subFgId] != dsg.getDataVector(subSgId).end());
//
//      // add grid point to subspace, mul with coeff
//      *it_sub[subFgId] += coeff * fullgridVector_[i];
//
//      ++it_sub[subFgId];
//    }
//  }
} /* namespace combigrid */
