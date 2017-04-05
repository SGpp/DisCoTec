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
#include "sgpp/distributedcombigrid/hierarchization/Hierarchization.hpp"
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
  } else if (signal == REINIT_TASK) {
    std::cout << "reinitializing a single task" << std::endl;

    Task* t;

    // local root receives task
    MASTER_EXCLUSIVE_SECTION {
      Task::receive(&t, theMPISystem()->getManagerRank(), theMPISystem()->getGlobalComm());
    }

    // broadcast task to other process of pgroup
    Task::broadcast(&t, theMPISystem()->getMasterRank(), theMPISystem()->getLocalComm());

    MPI_Barrier(theMPISystem()->getLocalComm());

    for (auto tt : tasks_){
      if (tt->getID() == t->getID()){
        currentTask_ = tt;
        break;
      }
    }

    // initalize task and set values to zero
    currentTask_->init( theMPISystem()->getLocalComm() );
    currentTask_->setZero();
    setCombinedSolutionUniform( currentTask_ );
    currentTask_->setFinished(true);
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
//  if (combinedUniDSG_ != NULL)
//    delete combinedUniDSG_;

  // erzeug dsg
//  combinedUniDSG_ = new DistributedSparseGridUniform<CombiDataType>(dim, lmax,
//      lmin, boundary,
//      theMPISystem()->getLocalComm());

  // todo: move to init function to avoid reregistering
  // register dsg in all dfgs
//  for (Task* t : tasks_) {
//    DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid();
//
//    dfg.registerUniformSG(*combinedUniDSG_);
//  }

  for (Task* t : tasks_) {

    DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid();

    // hierarchize dfg
//    DistributedHierarchization::hierarchize<CombiDataType>(dfg);

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

  // extract dfg from dsg
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

  std::vector<int> levelsSDC;
  compareSolutions( static_cast<int>(combiParameters_.getDim()), levelsSDC, method );

  int numLocalSDC = static_cast<int>(levelsSDC.size());
  int numGlobalSDC;
  MPI_Allreduce( &numLocalSDC, &numGlobalSDC, 1, MPI_INT, MPI_SUM, theMPISystem()->getLocalComm());
  if ( numGlobalSDC > 0 ){
    MPI_Status statusSDC;
    if( numLocalSDC > 0 ){
      if (!theMPISystem()->isMaster()){
        MPI_Send( &levelsSDC[0], numLocalSDC, MPI_INT, theMPISystem()->getMasterRank(), infoTag, theMPISystem()->getLocalComm() );
      }
    } else {
      MASTER_EXCLUSIVE_SECTION{
        levelsSDC.resize(tasks_.size());
        MPI_Recv( &levelsSDC[0], static_cast<int>(tasks_.size()), MPI_INT, MPI_ANY_SOURCE, infoTag, theMPISystem()->getLocalComm(), &statusSDC );
        MPI_Get_count( &statusSDC, MPI_INT, &numGlobalSDC );
        levelsSDC.resize(numGlobalSDC);
      }
    }
    MASTER_EXCLUSIVE_SECTION{
      status_ = PROCESS_GROUP_FAIL;
      MPI_Send( &levelsSDC[0], numGlobalSDC, MPI_INT, theMPISystem()->getManagerRank(), infoTag, theMPISystem()->getGlobalComm());
    }
  }
}

void ProcessGroupWorker::compareSolutions( int numNearestNeighbors, std::vector<int> &levelsSDC, SDCMethodType method ){

  // todo: delete old dsg
  if (combinedUniDSG_ != NULL)
    delete combinedUniDSG_;

  // erzeug dsg
  combinedUniDSG_ = new DistributedSparseGridUniform<CombiDataType>(
      combiParameters_.getDim(), combiParameters_.getLMax(), combiParameters_.getLMin(),
      combiParameters_.getBoundary(), theMPISystem()->getLocalComm());

  // We use this DSG to compare component solutions and search for SDC
  DistributedSparseGridUniform<CombiDataType>* SDCUniDSG = combinedUniDSG_;

  // We work in the hierarchical basis, so hierarchization is now done here
  for (auto t : tasks_){
    DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid();
    DistributedHierarchization::hierarchize<CombiDataType>(dfg);
    dfg.registerUniformSG( *SDCUniDSG );
  }

  if ( tasks_.size() == 1 )
    return;

  // Generate all pairs of grids
  std::vector<std::vector<Task*>> allPairs;
  generatePairs( numNearestNeighbors, allPairs );

  std::vector<CombiDataType> allBetas;
  std::vector<CombiDataType> allBetasSum;
  std::vector<LevelVector> allSubs;
  std::vector<size_t> allJs;

  // Iterate throuh all pairs of grids. For each pair, we store the gridpoint with
  // the largest difference
  for (auto pair : allPairs){

    DistributedFullGrid<CombiDataType>& dfg_t = pair[0]->getDistributedFullGrid();
    DistributedFullGrid<CombiDataType>& dfg_s = pair[1]->getDistributedFullGrid();

    dfg_t.addToUniformSG( *SDCUniDSG, 1.0 );
    dfg_s.addToUniformSG( *SDCUniDSG, -1.0 );

    CombiDataType localBetaMax(0.0);

    LevelVector subMax;
    size_t jMax = 0;

    for (size_t i = 0; i < SDCUniDSG->getNumSubspaces(); ++i){
      auto subData = SDCUniDSG->getData(i);
      auto subSize = SDCUniDSG->getDataSize(i);
      for (size_t j = 0; j < subSize; ++j){
        if (std::abs(subData[j]) >= std::abs(localBetaMax)){
          localBetaMax = subData[j];
          subMax = SDCUniDSG->getLevelVector(i);
          jMax = j;
        }
      }
    }

    allBetas.push_back( localBetaMax );
    allSubs.push_back( subMax );
    allJs.push_back( jMax );

    // Reset sparse grid to zero before continuing to next pair
    for (size_t i = 0; i < SDCUniDSG->getNumSubspaces(); ++i){
      auto subData = SDCUniDSG->getData(i);
      auto subSize = SDCUniDSG->getDataSize(i);
      for (size_t j = 0; j < subSize; ++j)
        subData[j] = 0.0;
    }
  }

  // We gather all beta values across processes in a group and store only the largest per pair
  std::vector<CombiDataType> allBetasReduced;
  allBetasReduced.resize(allBetas.size());

  CombiCom::BetasReduce( allBetas, allBetasReduced, theMPISystem()->getLocalComm() );

  // Largest beta in a given group -> potential candidate for SDC
  auto globalBetaMax = std::max_element(allBetasReduced.begin(), allBetasReduced.end(),
      [](CombiDataType a, CombiDataType b){ return std::abs(a) < std::abs(b); } );

  auto b = std::find( allBetas.begin(), allBetas.end(), *globalBetaMax );

  betas_.clear();
  subspaceValues_.clear();

  // The following is only done by the process that contains the suspicious value
  if(b != allBetas.end()) {

    size_t indMax = std::distance(allBetas.begin(), b);

    LevelVector subMax = allSubs[indMax];
    std::cout<< "Subspace with potential SDC = "<<subMax<<std::endl;

    size_t jMax = allJs[indMax];

    // First strategy: compare function values directly and try to fit a constant
    if ( method == COMPARE_VALUES ) {

      for (auto t: tasks_){

        LevelVector level = t->getLevelVector();

        DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid();

        dfg.addToUniformSG( *SDCUniDSG, 1.0 );

        auto subData = SDCUniDSG->getData(subMax);
        CombiDataType localValMax = subData[jMax];

        if ( subMax <= level )
          subspaceValues_[level] = localValMax;

        // Reset sparse grid to zero
        for (size_t i = 0; i < SDCUniDSG->getNumSubspaces(); ++i){
          auto subData = SDCUniDSG->getData(i);
          auto subSize = SDCUniDSG->getDataSize(i);
          for (size_t j = 0; j < subSize; ++j)
            subData[j] = 0.0;
        }
      }
      if( subspaceValues_.size() < 5 ){
        std::cout<<"To few measurements. Changing strategy to search pairs."<<std::endl;
        method = COMPARE_PAIRS;
      } else {
        robustRegressionValues( levelsSDC );
      }
    }

    // Second strategy: compare beta values and try to fit an error expansion model
    if ( method == COMPARE_PAIRS ) {

      for (auto pair : allPairs){

        DistributedFullGrid<CombiDataType>& dfg_t = pair[0]->getDistributedFullGrid();
        DistributedFullGrid<CombiDataType>& dfg_s = pair[1]->getDistributedFullGrid();

        LevelVector t_level = pair[0]->getLevelVector();
        LevelVector s_level = pair[1]->getLevelVector();

        dfg_t.addToUniformSG( *SDCUniDSG, 1.0 );
        dfg_s.addToUniformSG( *SDCUniDSG, -1.0 );

        auto subData = SDCUniDSG->getData(subMax);
        CombiDataType localBetaMax = subData[jMax];
        betas_[std::make_pair(t_level, s_level)] = localBetaMax;

        // Reset sparse grid to zero
        for (size_t i = 0; i < SDCUniDSG->getNumSubspaces(); ++i){
          auto subData = SDCUniDSG->getData(i);
          auto subSize = SDCUniDSG->getDataSize(i);
          for (size_t j = 0; j < subSize; ++j)
            subData[j] = 0.0;
        }
      }
    robustRegressionPairs( levelsSDC );
    }
  }
  // todo: Barrier needed?
  MPI_Barrier(theMPISystem()->getLocalComm());
}

void ProcessGroupWorker::computeStandardizedResiduals( gsl_multifit_robust_workspace* regressionWsp, gsl_vector* r_stud, gsl_vector* r_lms ){

  size_t p = regressionWsp->p;
  size_t n = regressionWsp->n;
  gsl_vector *r = gsl_vector_alloc( n );
  gsl_vector *r2 = gsl_vector_alloc( n );
  gsl_vector *r2_sorted = gsl_vector_alloc( n );
  gsl_vector *r_stand = gsl_vector_alloc( n );
  gsl_vector *weights = gsl_vector_alloc( n );
  gsl_vector *ones = gsl_vector_alloc( n );

  gsl_vector_memcpy(r, r_stud);
  gsl_vector_memcpy(r_stand, r_stud);
  gsl_vector_memcpy(r_lms, r_stud);

  for (size_t i = 0; i < r->size; ++i)
    gsl_vector_set(r2, i, std::pow(gsl_vector_get(r, i),2));

  gsl_vector_memcpy(r2_sorted, r2);
  gsl_sort_vector(r2_sorted);

  // Median of residuals squared
  double median_r2 = gsl_stats_median_from_sorted_data(r2_sorted->data, r2_sorted->stride, r2_sorted->size);

  // Preliminary scale estimate
  double s0 = 1.4826*( 1 + 5.0/(n-p-1))*(std::sqrt(median_r2));

  // Standardized residuals
  gsl_vector_scale(r_stand, 1.0/s0);

  for (size_t i = 0; i < r_stand->size; ++i)
    gsl_vector_set(r_stand, i, fabs(gsl_vector_get(r_stand, i)));

  // Threshold for outliers
  double eps = 2.5;

  // Weights for each residual: here we assume that the can
  // be multiple residuals
  for(size_t i = 0; i < r_stand->size; ++i){
    if(std::abs(r_stand->data[i]) <= eps)
      gsl_vector_set(weights, i, 1);
    else
      gsl_vector_set(weights, i, 0);
  }

  // Robust scale estimate
  double prod;
  gsl_blas_ddot( weights, r2, &prod );

  double sum;
  gsl_vector_set_all( ones, 1 );
  gsl_blas_ddot( weights, ones, &sum );
  double s_star = std::sqrt(prod/(sum-p));

  gsl_vector_scale( r_lms, 1.0/s_star );
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
      ){
        allPairs.push_back({currentPair[1],currentPair[0]});
        numPairs[s->getLevelVector()]++;
        numPairs[t->getLevelVector()]++;
        k++;
      }

      if (k == numNearestNeighbors)
        break;
    }
  }
  // Check if any grid was left out with fewer neighbors than it should
  for (Task* s : tasks_ ){

    int k = numPairs[s->getLevelVector()];
    if ( k < numNearestNeighbors ){

      std::sort(levels.begin(), levels.end(), [s](LevelVector const& a, LevelVector const& b) {
        return l1(a - s->getLevelVector()) < l1(b - s->getLevelVector());
      });

      for( size_t t_i = 1; t_i < levels.size(); ++t_i ){
        std::vector<Task*> currentPair;
        std::vector<Task*> currentPairBack;

        Task* t = *std::find_if(tasks_.begin(), tasks_.end(),
            [levels,t_i](Task* const &tt) -> bool { return tt->getLevelVector() == levels[t_i]; });

        currentPair.push_back(s);
        currentPair.push_back(t);
        currentPairBack.push_back(t);
        currentPairBack.push_back(s);

        if(std::find(allPairs.begin(), allPairs.end(), currentPair) == allPairs.end() &&
            std::find(allPairs.begin(), allPairs.end(), currentPairBack) == allPairs.end()){
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
}

void ProcessGroupWorker::robustRegressionPairs( std::vector<int> &levelsSDC ){

  auto dim  = combiParameters_.getDim();
  auto lmin = combiParameters_.getLMin();
  auto lmax = combiParameters_.getLMax();

  // Determine which indices appear in the set of multi-indices
  std::set<IndexType> indexSet;
   for ( auto t : tasks_ ){
     LevelVector key = t->getLevelVector();
     for( IndexType k_i : key )
       indexSet.insert( k_i );
   }

  std::map<IndexType, IndexType> indexMap;
  IndexType row = 0;
  for ( auto ii : indexSet ){
    indexMap[ii] = row;
    row++;
  }

  // Number of measurements (beta values)
  size_t n = betas_.size();

  // Number of unknowns (values of functions D_i)
  size_t numIndices = indexSet.size();
  size_t p = dim*numIndices;

  // Exponent in error splitting
  double ex = 2;

  if ( n < p ){
    std::cout<<"Too few measurements: SDC detection skipped."<<std::endl;
    return;
  }

  gsl_multifit_robust_workspace *regressionWsp = gsl_multifit_robust_alloc(gsl_multifit_robust_cauchy, n , p );

  gsl_matrix *X = gsl_matrix_alloc( n, p );
  gsl_vector *y = gsl_vector_alloc( n );
  gsl_vector *c = gsl_vector_alloc( p );
  gsl_vector *r_stud = gsl_vector_alloc( n );
  gsl_matrix *cov = gsl_matrix_alloc( p, p );

  // Initialize matrix with zeros
  gsl_matrix_set_zero( X );

  IndexType idx_Di;
  CombiDataType val;
  row  = 0;
  for( auto const &entry : betas_ ){
    LevelVector key_t = entry.first.first;
    LevelVector key_s = entry.first.second;
    CombiDataType beta = entry.second;

    for( size_t di = 0; di < dim; ++di){
      CombiDataType ht_i = 1.0/pow(2.0,static_cast<double>(key_t[di]));
      idx_Di = di*numIndices + indexMap[key_t[di]];
      gsl_matrix_set( X, row, idx_Di, pow(ht_i,ex) );

      CombiDataType hs_i = 1.0/pow(2.0,static_cast<double>(key_s[di]));
      idx_Di = di*numIndices + indexMap[key_s[di]];
      val = gsl_matrix_get( X, row, idx_Di );
      gsl_matrix_set( X, row, idx_Di, val - pow(hs_i,ex) );
    }

    gsl_vector_set( y, row, beta );

    row++;
  }

  gsl_set_error_handler_off();
  gsl_multifit_robust( X, y, c, cov, regressionWsp );

  gsl_multifit_robust_residuals(X, y, c, r_stud, regressionWsp);

  // Before checking the residuals, check the data for extremely large values
  double eps = 1e50;
  detectOutliers( y->data, levelsSDC, eps, COMPARE_PAIRS);

  // Now we can check for large residuals
  if(levelsSDC.size() == 0){
    eps = 2.5;
    detectOutliers( r_stud->data, levelsSDC, eps, COMPARE_PAIRS );
  }

  // Print betas and resiguals
  row = 0;
  for ( auto const &entry : betas_ ){
    LevelVector key_t = entry.first.first;
    LevelVector key_s = entry.first.second;
    CombiDataType beta = entry.second;
    CombiDataType res = r_stud->data[row];
    std::cout<<key_t <<","<< key_s <<","<<beta<<","<< res <<std::endl;
    row++;
  }

  gsl_matrix_free(X);
  gsl_matrix_free(cov);
  gsl_vector_free(y);
  gsl_vector_free(c);
  gsl_multifit_robust_free(regressionWsp);
}

void ProcessGroupWorker::robustRegressionValues( std::vector<int> &levelsSDC ){

  auto lmin = combiParameters_.getLMin();
  auto lmax = combiParameters_.getLMax();

  // Number of measurements (beta values)
  size_t n = subspaceValues_.size();

  // Number of unknowns
  size_t p = 1;

  gsl_multifit_robust_workspace *regressionWsp = gsl_multifit_robust_alloc(gsl_multifit_robust_cauchy, n , p );

  gsl_matrix *X = gsl_matrix_alloc( n, p );
  gsl_vector *y = gsl_vector_alloc( n );
  gsl_vector *c = gsl_vector_alloc( p );
  gsl_vector *r_lms = gsl_vector_alloc( n );
  gsl_vector *r_stud = gsl_vector_alloc( n );
  gsl_matrix *cov = gsl_matrix_alloc( p, p );

  // Initialize matrix with ones
  gsl_matrix_set_all( X, 1.0 );

  int row  = 0;
  for( auto const &entry : subspaceValues_){
    CombiDataType maxVal = entry.second;
    gsl_vector_set( y, row, maxVal );
    row++;
  }

  // Before checking the residuals, check the data for extremely large values
  double eps = 1e50;
  detectOutliers( y->data, levelsSDC, eps, COMPARE_VALUES );

  if ( n < 5 ){
    std::cout<<"Too few measurements: Robust regression skipped."<<std::endl;
    return;
  }

  gsl_set_error_handler_off();
  gsl_multifit_robust( X, y, c, cov, regressionWsp );

  gsl_multifit_robust_residuals(X, y, c, r_stud, regressionWsp);

  computeStandardizedResiduals( regressionWsp, regressionWsp->r, r_lms );

  // Now we can check for large residuals
  if(levelsSDC.size() == 0){
    eps = 2.5;
    detectOutliers( r_lms->data, levelsSDC, eps, COMPARE_VALUES, y->data[0] );
  }

  // Print data and residuals
  row = 0;
  for ( auto const &entry : subspaceValues_){
    LevelVector key = entry.first;
    CombiDataType maxVal = entry.second;
    CombiDataType res = r_lms->data[row];
    std::cout<<key <<","<<maxVal<<","<< res <<std::endl;
    row++;
  }

  gsl_matrix_free(X);
  gsl_matrix_free(cov);
  gsl_vector_free(y);
  gsl_vector_free(c);
  gsl_multifit_robust_free(regressionWsp);
}

void ProcessGroupWorker::detectOutliers( double* data, std::vector<int> &levelsSDC, double eps, SDCMethodType method, double y ){

  std::map<LevelVector,int> mapSDC;
  for ( auto t : tasks_ ){
    LevelVector key = t->getLevelVector();
    mapSDC[key] = 0;
  }

  if ( method == COMPARE_PAIRS ) {
    size_t numSDCPairs = 0;
    size_t row = 0;
    for( auto const &entry : betas_ ){
      LevelVector key_t = entry.first.first;
      LevelVector key_s = entry.first.second;
      CombiDataType beta = entry.second;
      if ( std::abs(data[row]) > eps && beta != 0 ){
        mapSDC[key_t]++;
        mapSDC[key_s]++;
        numSDCPairs++;
      }
      row++;
    }
    std::cout<< "SDC grid: " << std::endl;
    for (auto s : mapSDC){
      if (s.second >= 2 || (s.second == 1 && numSDCPairs == 1)){
        std::cout<<s.first<<std::endl;
        int id = combiParameters_.getID(s.first);
        levelsSDC.push_back(id);
      }
    }
  }
  if ( method == COMPARE_VALUES ) {
    size_t row = 0;
    for( auto const &entry : subspaceValues_ ){
      LevelVector key = entry.first;
      if ( std::abs(data[row]) > eps ){
        mapSDC[key]++;
      }
      row++;
    }
    std::cout<< "SDC grid: " << std::endl;
    for (auto s : mapSDC){
      if (s.second == 1 ){
        std::cout<<s.first<<std::endl;
        int id = combiParameters_.getID(s.first);
        levelsSDC.push_back(id);
      }
    }
    // Check for false positives and remove them if it's the case
    removeFalsePositives( levelsSDC, y );
  }
}

void ProcessGroupWorker::removeFalsePositives( std::vector<int>& faultsID, double u_robust ) {

  std::map<int, LevelVector> IDsToLevels = combiParameters_.getLevelsDict();
  int dim = static_cast<int>(combiParameters_.getDim());

  const std::string prob_name = "interpolation based optimization";

  LevelVectorList lvlminmax;
  lvlminmax.push_back(combiParameters_.getLMin());
  lvlminmax.push_back(combiParameters_.getLMax());

  CombigridDict given_dict = combiParameters_.getCombiDict();

  std::vector<int> faultsIDCopy = faultsID;

  // Iterate through tasks currently marked as outliers
  for (auto id : faultsIDCopy){

    LevelVectorList faultLevelVectors;
    faultLevelVectors.push_back(IDsToLevels[id]);

    // Find alternative combination coefficients to exclude suspect value
    LP_OPT_INTERP opt_interp(lvlminmax,dim,GLP_MAX,given_dict,faultLevelVectors);
    opt_interp.init_opti_prob(prob_name);
    opt_interp.set_constr_matrix();
    opt_interp.solve_opti_problem();

    LevelVectorList recomputeLevelVectors;
    CombigridDict new_dict = opt_interp.get_results(recomputeLevelVectors);

    CombiDataType u_c = 0.0;
    CombiDataType ut_c = 0.0;
    CombiDataType rel_err_max = 0.10;

    // Combine values with and without the suspect
    for( auto const &entry : given_dict){
      LevelVector lvl = entry.first;
      CombiDataType u_i;
      if ( subspaceValues_.find(lvl) == subspaceValues_.end() )
        u_i = u_robust;
      else
        u_i = subspaceValues_[lvl];
      u_c += given_dict[lvl]*u_i;
      ut_c += new_dict[lvl]*u_i;
    }

    // Calculate relative error
    CombiDataType rel_err = std::abs( u_c - ut_c )/u_robust;

    // If relative error is small, value was a false positive (not an outlier)
    // -> Remove ID from list of SDC IDs
    if ( rel_err < rel_err_max ){
      faultsID.erase(std::remove(faultsID.begin(), faultsID.end(), id), faultsID.end());
      std::cout<<"Removing id = "<<id<<std::endl;
    }
  }
}
} /* namespace combigrid */
