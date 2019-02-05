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

namespace combigrid {

ProcessGroupWorker::ProcessGroupWorker() :
    currentTask_( NULL),
    status_(PROCESS_GROUP_WAIT),
    combinedFG_( NULL),
    combinedUniDSG_(NULL),
    thirdLevel_(NULL),
    combinedFGexists_(false),
    combiParameters_(),
    combiParametersSet_(false)
{
}

ProcessGroupWorker::~ProcessGroupWorker() {
  delete combinedFG_;
  delete thirdLevel_;
}

SignalType ProcessGroupWorker::wait() {
  if (status_ != PROCESS_GROUP_WAIT)
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
    Stats::startEvent("worker init");
    currentTask_->init(theMPISystem()->getLocalComm());
    Stats::stopEvent("worker init");

    // execute task
    Stats::startEvent("worker run first");
    currentTask_->run(theMPISystem()->getLocalComm());
    Stats::stopEvent("worker run first");

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
    Stats::startEvent("worker run");
    currentTask_->run(theMPISystem()->getLocalComm());
    Stats::stopEvent("worker run");

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

    Stats::startEvent("combine");
    combineUniform();
    Stats::stopEvent("combine");

//    if (theMPISystem()->getWorldRank() == 0) {
//      std::cout << "DSG on rank " << theMPISystem()->getWorldRank() << std::endl
//                << "-------------" << std::endl;
//      for (size_t i = 0; i < combinedUniDSG_->getNumSubspaces(); i++) {
//        std::cout << "Subspace with level " << combinedUniDSG_->getLevelVector(i)
//                  << "has " << combinedUniDSG_->getDataVector(i).size()
//                  << "elements" << std::endl;
//      }
//    }

  } else if (signal == COMBINE_UNIFORM_THIRD_LEVEL) {

    Stats::startEvent("combineUniformThirdLevelRecvFirst");
    combineUniformThirdLevelRecvFirst<CombiDataType>();
    Stats::stopEvent("combineUniformThirdLevelRecvFirst");

  } else if (signal == EXCHANGE_COMMON_SS) {

    Stats::startEvent("combineUniformThirdLevelSendFirst");
    combineUniformThirdLevelSendFirst<CombiDataType>();
    Stats::stopEvent("combineUniformThirdLevelSendFirst");

  } else if (signal == COMBINE_LOCAL_AND_GLOBAL) {

    Stats::startEvent("combineLocalAndGlobal");
    combineLocalAndGlobal();
    Stats::stopEvent("combineLocalAndGlobal");

  } else if (signal == INTEGRATE_COMMON_SS) {

    Stats::startEvent("integrateCommonSS");
    integrateCommonSS<CombiDataType>();
    Stats::stopEvent("integrateCommonSS");

  } else if (signal == GRID_EVAL) {

    Stats::startEvent("eval");
    gridEval();
    Stats::stopEvent("eval");

    return signal;

  } else if (signal == COMBINE_FG) {

    combineFG();

  } else if (signal == UPDATE_COMBI_PARAMETERS) {

    updateCombiParameters();

  } else if( signal == PARALLEL_EVAL ){

    Stats::startEvent("parallel eval");
    parallelEval();
    Stats::stopEvent("parallel eval");

  }

  // special solution for GENE
  // todo: find better solution and remove this
  if( ( signal == RUN_FIRST || signal == RUN_NEXT ) && omitReadySignal )
    return signal;

  ready();

  return signal;
}

void ProcessGroupWorker::ready() {
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
      }
    }

    // all tasks finished -> group waiting
    status_ = PROCESS_GROUP_WAIT;
  }

  // send ready status to manager
  MASTER_EXCLUSIVE_SECTION{
    StatusType status = status_;

    MPI_Send(&status, 1, MPI_INT, theMPISystem()->getManagerRank(), statusTag, theMPISystem()->getGlobalComm());
  }

  // reset current task
  currentTask_ = NULL;
}

void ProcessGroupWorker::combine() {
  assert( false && "not properly implemented" );

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
  // each pgrouproot must call reduce function
  assert(tasks_.size() > 0);

  assert( combiParametersSet_ );
  DimType dim = combiParameters_.getDim();
  const LevelVector& lmin = combiParameters_.getLMin();
  LevelVector lmax = combiParameters_.getLMax();
  const std::vector<bool>& boundary = combiParameters_.getBoundary();

  // the dsg can be smaller than lmax because the highest subspaces do not have
  // to be exchanged
  // todo: use a flag to switch on/off optimized combination
  /*
  for (size_t i = 0; i < lmax.size(); ++i)
    if (lmax[i] > lmin[i])
      lmax[i] -= 1;
      */

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
    DistributedHierarchization::hierarchize<CombiDataType>(
        dfg, combiParameters_.getHierarchizationDims() );

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
    DistributedHierarchization::dehierarchize<CombiDataType>(
        dfg, combiParameters_.getHierarchizationDims() );
  }

}

void ProcessGroupWorker::combineLocalAndGlobal() {
  // each pgrouproot must call reduce function
  assert(tasks_.size() > 0);

  assert( combiParametersSet_ );
  DimType dim = combiParameters_.getDim();
  const LevelVector& lmin = combiParameters_.getLMin();
  LevelVector lmax = combiParameters_.getLMax();
  const std::vector<bool>& boundary = combiParameters_.getBoundary();

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
    DistributedHierarchization::hierarchize<CombiDataType>(
        dfg, combiParameters_.getHierarchizationDims() );

    // lokales reduce auf sg
    dfg.addToUniformSG( *combinedUniDSG_, combiParameters_.getCoeff( t->getID() ) );
  }

  // global reduce
  CombiCom::distributedGlobalReduce( *combinedUniDSG_ );
}

template <typename FG_ELEMENT>
void ProcessGroupWorker::integrateCommonSS() {
  const std::vector<LevelVector> commonSS = combiParameters_.getCommonSubspaces();
  const size_t numCommonSS = commonSS.size();

  // broadcast common subspaces to global comm
  MPI_Datatype dtype = abstraction::getMPIDatatype(
                         abstraction::getabstractionDataType<FG_ELEMENT>());
  for (size_t ss = 0; ss < numCommonSS; ss++) {
    size_t ssSize = combinedUniDSG_->getDataSize(commonSS[ss]);
    FG_ELEMENT* ssData = combinedUniDSG_->getData(commonSS[ss]);
    MPI_Bcast(ssData, ssSize, dtype, 0, theMPISystem()->getGlobalReduceComm());
  }

  // extract dfg and dehierarchize
  for (Task* t : tasks_) {
    // get handle to dfg
    DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid();

    // extract dfg vom dsg
    dfg.extractFromUniformSG( *combinedUniDSG_ );

    // dehierarchize dfg
    DistributedHierarchization::dehierarchize<CombiDataType>(
        dfg, combiParameters_.getHierarchizationDims() );
  }
}

template <typename FG_ELEMENT>
void ProcessGroupWorker::combineUniformThirdLevelSendFirst() {
  const std::vector<LevelVector> commonSS = combiParameters_.getCommonSubspaces();
  const size_t numCommonSS = commonSS.size();
  assert(theMPISystem()->getThirdLevelComms().size() == 1 && "init thirdLevel communicator failed");
  const CommunicatorType& comm = theMPISystem()->getThirdLevelComms()[0];
  const RankType& thirdLevelManager = theMPISystem()->getThirdLevelManagerRank();

  // send sizes of common subspace parts
  std::vector<int> commonSSPartSizes(0);
  for (size_t ss = 0; ss < numCommonSS; ss++)
    commonSSPartSizes[ss] = static_cast<int>(combinedUniDSG_->getDataSize(commonSS[ss]));
  MPI_Send(commonSSPartSizes.data(), static_cast<int>(numCommonSS), MPI_INT,
      thirdLevelManager, 0, comm);

  // send common subspaces parts sequentially
  MPI_Datatype dtype = abstraction::getMPIDatatype(
                         abstraction::getabstractionDataType<FG_ELEMENT>());
  for (size_t ss = 0; ss < numCommonSS; ss++)
    MPI_Send(combinedUniDSG_->getData(commonSS[ss]), commonSSPartSizes[ss],
        dtype, thirdLevelManager, 0, comm);

  // receive and integrate combined common subspace parts from remote
  for (size_t ss = 0; ss < numCommonSS; ss++)
    MPI_Recv(combinedUniDSG_->getData(commonSS[ss]), commonSSPartSizes[ss],
        dtype, thirdLevelManager, 0, comm, MPI_STATUS_IGNORE);
}

template <typename FG_ELEMENT>
void ProcessGroupWorker::combineUniformThirdLevelRecvFirst() {
  const std::vector<LevelVector> commonSS = combiParameters_.getCommonSubspaces();
  const size_t numCommonSS = commonSS.size();
  assert(theMPISystem()->getThirdLevelComms().size() == 1 && "init thirdLevel communicator failed");
  const CommunicatorType& comm = theMPISystem()->getThirdLevelComms()[0];
  const RankType& thirdLevelManager = theMPISystem()->getThirdLevelManagerRank();

  // send sizes of common subspace parts
  std::vector<int> commonSSPartSizes(0);
  for (size_t ss = 0; ss < numCommonSS; ss++)
    commonSSPartSizes[ss] = static_cast<int>(combinedUniDSG_->getDataSize(commonSS[ss]));
  MPI_Send(commonSSPartSizes.data(), static_cast<int>(numCommonSS), MPI_INT,
      thirdLevelManager, 0, comm);

  // allreduce common subspace parts from remote with local
  MPI_Datatype dtype = abstraction::getMPIDatatype(
                         abstraction::getabstractionDataType<FG_ELEMENT>());
  for (size_t ss = 0; ss < numCommonSS; ss++)
    MPI_Allreduce(MPI_IN_PLACE, combinedUniDSG_->getData(commonSS[ss]), commonSSPartSizes[ss], dtype, MPI_SUM, comm);
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
  const int dim = static_cast<int>( combiParameters_.getDim() );

  // combine must have been called before this function
  assert( combinedUniDSG_ != NULL && "you must combine before you can eval" );

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


  // create dfg
  DistributedFullGrid<CombiDataType> dfg( dim, leval,
                                          combiParameters_.getApplicationComm(),
                                          combiParameters_.getBoundary(),
                                          combiParameters_.getParallelization(),
                                          false
                                          );

  // register dsg
  dfg.registerUniformSG(*combinedUniDSG_);

  // fill dfg with hierarchical coefficients from distributed sparse grid
  dfg.extractFromUniformSG( *combinedUniDSG_ );

  // dehierarchize dfg
  DistributedHierarchization::dehierarchize<CombiDataType>(
      dfg, combiParameters_.getHierarchizationDims() );

  // save dfg to file with MPI-IO
  dfg.writePlotFile( filename.c_str() );
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

} /* namespace combigrid */
