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

#include "sgpp/distributedcombigrid/mpi_fault_simulator/MPI-FT.h"

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
  MASTER_EXCLUSIVE_SECTION {
    std::string fname ="out/all-betas-"+std::to_string(theMPISystem()->getGlobalRank())+".txt";
    //betasFile_ = std::ofstream( fname, std::ofstream::out );
  }
}

ProcessGroupWorker::~ProcessGroupWorker() {
  delete combinedFG_;
}

SignalType ProcessGroupWorker::wait() {
  if(status_ == PROCESS_GROUP_FAIL){ //exit all procs in failed group
    return EXIT;
  }
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
  //std::cout << theMPISystem()->getWorldRank() << " waits for signal " << signal << " \n";
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
    t->init( theMPISystem()->getLocalComm());

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
  } else if ( signal ==  RECOVER_COMM ){
    theMPISystem()->recoverCommunicators( true );
    return signal;
  } else if( signal == PARALLEL_EVAL ){

    parallelEval();
  }

  // special solution for GENE
  // todo: find better solution and remove this
  if( ( signal == RUN_FIRST || signal == RUN_NEXT ) && omitReadySignal )
    return signal;

  // in the general case: send ready signal.
  //if(!omitReadySignal)
  ready();

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
        currentTask_->run(theMPISystem()->getLocalComm());

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
        return;
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
    if( status_ == PROCESS_GROUP_FAIL )
      theMPISystem()->recoverCommunicators( false );
  }
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
  // early exit if no tasks available
  // todo: doesnt work, each pgrouproot must call reduce function
  assert(tasks_.size() > 0);

  assert( combiParametersSet_ );
  DimType dim = combiParameters_.getDim();
  const LevelVector& lmin = combiParameters_.getLMin();
  LevelVector lmax = combiParameters_.getLMax();
  const std::vector<bool>& boundary = combiParameters_.getBoundary();

  // the dsg can be smaller than lmax because the highest subspaces do not have
  // to be exchanged
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

  real localMax(0.0);

  for (Task* t : tasks_) {

    DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid();

    // compute max norm
    /*
    real max = dfg.getLpNorm(0);
    if( max > localMax )
      localMax = max;
      */

    // hierarchize dfg
    DistributedHierarchization::hierarchize<CombiDataType>(
        dfg, combiParameters_.getHierarchizationDims() );

    // lokales reduce auf sg ->
    dfg.addToUniformSG( *combinedUniDSG_, combiParameters_.getCoeff( t->getID() ) );
  }

  // compute global max norm
  /*
  real globalMax_tmp;
  MPI_Allreduce(  &localMax, &globalMax_tmp, 1, MPI_DOUBLE,
                  MPI_MAX, theMPISystem()->getGlobalReduceComm() );

  real globalMax;
  MPI_Allreduce(  &globalMax_tmp, &globalMax, 1, MPI_DOUBLE,
                    MPI_MAX, theMPISystem()->getLocalComm() );
                    */

  CombiCom::distributedGlobalReduce( *combinedUniDSG_ );

  for (Task* t : tasks_) {
    // get handle to dfg
    DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid();

    // extract dfg vom dsg
    dfg.extractFromUniformSG( *combinedUniDSG_ );

    // dehierarchize dfg
    DistributedHierarchization::dehierarchize<CombiDataType>(
        dfg, combiParameters_.getHierarchizationDims() );

    // if exceeds normalization limit, normalize dfg with global max norm
    /*
    if( globalMax > 1000 ){
      dfg.mul( 1.0 / globalMax );
      std::cout << "normalized dfg with " << globalMax << std::endl;
    }
    */

  }

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
      // todo: remove
      // save fg file for debugging
      //std::string tmp( "fg.dat" );
      //fg.save( tmp );

      // todo: remove. works only for small grids. also hard coded path
      /*
      int id = t->getID();
      std::string path = std::string("/home/heenemo/workspace/combi-gene/distributedcombigrid/examples/gene_distributed/")
                         + "ginstance" + boost::lexical_cast<std::string>( id )
                         + "/fg.dat";
      FullGrid<CombiDataType> fg_plot(dim, leval, boundary);
      fg_plot.add( fg, 1.0 );
      fg_plot.writePlotFile( path.c_str() );
      */

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
