/*
 * GeneTask.cpp
 *
 *  Created on: Jul 10, 2014
 *      Author: heenemo
 */

#include "GeneTask.hpp"
#include <sstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <unistd.h>
#include <fstream>
#include "CombiGeneConverter.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupSignals.hpp"
#include <math.h>

//#include "sgpp/distributedcombigrid/utils/StatsContainer.hpp"

namespace combigrid
{

GeneTask::GeneTask( DimType dim, LevelVector& l,
                    std::vector<bool>& boundary, real coeff, LoadModel* loadModel,
                    std::string& path, real dt, size_t nsteps,
                    real shat, real kymin, real lx, int ky0_ind,
                    IndexVector p , FaultCriterion *faultCrit,
                    IndexType numSpecies, bool GENE_Global, bool GENE_Linear)
    : Task( dim, l, boundary, coeff, loadModel,faultCrit),
      path_( path ),
      dt_( dt ),
      nsteps_( nsteps ),
      stepsTotal_(0),
      combiStep_(0),
      shat_( shat ),
      kymin_( kymin ),
      lx_( lx ),
      ky0_ind_( ky0_ind ),
      x0_(lx/2.0),
      p_(p),
      checkpoint_(), initialized_(false),
      checkpointInitialized_(false),
      nspecies_(numSpecies),
      _GENE_Global(GENE_Global),
      _GENE_Linear(GENE_Linear)
{

// theres only one boundary configuration allowed at the moment
assert( boundary[0] == true );//x
assert( boundary[1] == false );//y
assert( boundary[2] == true );//z
assert( boundary[3] == true );//v
assert( boundary[4] == true );//w
assert( boundary[5] == false );//nspec
}

GeneTask::GeneTask() :
    checkpoint_(),
    dfgVector_(0),
    nrg_(0.0),
    initialized_(false),
    checkpointInitialized_(false),
    _GENE_Global(false),
    _GENE_Linear(true)
{
  ;
}

GeneTask::~GeneTask()
{
  for(unsigned int g = 0; g < dfgVector_.size(); g++){
    delete dfgVector_[g];
  }
  dfgVector_.clear();
}

/**
 * This routine is used to initialize everything needed for running GENE.
 * This routine does NOT actually run GENE. It only changes directories
 */
void
GeneTask::run( CommunicatorType lcomm )
{
  using namespace std::chrono;

  // change dir to wdir
  if( chdir( path_.c_str() ) ){


    printf( "could not change to directory %s \n", path_.c_str() );
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }

//  char cwd[1024];
//  getcwd(cwd, sizeof(cwd));

  // it is more save to wait here until all procs are in the right directory
  MPI_Barrier( lcomm );

  MASTER_EXCLUSIVE_SECTION{
    std::cout << "run task " << this->getID() << std::endl;
  }
//  int globalRank;
//  // MPI_Comm_rank(lcomm, &lrank);
//  MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);
//  if(combiStep_ != 0){
//    //theStatsContainer()->setTimerStop("computeIterationRank" + std::to_string(globalRank));
//    //theStatsContainer()->setValue("computeIterationRank" + std::to_string(globalRank),0.0);
//  }
//  //startTimeIteration_ = high_resolution_clock::now();
//
//  //theStatsContainer()->setTimerStart("computeIterationRank" + std::to_string(globalRank));
//
//  //printf("running gene!!! \n");

}
/**
 * This routine is used to change the directory to the directory of the current task
 */
void GeneTask::changeDir(CommunicatorType lcomm){
  // change dir to wdir
  if( chdir( path_.c_str() ) ){
    printf( "could not change to directory %s \n", path_.c_str() );
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }

//  char cwd[1024];
//  getcwd(cwd, sizeof(cwd));
  // it is more save to wait here until all procs are in the right directory
  MPI_Barrier( lcomm );

  MASTER_EXCLUSIVE_SECTION{
    std::cout << "changed to task " << this->getID() << std::endl;
  }
}
/**
 * This method is used to kill a process according to the faultCriterion.
 * Failing processes call the Sim_FT_kill_me() function
 */
void GeneTask::decideToKill(){ //toDo check if combiStep should be included in task and sent to process groups in case of reassignment
  using namespace std::chrono;

  int globalRank;
  // MPI_Comm_rank(lcomm, &lrank);
  MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);
  //theStatsContainer()->setTimerStop("computeIterationRank" + std::to_string(globalRank));
  //duration<real> dur = high_resolution_clock::now() - startTimeIteration_;
  //real t_iter = dur.count();
  //std::cout << "Current iteration took " << t_iter << "\n";

  //theStatsContainer()->setTimerStart("computeIterationRank" + std::to_string(globalRank));


  //check if killing necessary
  //std::cout << "failNow result " << failNow(globalRank) << " at rank: " << globalRank <<" at step " << combiStep_ << "\n" ;
  //real t = dt_ * nsteps_ * combiStep_;
  if (combiStep_ != 0 && faultCriterion_->failNow(combiStep_, -1.0, globalRank)){
        std::cout<<"Rank "<< globalRank <<" failed at iteration "<<combiStep_<<std::endl;
        StatusType status=PROCESS_GROUP_FAIL;
        MASTER_EXCLUSIVE_SECTION{
          simft::Sim_FT_MPI_Send( &status, 1, MPI_INT,  theMPISystem()->getManagerRank(), statusTag,
                            theMPISystem()->getGlobalCommFT() );
        }
        theMPISystem()->sendFailedSignal();
        simft::Sim_FT_kill_me();
  }
  combiStep_++;
}
/**
 * This routine initializes GeneTask; currently it only sets a bool value.
 */
void GeneTask::init(CommunicatorType lcomm, std::vector<IndexVector> decomposition){
//  if( dfg_ == NULL ){
//      dfg_ = new DistributedFullGrid<CombiDataType>( dim_, l_, lcomm,
//          this->getBoundary(), p_, false);
//  }
  initialized_ = true;
}

/**
 * This routine writes the grid values from gene (data) to a checkpoint file
 * @param data GENE grid
 * @param size size of the data grid (total)
 * @param sizes size of the data grid in each dimension
 * @param bounds bounds of the grid in distributed implementation
 */
void
GeneTask::writeLocalCheckpoint( GeneComplex* data, size_t size,
                                std::vector<size_t>& sizes,
                                std::vector<size_t>& bounds )
{
  // todo: doing it like this will require two times copying
  for(unsigned int i= 0; i < sizes.size(); i++){
//    std::cout << i << " size[i]: " << sizes[i] << "\n";
    int index_l = sizes.size()- 1 - i ; // sizes is reversed order of l; i.e. l is x y z v w spec and sizes spec, w, v, z, y, x
//    std::cout << index_l << " l[i]: " << pow(2,l_[index_l]) << "\n";
    if(i==0){
      assert(sizes[0] == nspecies_);
    }else{
      if(i == 4 && _GENE_Linear){//we have only 1 coordinate in y direction in linear scenarios
        assert(sizes[i] == 1);
      }
      else{
        assert(sizes[i] == pow(2,l_[index_l]));
      }
    }
  }
  checkpoint_.writeCheckpoint( data, size, sizes, bounds );
  checkpointInitialized_= true;

}
/**
 * Initializes local checkpoint.
 * This method is needed in case tasks are redistributed or recomputed
 * and the checkpoint is not initialized during read memory
 * @param size size of data (total)
 * @param sizes size of data in each dimension
 * @param bounds bounds of data in distributed array
 */
void GeneTask::InitLocalCheckpoint(size_t size,
    std::vector<size_t>& sizes,
    std::vector<size_t>& bounds ){
  for(unsigned int i= 0; i < sizes.size(); i++){
//    std::cout << i << " size[i]: " << sizes[i] << "\n";
    int index_l = sizes.size()- 1 - i ; // sizes is reversed order of l; i.e. l is x y z v w spec and sizes spec, w, v, z, y, x
//    std::cout << index_l << " l[i]: " << pow(2,l_[index_l]) << "\n";
    if(i==0){
      assert(sizes[0] == nspecies_); // we have nspecies elements in this dimension
    }
    else{
      if(i == 4){ //we have only 1 coordinate in y direction
        assert(sizes[i] == 1);
      }
      else{
        assert(sizes[i] == pow(2,l_[index_l])); //check if parameters match
      }
    }
  }
  checkpoint_.initCheckpoint(size,sizes,bounds);
  checkpointInitialized_= true;
}

/**
 * This routine returns the complete fullgrid on rank lroot
 * Therefore the fullgrid is collected from the other ranks of the local communicator
 */
void GeneTask::getFullGrid( FullGrid<CombiDataType>& fg, RankType lroot,
                            CommunicatorType lcomm, int species)
{
  dfgVector_[species]->gatherFullGrid(fg, lroot);
}

/**
 * This routine returns the local part of the fullgrid
 */
DistributedFullGrid<complex>& GeneTask::getDistributedFullGrid(int specie){
  return *dfgVector_[specie];
}


void
GeneTask::setLocalCheckpoint( FullGrid<complex>& fg, CommunicatorType lcomm,
                              RankType localRootID )
{
// todo: check whether local cp has right size

int lsize;
RankType localID;
MPI_Comm_size( lcomm, &lsize );
MPI_Comm_rank( lcomm, &localID );

// collect checkpoint on local root
GeneLocalCheckpoint& lcp = this->getLocalCheckpoint();

// gather bounds on local root
std::vector<int> sendbounds;
std::vector<int> recvbounds;
sendbounds.assign( lcp.getBounds().begin(), lcp.getBounds().end() );

if( localID == localRootID )
  recvbounds.resize( lsize * sendbounds.size() );

MPI_Gather( &sendbounds[0], static_cast<int>( sendbounds.size() ), MPI_INT,
            &recvbounds[0], static_cast<int>( sendbounds.size() ), MPI_INT,
            localRootID, lcomm );

// create buffer for data
std::vector<GeneComplex> recvdata;
if( localID == localRootID )
  recvdata.resize( static_cast<size_t>( lsize ) * lcp.getSize() );

double tstart_scatterCP = -1;

if( localID == localRootID ){
  std::vector<int> levels( this->getLevelVector().begin(),
                           this->getLevelVector().end() );

  // start time convert fg to global cp
  double tstart_convertFGtoCP = MPI_Wtime();

  // create temporary global GeneGrid
  const std::vector<size_t>& globalSize = this->getLocalCheckpoint().getSizes();
  GeneGrid global_tmp(
      boost::extents[globalSize[0]][globalSize[1]][globalSize[2]][globalSize[3]][globalSize[4]][globalSize[5]] );

  CombiGeneConverter::FullGridToGeneGrid( fg, global_tmp );

  std::cout << "task " << this->getID() << " conversion finished" << std::endl;

  // start time scatter global cp to local cp
  tstart_scatterCP = MPI_Wtime();

  // copy back global grid data with checkpoint data
  for( size_t r = 0; r < static_cast<size_t>( lsize ); ++r ){
    // bounds of subdomain
    int* bounds = &recvbounds[r * 12];

    // data of subdomain
    GeneComplex* data = &recvdata[r * lcp.getSize()];

    // create multiarray ref to subdomain
    int ln_spec = bounds[1] - bounds[0];
    int lnw0 = bounds[3] - bounds[2];
    int lnv0 = bounds[5] - bounds[4];
    int lnz0 = bounds[7] - bounds[6];
    int lnj0 = bounds[9] - bounds[8];
    int lni0 = bounds[11] - bounds[10];

    GeneGridRef local( data,
                       boost::extents[ln_spec][lnw0][lnv0][lnz0][lnj0][lni0] );

    // set view to global array
    typedef GeneGrid::index_range range;
    GeneGrid::array_view<6>::type lview = global_tmp[boost::indices[range(
        bounds[0], bounds[1] )][range( bounds[2], bounds[3] )][range(
        bounds[4], bounds[5] )][range( bounds[6], bounds[7] )][range(
        bounds[8], bounds[9] )][range( bounds[10], bounds[11] )]];

    // copy local data to global array
    typedef GeneGrid::index index;
    for( index n = 0; n < ln_spec; ++n )
      for( index m = 0; m < lnw0; ++m )
        for( index l = 0; l < lnv0; ++l )
          for( index k = 0; k < lnz0; ++k )
            for( index j = 0; j < lnj0; ++j )
              for( index i = 0; i < lni0; ++i ){
                local[n][m][l][k][j][i].r = lview[n][m][l][k][j][i].r;
                local[n][m][l][k][j][i].i = lview[n][m][l][k][j][i].i;
              }
  }
} // local root section

// distribute global grid to processes
// global grid has been copied to recvdata on root
std::vector<GeneComplex> lcp_tmp_data( lcp.getSize() );
MPI_Scatter( &recvdata[0], static_cast<int>( lcp_tmp_data.size() ),
MPI_DOUBLE_COMPLEX,
             &lcp_tmp_data[0], static_cast<int>( lcp_tmp_data.size() ),
             MPI_DOUBLE_COMPLEX,
             localRootID, lcomm );

// create multiarray refs
const std::vector<int> bounds( lcp.getBounds().begin(), lcp.getBounds().end() );
int ln_spec = bounds[1] - bounds[0];
int lnw0 = bounds[3] - bounds[2];
int lnv0 = bounds[5] - bounds[4];
int lnz0 = bounds[7] - bounds[6];
int lnj0 = bounds[9] - bounds[8];
int lni0 = bounds[11] - bounds[10];

GeneComplex* ldata = lcp.getData();
GeneComplex* ldata_tmp = &lcp_tmp_data[0];
GeneGridRef local( ldata,
                   boost::extents[ln_spec][lnw0][lnv0][lnz0][lnj0][lni0] );
GeneGridRef local_tmp( ldata_tmp,
                       boost::extents[ln_spec][lnw0][lnv0][lnz0][lnj0][lni0] );

// copy data (combination solution) to local cp
typedef GeneGrid::index index;
if( bounds[4] == 0 ){
  for( index n = 0; n < ln_spec; ++n )
    for( index m = 0; m < lnw0; ++m )
      for( index l = 0; l < lnv0; ++l )
        for( index k = 0; k < lnz0; ++k )
          for( index j = 0; j < lnj0; ++j )
            for( index i = 0; i < lni0; ++i ){
              local[n][m][l][k][j][i].r = local_tmp[n][m][l][k][j][i].r;
              local[n][m][l][k][j][i].i = local_tmp[n][m][l][k][j][i].i;
            }
}

}

void
GeneTask::saveCheckpoint( FullGrid<CombiDataType>& fg, const char* filename )
{
// first make sure we have a valid fullgrid
assert( fg.isGridCreated() );
assert( fg.getDimension() == 6 );
const std::vector<bool>& boundary = fg.returnBoundaryFlags();
assert( boundary[0] == true ); //x
assert( boundary[1] == false ); //y
assert( boundary[2] == true ); // z
assert( boundary[3] == true ); // v
assert( boundary[4] == true ); // w
assert( boundary[5] == false ); // nspec
assert( fg.getLevels()[1] == 1 ); // y
assert( fg.getLevels()[5] == 1 ); // nspec

size_t nspec = 1;
size_t nw = static_cast<size_t>( std::pow( 2.0, fg.getLevels()[4] ) );
size_t nv = static_cast<size_t>( std::pow( 2.0, fg.getLevels()[3] ) );
size_t nz = static_cast<size_t>( std::pow( 2.0, fg.getLevels()[2] ) );
size_t ny = 1;
size_t nx = static_cast<size_t>( std::pow( 2.0, fg.getLevels()[0] ) );

// create MultiArray to store Gene grid temporarily
// note reverse notation of shape vector
IndexVector geneShape = { nspec, nw, nv, nz, ny, nx };
MultiArray<GeneComplex,6> geneGrid = createMultiArray<GeneComplex,6>( geneShape );

// create MultiArrayRef to fg
MultiArrayRef6 fgData = createMultiArrayRef<CombiDataType,6>( fg );

// copy data to gene grid
// note that the last grid points in x,z,v,w dimension are ignored
for( size_t n=0; n < geneShape[0]; ++n ){ //n_spec
  for( size_t m=0; m < geneShape[1]; ++m ){ //w
    for( size_t l=0; l < geneShape[2]; ++l ){ //v
      for( size_t k=0; k < geneShape[3]; ++k ){ //z
        for( size_t j=0; j < geneShape[4]; ++j ){ //y
          for( size_t i=0; i < geneShape[5]; ++i ){ //x
            geneGrid[n][m][l][k][j][i].r = fgData[n][m][l][k][j][i].real();
            geneGrid[n][m][l][k][j][i].i = fgData[n][m][l][k][j][i].imag();
          }
        }
      }
    }
  }
}

// save to file
std::ofstream cpFile( filename, std::ofstream::binary );

std::cout << "saving full grid as GENE checkpoint" << std::endl;

// write prec flag
const char* prec = "DOUBLE";
cpFile.write( prec, 6 );
std::cout << "prec:(6)" << prec << std::endl;

// write time and dt
// todo: use real values here
double mytime( 0 ), dt( 0 );
cpFile.write( (char*) &mytime, sizeof(double) );
cpFile.write( (char*) &dt, sizeof(double) );
std::cout << "mytime(" << sizeof(double) << "): " << mytime << std::endl;
std::cout << "dt(" << sizeof(double) << "): " << dt << std::endl;

// write resolution
int res[6];
res[0] = static_cast<int>(nx);
res[1] = static_cast<int>(ny);
res[2] = static_cast<int>(nz);
res[3] = static_cast<int>(nv);
res[4] = static_cast<int>(nw);
res[5] = static_cast<int>(nspec);
cpFile.write( (char*) res, sizeof(int) * 6 );
std::cout << "resolution:" << "\n" << "\t nx(" << sizeof(int) << ") " << nx
    << "\n" << "\t ny(" << sizeof(int) << ") " << ny << "\n" << "\t nz("
    << sizeof(int) << ") " << nz << "\n" << "\t nv(" << sizeof(int) << ") "
    << nv << "\n" << "\t nw(" << sizeof(int) << ") " << nw << "\n"
    << "\t n_spec(" << sizeof(int) << ") " << nspec << std::endl;

// write data
size_t dsize( geneGrid.num_elements() );
cpFile.write( (char*) geneGrid.origin(), sizeof(GeneComplex) * dsize );

cpFile.close();
}

/**
 * This routine sets the dfg to zero
 */
void GeneTask::setZero(){
  if(dfgVector_.size() != 0){
    for(int i=0; i< dfgVector_.size(); i++){
      std::vector<CombiDataType>& data = dfgVector_[i]->getElementVector();

      for( size_t i=0; i<data.size(); ++i ){
        data[i] = complex(0.0);
      }
    }
  }
}

/**
 * This routine initializes the dfg only if it doesn't exist so far
 * The content of the dfg is set to zero in case it is initialized.
 */
void GeneTask::initDFG( CommunicatorType comm,
                        std::vector<IndexVector>& decomposition ){
  /*
  // this is the clean version. however requires creation of dfg before each
  // combination step
  for(auto d:decomposition){
    std::cout << d << " ,";
  }
  std::cout << "\n";
  if( dfg_ != NULL )
    delete dfg_;

  dfg_ = new DistributedFullGrid<CombiDataType>( dim_, l_, comm,
      this->getBoundary(), p_, false, decomposition );

  */
  // todo: keep in mind
  // in this version the dfg is only created once. this only works if always exactly
  // the same set of processes is used by gene
  // this will probably not work, when the task is moved to another group.
  if( dfgVector_.size() == 0 ){
    dfgVector_.resize(nspecies_,NULL);
    for(int i=0; i<nspecies_; i++){
      dfgVector_[i] = new DistributedFullGrid<CombiDataType>( dim_, l_, comm,
          this->getBoundary(), p_, false, decomposition );
    }
    setZero();

  }
  //std::cout << "initDFG \n";

}
/**
 * This routine creates a new dfg (and initializes it). This is needed when a task is redistributed or recomputed
 * The content of the dfg is set to zero in case it is initialized.
 */
void GeneTask::initDFG2( CommunicatorType comm,
                        std::vector<IndexVector>& decomposition ){
  // this is the clean version. however requires creation of dfg before each
  // combination step
  for(auto d:decomposition){
    std::cout << d << " ,";
  }
  std::cout << "\n";
  if(dfgVector_.size() != nspecies_){
    dfgVector_.resize(nspecies_,NULL);
  }
  for(int i=0; i<nspecies_; i++){
    if(dfgVector_[i] != NULL ){
      delete dfgVector_[i];
    }
    dfgVector_[i] = new DistributedFullGrid<CombiDataType>( dim_, l_, comm,
        this->getBoundary(), p_, false, decomposition );
  }
  setZero();
/*
  // todo: keep in mind
  // in this version the dfg is only created once. this only works if always exactly
  // the same set of processes is used by gene
  // this will probably not work, when the task is moved to another group.
  if( dfg_ == NULL ){
    dfg_ = new DistributedFullGrid<CombiDataType>( dim_, l_, comm,
        this->getBoundary(), p_, false, decomposition );
  }
  //std::cout << "initDFG \n";
   */
}


/** copy and convert data from local checkpoint to dfg
 * This is used for importing data from GENE to our framework
 * */
void GeneTask::setDFG(){
  // step one: copy data of localcheckpoint into dfg

  const std::vector<size_t>& b = checkpoint_.getBounds();
  IndexVector lowerBounds = { b[0], b[2], b[4], b[6], b[8], b[10] };
  IndexVector upperBounds = { b[1], b[3], b[5], b[7], b[9], b[11] };
  IndexVector lcpSizes = upperBounds - lowerBounds;

  MultiArrayRef<GeneComplex,6> lcpData =
    createMultiArrayRef<GeneComplex,6>( checkpoint_.getData(), lcpSizes );

  for(int species=0; species<nspecies_; species++){

    MultiArrayRef6 dfgData =
        createMultiArrayRef<CombiDataType,6>( *dfgVector_[species] );

    // so far, parallelization in x direction is not supported. since this is
    // the innermost dimension, this is not sensible anyway.
    // parallelization in this dimension would require very expensive communication
    // in order to change the ordering in this dimension
    // Parallelization in x required for global cases
    if(!_GENE_Global){
      assert( dfgVector_[species]->getParallelization()[0] == 1 );
    }
    // some checks
    const IndexVector p( dfgVector_[species]->getParallelization().rbegin(),
                          dfgVector_[species]->getParallelization().rend() );
    IndexVector tmp( dfgVector_[species]->getDimension() );
    dfgVector_[species]->getPartitionCoords( tmp );

    // get partition coords of the current process in the ordering as
    // the shape of the multi arrays, i.e. spec, w, v, z, y, x
    IndexVector coords( tmp.rbegin(), tmp.rend() );

    const size_t* dfgShape = dfgData.shape();
    const size_t* lcpShape = lcpData.shape();

    for( DimType d=0; d<dfgVector_[species]->getDimension(); ++d ){
      // for the last rank in the dimension w(d=1), v(2), z(3), x(5) the number of elements in
      // dfg and lcp differ
      if( coords[d] == p[d] - 1 && ( d==1 || d == 2 || d == 3 || d==5 ) ){
        assert( dfgShape[d] == lcpShape[d] + 1 );
      } else{
        if(d==0){
          assert(lcpShape[0] == nspecies_ && dfgShape[0] == 1); //dfg has always 1 coordinate in species direction
        }
        else{
          assert( dfgShape[d] == lcpShape[d] );
        }
      }
    }

    // copy data from local checkpoint to dfg
    // note that on the last process in some dimensions dfg is larger than the
    // local checkpoint
    for( size_t m=0; m < lcpShape[1]; ++m ){ //w
      for( size_t l=0; l < lcpShape[2]; ++l ){ //v
        for( size_t k=0; k < lcpShape[3]; ++k ){ //z
          for( size_t j=0; j < lcpShape[4]; ++j ){ //y
            for( size_t i=0; i < lcpShape[5]; ++i ){ //x
              dfgData[0][m][l][k][j][i].real( lcpData[species][m][l][k][j][i].r );
              dfgData[0][m][l][k][j][i].imag( lcpData[species][m][l][k][j][i].i );
            }
          }
        }
      }
    }



  /*
  // check if last grid point of x is zero
  if( coords[1] == p[1] - 1 ){
    for( size_t n=0; n < dfgShape[0]; ++n ) //n_spec
      for( size_t m=0; m < dfgShape[1]; ++m ) //w
        for( size_t l=0; l < dfgShape[2]; ++l ) //v
          for( size_t k=0; k < dfgShape[3]; ++k ) //z
            for( size_t j=0; j < dfgShape[4]; ++j ){ //y
                assert( dfgData[n][ m ][l][k][j][ dfgShape[5]-1 ]
                        == complex(0.0) );
              }
  }

  // check if last grid point of w is zero
  if( coords[1] == p[1] - 1 ){
    for( size_t n=0; n < dfgShape[0]; ++n ) //n_spec
        for( size_t l=0; l < dfgShape[2]; ++l ) //v
          for( size_t k=0; k < dfgShape[3]; ++k ) //z
            for( size_t j=0; j < dfgShape[4]; ++j ) //y
              for( size_t i=0; i < dfgShape[5]; ++i ){ //x
                assert( dfgData[n][ dfgShape[1]-1 ][l][k][j][i]
                        == complex(0.0) );
              }
  }

  // check if last grid point of v is zero
  if( coords[2] == p[2] - 1 ){
    for( size_t n=0; n < dfgShape[0]; ++n ) //n_spec
      for( size_t m=0; m < dfgShape[1]; ++m ) //w
          for( size_t k=0; k < dfgShape[3]; ++k ) //z
            for( size_t j=0; j < dfgShape[4]; ++j ) //y
              for( size_t i=0; i < dfgShape[5]; ++i ){ //x
                assert( dfgData[n][ m ][ dfgShape[2]-1 ][k][j][i]
                        == complex(0.0) );
              }
  }
  */

  // adapt z-boundary
  adaptBoundaryZ(species);

  // normalize
  normalizeDFG(species);
 }
}

/*
 * This routine applies the z boundary condition
 */
void GeneTask::adaptBoundaryZ(int species){
  //different implementation whether z is parallelized
  //local <-> no parallelization; global <-> parallelized in z
  if( dfgVector_[species]->getParallelization()[2] > 1 )
    adaptBoundaryZglobal(species);
  else
    adaptBoundaryZlocal(species);
}


/** the term local may be misleading in the context of GENE
 *  this simply means that the copy operations are done locally without any
 *  MPI communication. this is only possible in the case that there is only
 *  one process in z direction (and in x, what we require anyway)
 */
void GeneTask::adaptBoundaryZlocal(int species){
  // make sure no parallelization in z and x
  assert( dfgVector_[species]->getParallelization()[0] == 1 );
  assert( dfgVector_[species]->getParallelization()[2] == 1 );

  MultiArrayRef6 dfgData = createMultiArrayRef<CombiDataType,6>( *dfgVector_[species] );

  adaptBoundaryZKernel(dfgData,dfgData,species);
}
/**
 * This method updates the last z component of target data according to the z boundary conditions.
 * The value at z[0] is taken from source data which might be a separate array in case of domain decomposition in z direction
 * @param sourceData contains the z[0] elements at z[0]
 * @param targetData needs to be updated at last z coordinate using boundary condition
 *
 */
void GeneTask::adaptBoundaryZKernel(MultiArrayRef6& sourceData, MultiArrayRef6& targetData, int species){
  // calculate offset and factor
    const IndexVector p( dfgVector_[species]->getParallelization().rbegin(),
                          dfgVector_[species]->getParallelization().rend() );
    IndexVector tmp( dfgVector_[species]->getDimension() );
    dfgVector_[species]->getPartitionCoords( tmp );

    // get partition coords of the current process in the ordering as
    // the shape of the multi arrays, i.e. spec, w, v, z, y, x
    IndexVector coords( tmp.rbegin(), tmp.rend() );
    bool xBorder = coords[5] == p[5] - 1; //check if we are at the upper x border
    IndexType xoffset;
    CombiDataType factor;
    getOffsetAndFactor( xoffset, factor );
    std::cout <<"lx: " <<lx_ <<  "s: " << shat_ << " kymin: " << kymin_ << "\n";

    MASTER_EXCLUSIVE_SECTION{
      std::cout << "xoffset: " << xoffset
                << "factor: " << factor << std::endl;
    }

    const size_t* targetShape = targetData.shape();

    // we ignore the last point in x direction
    size_t nkx; //number of points in x direction in local dfg
    if(xBorder){
      nkx= targetShape[5]-1; //toDo what happens when x not fourier
    }
    else{
      nkx = targetShape[5];
    }
    size_t nkxGlobal = dfgVector_[species]->getGlobalSizes()[0] - 1; //here x is at position 0
    assert(nkxGlobal >= nkx);
    std::cout << "local number of x points: " << nkx << " global number of x points: " << nkxGlobal << "\n";
    // make sure this value is even (actually should be power of two)
    // because we assume the highest kx is always zero
    assert( nkx%2 == 0 ); //toDo global case

    for( size_t n=0; n < targetShape[0]; ++n ){ //n_spec
      for( size_t m=0; m < targetShape[1]; ++m ){ //w
        for( size_t l=0; l < targetShape[2]; ++l ){ //v
          for( size_t j=0; j < targetShape[4]; ++j ){ //y
            for( size_t i=0; i < nkx; ++i ){ //x
              // ignore highest mode, because it is always zero toDo is this only valid in local case?
              if(!_GENE_Global){ //toDo is this right?
                if( i == nkxGlobal/2 )
                  continue;
              }
              if(!_GENE_Global){
                getOffsetAndFactor( xoffset, factor,j+1 );
              }
              else{
                getOffsetAndFactor( xoffset, factor,j+1, (i*1.0)/(nkxGlobal*1.0) * lx_ );
              }
              // calc kx_star
              IndexType kx_star = (i + xoffset)%nkxGlobal; //it might be problematic if kx_star is not on the same process

              // if kx* is the highest mode, increase by one
              //if( kx_star == nkx/2 )
              //  kx_star += 1;

              // this should never happen
              assert(xoffset >= 0);
              assert( kx_star < nkx);

              targetData[n][m][l][ targetShape[3]-1 ][j][i] =
                  sourceData[n][m][l][0][j][ kx_star ] * factor;
            }
          }
        }
      }
    }
}

/**
 * This function calculates the offset in x coordinate and the scaling factor for the z-boundary condition.
 * See diss of Christoph Kowitz for more information.
 * @param xoffset variable that stores the xoffset after call
 * @param factor variable that stores the scaling factors after call
 * @param l defines mode ky=kymin_*l. l is always integer valued
 */
void GeneTask::getOffsetAndFactor( IndexType& xoffset, CombiDataType& factor, IndexType l, real x ){
  // calculate x offset and factor
  if(!_GENE_Global){
  int N = int( round( shat_ * kymin_ * lx_ ) );
  int ky0_ind = 1;
  assert( N == 1);

  xoffset = l*N;
  factor = std::pow(-1.0,N * l);

  // i think this function is only right if nky=1
  assert( l_[1] == 1 && boundary_[1] == false );
  }
  else{
    xoffset = 0;
    double ky = kymin_ * l;
    double pi = boost::math::constants::pi<double>();
    double angle = 2*pi*ky*(x0_ + shat_ *(x - x0_));
    factor = complex(cos(angle),sin(angle));
  }
}


/* name may be misleading, this is for local iv compuations
 * global refers to a global adaption of the z-boundary conditions, which is
 * the case when there is more than one process in this direction
 */
void GeneTask::adaptBoundaryZglobal(int species){
  // make sure at least two processes in z-direction
  assert( dfgVector_[species]->getParallelization()[2] > 1 );
  //make sure that no parallelization in x is performed
  assert(dfgVector_[species]->getParallelization()[0] == 1);
  // everything in normal dfg notation here except MPI subarrays!

  // get parallelization and coordinates of process
  const IndexVector& p = dfgVector_[species]->getParallelization();
  IndexVector coords( dfgVector_[species]->getDimension() );
  dfgVector_[species]->getPartitionCoords( coords );

  // x range of current process in global coordinates
  IndexVector lBounds( dfgVector_[species]->getLowerBounds() );
  IndexVector uBounds( dfgVector_[species]->getUpperBounds() );

  // todo: works only for nky = 1. for linear simulations we do not need this
  // for non-linear simulations the boundary treatment is different
  // for global simulations it is probably different too
  assert( dfgVector_[species]->getGlobalSizes()[1] == 1 );

  IndexType xoffset;
  CombiDataType factor;

  bool sendingProc = false;
  if( coords[2] == 0 ){
    sendingProc = true;
  }
  else if( coords[2] == p[2]-1 ){
    sendingProc = false;
  } else{
    // other processes don't have to do anything here
    return;
  }
/*
  // for each remote process there may be up to two blocks of data to send
  // or receive
  std::vector< std::vector<IndexType> > transferIndicesBlock1( dfg_->getCommunicatorSize() );
  std::vector< std::vector<IndexType> > transferIndicesBlock2( dfg_->getCommunicatorSize() );

  // we ignore the last point in x direction
  IndexType numx = dfg_->getGlobalSizes()[0] - 1;
  assert( xoffset >= 0 && xoffset <= numx ); // otherwise index calculations may not work

  for( IndexType xi=lBounds[0]; xi<uBounds[0]; ++xi ){
    // ignore last global index of x
    if( xi >= numx )
      continue;

    // create a valid global index vector which has xi
    IndexVector gix ( lBounds );

    int whichBlock = 0;

    if( sendingProc ){
      gix[0] = xi - xoffset;

      if( gix[0] < 0 ){
        whichBlock = 2;
        gix[0] += numx; // avoids unclear definition of modulo operation for negative numbers
      } else{
        whichBlock = 1;
      }

      gix[2] = dfg_->getGlobalSizes()[2] - 1;
    }
    else{
      gix[0] = xi + xoffset;

      if( gix[0] >= numx ){
        whichBlock = 2;
        gix[0] = ( xi + xoffset )%numx;
      } else{
        whichBlock = 1;
      }

      gix[2] = 0;
    }

    // check if index out of bounds
    assert( gix[0] >= 0 && gix[0] < numx );

    // find process which has this point (receiving process)
    IndexVector scoords( dfg_->getDimension() );
    dfg_->getPartitionCoords( gix, scoords );
    RankType r = dfg_->getRank( scoords );

    // append xi to process' list of x coordinates
    assert( whichBlock == 1 || whichBlock == 2 );

    if( whichBlock == 1 )
      transferIndicesBlock1[r].push_back( xi );
    else if( whichBlock == 2 )
      transferIndicesBlock2[r].push_back( xi );
  }


  std::vector< MPI_Request > requests;

  // first exchange block1 then exchange block2
  for( size_t block=1; block<=2; ++block ){
    std::vector< std::vector<IndexType> >& transferIndices =
        (block==1) ? transferIndicesBlock1 : transferIndicesBlock2;

    for( RankType r=0; r<dfg_->getCommunicatorSize(); ++r ){
      if( transferIndices[r].size() == 0 )
        continue;

      // set lower bounds of subarray
      IndexVector subarrayLowerBounds = dfg_->getLowerBounds();
      subarrayLowerBounds[0] = transferIndices[r][0];

      // set upper bounds of subarray
      IndexVector subarrayUpperBounds = dfg_->getUpperBounds();
      subarrayUpperBounds[0] = transferIndices[r].back() + 1;

      if( sendingProc ){
        subarrayLowerBounds[2] = 0;
        subarrayUpperBounds[2] = 1;
      } else{
        subarrayLowerBounds[2] = dfg_->getGlobalSizes()[2] - 1;
        subarrayUpperBounds[2] = dfg_->getGlobalSizes()[2];
      }

      // create MPI datatype
      IndexVector sizes( dfg_->getLocalSizes() );
      IndexVector subsizes = subarrayUpperBounds - subarrayLowerBounds;
      // the starts are local indices
      IndexVector starts = subarrayLowerBounds - dfg_->getLowerBounds();

      // convert to mpi notation
      // also we have to use int as type for the indizes
      std::vector<int> csizes(sizes.rbegin(), sizes.rend());
      std::vector<int> csubsizes(subsizes.rbegin(), subsizes.rend());
      std::vector<int> cstarts(starts.rbegin(), starts.rend());

      // create subarray view on data
      MPI_Datatype mysubarray;
      MPI_Type_create_subarray(static_cast<int>(dfg_->getDimension()),
                               &csizes[0], &csubsizes[0], &cstarts[0],
                               MPI_ORDER_C, dfg_->getMPIDatatype(), &mysubarray);
      MPI_Type_commit(&mysubarray);

      MPI_Request req;

      if( sendingProc ){
        MPI_Isend( dfg_->getData(), 1, mysubarray, r, block, dfg_->getCommunicator(),
                   &req);

      } else{
        MPI_Irecv( dfg_->getData(), 1, mysubarray, r, block, dfg_->getCommunicator(),
                   &req);
      }

      requests.push_back(req);
    }
  }
*/
  if(species==0){
    requestArray_ = new MPI_Request[nspecies_];
    receiveBufferArray_.resize(nspecies_);
  }
  // set lower bounds of subarray
  IndexVector subarrayLowerBounds = dfgVector_[species]->getLowerBounds();

  // set upper bounds of subarray
  IndexVector subarrayUpperBounds = dfgVector_[species]->getUpperBounds();

  if( sendingProc ){
    subarrayLowerBounds[2] = 0;
    subarrayUpperBounds[2] = 1;
  } else{
    subarrayLowerBounds[2] = dfgVector_[species]->getGlobalSizes()[2] - 1;
    subarrayUpperBounds[2] = dfgVector_[species]->getGlobalSizes()[2];
  }

  // create MPI datatype
  IndexVector sizes( dfgVector_[species]->getLocalSizes() );
  IndexVector subsizes = subarrayUpperBounds - subarrayLowerBounds;
  // the starts are local indices
  IndexVector starts = subarrayLowerBounds - dfgVector_[species]->getLowerBounds();

  // convert to mpi notation
  // also we have to use int as type for the indizes
  std::vector<int> csizes(sizes.rbegin(), sizes.rend());
  std::vector<int> csubsizes(subsizes.rbegin(), subsizes.rend());
  std::vector<int> cstarts(starts.rbegin(), starts.rend());

  // create subarray view on data
  MPI_Datatype mysubarray;
  MPI_Type_create_subarray(static_cast<int>(dfgVector_[species]->getDimension()),
                           &csizes[0], &csubsizes[0], &cstarts[0],
                           MPI_ORDER_C, dfgVector_[species]->getMPIDatatype(), &mysubarray);
  MPI_Type_commit(&mysubarray);
  if(sendingProc){
    coords[2] = p[2]-1 ;
  }
  else{
    coords[2] = 0;
  }
  RankType r = dfgVector_[species]->getRank( coords );
  if( sendingProc ){
    MPI_Isend( dfgVector_[species]->getData(), 1, mysubarray, r, 1000, dfgVector_[species]->getCommunicator(),
               &requestArray_[species]);

  } else{
    IndexType numElements = 1;
    for(int i=0; i<subsizes.size(); i++){
      numElements *= subsizes[i];
    }
    receiveBufferArray_[species] = new CombiDataType[numElements];
    MPI_Irecv( receiveBufferArray_[species], 1, mysubarray, r, 1000, dfgVector_[species]->getCommunicator(),
               &requestArray_[species]);
  }



  if(species == nspecies_ - 1){
    MPI_Waitall(nspecies_, requestArray_, MPI_STATUSES_IGNORE );
    //toDo check if correct
    // nur letze processe in z-richtung
    if(!sendingProc){
      for(int i=0; i<nspecies_; i++){
        MultiArrayRef6 dfgData = createMultiArrayRef<CombiDataType,6>( *dfgVector_[i] );
        MultiArrayRef6 receivedData =  MultiArrayRef<CombiDataType,6>(receiveBufferArray_[i],subsizes);
        adaptBoundaryZKernel(receivedData,dfgData,i);
        delete receiveBufferArray_[i];
      }
    }
    delete[] requestArray_;
    receiveBufferArray_.clear();
  }
}


/* copy dfg data to local checkpoint
 * This is used to get the data from our framework to GENE
 * */
void GeneTask::getDFG(){
  // get multiarrayref to lcp
  // todo: put this into a function as it is used multiple times
  const std::vector<size_t>& b = checkpoint_.getBounds();
  IndexVector lowerBounds = { b[0], b[2], b[4], b[6], b[8], b[10] };
  IndexVector upperBounds = { b[1], b[3], b[5], b[7], b[9], b[11] };
  IndexVector lcpSizes = upperBounds - lowerBounds;
  MultiArrayRef<GeneComplex,6> lcpData =
    createMultiArrayRef<GeneComplex,6>( checkpoint_.getData(), lcpSizes );



  // copy data back to lcp
  // note that the last grid points in x,z,v,w dimension are ignored
  const size_t* lcpShape = lcpData.shape();
  for( size_t n=0; n < lcpShape[0]; ++n ){ //n_spec
    // get multiarrayref to dfg
    MultiArrayRef6 dfgDataN =
    createMultiArrayRef<CombiDataType,6>( *dfgVector_[n] );
    for( size_t m=0; m < lcpShape[1]; ++m ){ //w
      for( size_t l=0; l < lcpShape[2]; ++l ){ //v
        for( size_t k=0; k < lcpShape[3]; ++k ){ //z
          for( size_t j=0; j < lcpShape[4]; ++j ){ //y
            for( size_t i=0; i < lcpShape[5]; ++i ){ //x
              lcpData[n][m][l][k][j][i].r = dfgDataN[0][m][l][k][j][i].real();
              lcpData[n][m][l][k][j][i].i = dfgDataN[0][m][l][k][j][i].imag();
            }
          }
        }
      }
    }
  }
}


/* we normalize the magnitude by the square root of nrg0 (first entry in nrg.dat)
 * we also normalize the average phase to zero
 */
void GeneTask::normalizeDFG(int species){
  const bool normalizePhase = false;
  const bool normalizeAmplitude = false;

  // 1: normalize phase

  if( normalizePhase ){
    // compute local mean value of dfg
    CombiDataType localMean(0.0);
    std::vector<CombiDataType>& data = dfgVector_[species]->getElementVector();
    for( size_t i=0; i<data.size(); ++i )
      localMean += data[i];

    // allreduce mean to get global mean value
    CombiDataType globalMean(0.0);
    MPI_Allreduce( &localMean, &globalMean, 1,
                   dfgVector_[species]->getMPIDatatype(), MPI_SUM,
                   theMPISystem()->getLocalComm() );


    // subtract phase of local mean from all values of dfg
    real phaseShift = std::arg(globalMean);
    for( size_t i=0; i<data.size(); ++i )
      data[i] *= std::exp( - complex(0,1) * phaseShift );

    MASTER_EXCLUSIVE_SECTION{
      std::cout << "phase shift = " << phaseShift << std::endl;
    }
  }


  // 2: normalize magnitude

  if(normalizeAmplitude){
    // get squared value of nrg0 (only available on master process)
    real factor = 1.0 / std::sqrt(nrg_);

    // broadcast to other procs
    MPI_Bcast( &factor
        , 1, MPI_DOUBLE,
               theMPISystem()->getMasterRank(),
               theMPISystem()->getLocalComm() );

    // divide values of dfg
    std::vector<CombiDataType>& data = dfgVector_[species]->getElementVector();
    for( size_t i=0; i<data.size(); ++i )
      data[i] *= factor;
  }

}

/*
inline bool GeneTask::failNow( const int& globalRank ){
  FaultsInfo faultsInfo = faultsInfo_;
  IndexVector iF = faultsInfo_.iterationFaults_;
  IndexVector rF = faultsInfo_.globalRankFaults_;

  std::vector<IndexType>::iterator it;
  it = std::find(iF.begin(), iF.end(), combiStep_);
  IndexType idx = std::distance(iF.begin(),it);
  //std::cout << "faultInfo" << iF[0] << " " << rF[0] << "\n";
  // Check if current iteration is in iterationFaults_
  while (it!=iF.end()){
    // Check if my rank is the one that fails
    if (globalRank == rF[idx])
      return true;
    it = std::find(++it, iF.end(), combiStep_);
    idx = std::distance(iF.begin(),it);
  }
  return false;
}
*/
} /* namespace combigrid */
