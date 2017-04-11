/*
 * GeneTask.cpp
 *
 *  Created on: Jul 10, 2014
 *      Author: heenemo
 */

#include "GeneTask.hpp"
#include <sstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <unistd.h> //chdir
#include <fstream>
#include "CombiGeneConverter.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/fullgrid/MultiArray.hpp"
//#include "sgpp/distributedcombigrid/utils/StatsContainer.hpp"

namespace combigrid
{

GeneTask::GeneTask( DimType dim, LevelVector& l,
                    std::vector<bool>& boundary, real coeff, LoadModel* loadModel,
                    std::string& path, real dt, size_t nsteps,
                    real shat, real kymin, real lx, int ky0_ind,
                    IndexVector p , FaultCriterion *faultCrit )
    : Task( dim, l, boundary, coeff, loadModel),
      path_( path ),
      dt_( dt ),
      nsteps_( nsteps ),
      stepsTotal_(0),
      combiStep_(0),
      shat_( shat ),
      kymin_( kymin ),
      lx_( lx ),
      ky0_ind_( ky0_ind ),
      p_(p),
      checkpoint_(), initialized_(false),
      checkpointInitialized_(false),
      faultCriterion_(faultCrit)
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
    dfg_(NULL),
    nrg_(0.0),
    initialized_(false),
    checkpointInitialized_(false), faultCriterion_((new FaultCriterion()))
{
  ;
}

GeneTask::~GeneTask()
{
// TODO Auto-generated destructor stub
}

void
GeneTask::run( CommunicatorType lcomm )
{
  using namespace std::chrono;

  // change dir to wdir
  if( chdir( path_.c_str() ) ){


    printf( "could not change to directory %s \n", path_.c_str() );
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }

  char cwd[1024];
  getcwd(cwd, sizeof(cwd));

  MASTER_EXCLUSIVE_SECTION{
    std::cout << "run task " << this->getID() << std::endl;
  }
  int globalRank;
  // MPI_Comm_rank(lcomm, &lrank);
  MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);
  if(combiStep_ != 0){
    //theStatsContainer()->setTimerStop("computeIterationRank" + std::to_string(globalRank));
    //theStatsContainer()->setValue("computeIterationRank" + std::to_string(globalRank),0.0);
  }
  startTimeIteration_ = high_resolution_clock::now();

  //theStatsContainer()->setTimerStart("computeIterationRank" + std::to_string(globalRank));

  //printf("running gene!!! \n");

}

void GeneTask::decideToKill(){
  using namespace std::chrono;

  int globalRank;
  // MPI_Comm_rank(lcomm, &lrank);
  MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);
  //theStatsContainer()->setTimerStop("computeIterationRank" + std::to_string(globalRank));
  duration<real> dur = high_resolution_clock::now() - startTimeIteration_;
  real t_iter = dur.count();
  //std::cout << "Current iteration took " << t_iter << "\n";

  //theStatsContainer()->setTimerStart("computeIterationRank" + std::to_string(globalRank));


  //check if killing necessary
  //std::cout << "failNow result " << failNow(globalRank) << " at rank: " << globalRank <<" at step " << combiStep_ << "\n" ;
  //real t = dt_ * nsteps_ * combiStep_;
  if (faultCriterion_->failNow(combiStep_, t_iter, globalRank)){
        std::cout<<"Rank "<< globalRank <<" failed at iteration "<<combiStep_<<std::endl;
        simft::Sim_FT_kill_me();
  }
  combiStep_++;
}
void GeneTask::init(CommunicatorType lcomm, std::vector<IndexVector> decomposition){
  if( dfg_ == NULL ){
      dfg_ = new DistributedFullGrid<CombiDataType>( dim_, l_, lcomm,
          this->getBoundary(), p_, false);
  }
  initialized_ = true;
}

void
GeneTask::writeLocalCheckpoint( GeneComplex* data, size_t size,
                                std::vector<size_t>& sizes,
                                std::vector<size_t>& bounds )
{
  // todo: doing it like this will require two times copying
  checkpoint_.writeCheckpoint( data, size, sizes, bounds );
  checkpointInitialized_= true;

}

void GeneTask::InitLocalCheckpoint(size_t size,
    std::vector<size_t>& sizes,
    std::vector<size_t>& bounds ){
  checkpoint_.initCheckpoint(size,sizes,bounds);
  checkpointInitialized_= true;
}


void GeneTask::getFullGrid( FullGrid<CombiDataType>& fg, RankType lroot,
                            CommunicatorType lcomm)
{
  dfg_->gatherFullGrid(fg, lroot);
}


DistributedFullGrid<complex>& GeneTask::getDistributedFullGrid(){
  return *dfg_;
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


void GeneTask::setZero(){
  std::vector<CombiDataType>& data = dfg_->getElementVector();

  for( size_t i=0; i<data.size(); ++i )
    data[i] = std::complex<real>(0);
}


void GeneTask::initDFG( CommunicatorType comm,
                        std::vector<IndexVector>& decomposition ){
  // this is the clean version. however requires creation of dfg before each
  // combination step
/*
  if( dfg_ != NULL )
    delete dfg_;

  dfg_ = new DistributedFullGrid<CombiDataType>( dim_, l_, comm,
      this->getBoundary(), p_, false, decomposition );
*/

  // todo: keep in mind
  // in this version the dfg is only created once. this only works if always exactly
  // the same set of processes is used by gene
  // this will probably not work, when the task is moved to another group.
  if( dfg_ == NULL ){
    dfg_ = new DistributedFullGrid<CombiDataType>( dim_, l_, comm,
        this->getBoundary(), p_, false, decomposition );
  }
  //std::cout << "initDFG \n";
}


/** copy and convert data from local checkpoint to dfg */
void GeneTask::setDFG(){
  // step one: copy data of localcheckpoint into dfg

  const std::vector<size_t>& b = checkpoint_.getBounds();
  IndexVector lowerBounds = { b[0], b[2], b[4], b[6], b[8], b[10] };
  IndexVector upperBounds = { b[1], b[3], b[5], b[7], b[9], b[11] };
  IndexVector lcpSizes = upperBounds - lowerBounds;

  MultiArrayRef<GeneComplex,6> lcpData =
    createMultiArrayRef<GeneComplex,6>( checkpoint_.getData(), lcpSizes );

  MultiArrayRef6 dfgData =
      createMultiArrayRef<CombiDataType,6>( *dfg_ );

  // so far, parallelization in x direction is not supported. since this is
  // the innermost dimension, this is not sensible anyway.
  // parallelization in this dimension would require very expensive communication
  // in order to change the ordering in this dimension
  assert( dfg_->getParallelization()[0] == 1 );

  // some checks
  const IndexVector p( dfg_->getParallelization().rbegin(),
                        dfg_->getParallelization().rend() );
  IndexVector tmp( dfg_->getDimension() );
  dfg_->getPartitionCoords( tmp );

  // get partition coords of the current process in the ordering as
  // the shape of the multi arrays, i.e. spec, w, v, z, y, x
  IndexVector coords( tmp.rbegin(), tmp.rend() );

  const size_t* dfgShape = dfgData.shape();
  const size_t* lcpShape = lcpData.shape();

  for( DimType d=0; d<dfg_->getDimension(); ++d ){
    // for the last rank in the dimension w(d=1), v(2), z(3), x(5) the number of elements in
    // dfg and lcp differ
    if( coords[d] == p[d] - 1 && ( d==1 || d == 2 || d == 3 || d==5 ) ){
      assert( dfgShape[d] == lcpShape[d] + 1 );
    } else{
      assert( dfgShape[d] == lcpShape[d] );
    }
  }

  // copy data from local checkpoint to dfg
  // note that on the last process in some dimensions dfg is larger than the
  // local checkpoint
  for( size_t n=0; n < lcpShape[0]; ++n ){ //n_spec
    for( size_t m=0; m < lcpShape[1]; ++m ){ //w
      for( size_t l=0; l < lcpShape[2]; ++l ){ //v
        for( size_t k=0; k < lcpShape[3]; ++k ){ //z
          for( size_t j=0; j < lcpShape[4]; ++j ){ //y
            for( size_t i=0; i < lcpShape[5]; ++i ){ //x
              dfgData[n][m][l][k][j][i].real( lcpData[n][m][l][k][j][i].r );
              dfgData[n][m][l][k][j][i].imag( lcpData[n][m][l][k][j][i].i );
            }
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
  adaptBoundaryZ();

  // normalize
  normalizeDFG();
}


void GeneTask::adaptBoundaryZ(){
  if( dfg_->getParallelization()[2] > 1 )
    adaptBoundaryZglobal();
  else
    adaptBoundaryZlocal();
}


/** the term local may be misleading in the context of GENE
 *  this simply means that the copy operations are done locally without any
 *  MPI communication. this is only possible in the case that there is only
 *  one process in z direction (and in x, what we require anyway)
 */
void GeneTask::adaptBoundaryZlocal(){
  // make sure no parallelization in z and x
  assert( dfg_->getParallelization()[0] == 1 );
  assert( dfg_->getParallelization()[2] == 1 );

  MultiArrayRef6 dfgData = createMultiArrayRef<CombiDataType,6>( *dfg_ );

  // calculate offset and factor
  IndexType xoffset;
  CombiDataType factor;
  getOffsetAndFactor( xoffset, factor );

  MASTER_EXCLUSIVE_SECTION{
    std::cout << "xoffset: " << xoffset
              << "factor: " << factor << std::endl;
  }

  const size_t* dfgShape = dfgData.shape();

  // we ignore the last point in x direction
  size_t nkx = dfgShape[5]-1;

  // make sure this value is even (actually should be power of two)
  // because we assume the highest kx is always zero
  assert( nkx%2 == 0 );

  for( size_t n=0; n < dfgShape[0]; ++n ) //n_spec
    for( size_t m=0; m < dfgShape[1]; ++m ) //w
      for( size_t l=0; l < dfgShape[2]; ++l ) //v
          for( size_t j=0; j < dfgShape[4]; ++j ) //y
            for( size_t i=0; i < nkx; ++i ){ //x
              // ignore highest mode, because it is always zero
              if( i == nkx/2 )
                continue;

              // calc kx_star
              IndexType kx_star = (i + xoffset)%nkx;

              // if kx* is the highest mode, increase by one
              //if( kx_star == nkx/2 )
              //  kx_star += 1;

              // this should never happen
              assert( kx_star < nkx);

              dfgData[n][m][l][ dfgShape[3]-1 ][j][i] =
                  dfgData[n][m][l][0][j][ kx_star ] * factor;
            }
}


void GeneTask::getOffsetAndFactor( IndexType& xoffset, CombiDataType& factor ){
  // calculate x offset and factor
  int N = int( round( shat_ * kymin_ * lx_ ) );
  int ky0_ind = 1;
  assert( N == 1);

  xoffset = N;
  factor = std::pow(-1.0,N);

  // i think this function is only right if nky=1
  assert( l_[1] == 1 && boundary_[1] == false );
}


/* name may be misleading, this is for local iv compuations
 * global refers to a global adaption of the z-boundary conditions, which is
 * the case when there is more than one process in this direction
 */
void GeneTask::adaptBoundaryZglobal(){
  // make sure at least two processes in z-direction
  assert( dfg_->getParallelization()[2] > 1 );

  // everything in normal dfg notation here except MPI subarrays!

  // get parallelization and coordinates of process
  const IndexVector& p = dfg_->getParallelization();
  IndexVector coords( dfg_->getDimension() );
  dfg_->getPartitionCoords( coords );

  // x range of current process in global coordinates
  IndexVector lBounds( dfg_->getLowerBounds() );
  IndexVector uBounds( dfg_->getUpperBounds() );

  // todo: works only for nky = 1. for linear simulations we do not need this
  // for non-linear simulations the boundary treatment is different
  // for global simulations it is probably different too
  assert( dfg_->getGlobalSizes()[1] == 1 );

  IndexType xoffset;
  CombiDataType factor;
  getOffsetAndFactor( xoffset, factor );

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

  MPI_Waitall( requests.size(), &requests[0], MPI_STATUSES_IGNORE );

  // todo: multiply with factor
  assert( false && "no factor correction" );

  // nur letze processe in z-richtung

  // loop über x,y,v,w
  // multipliziere
}


/* copy dfg data to local checkpoint */
void GeneTask::getDFG(){
  // get multiarrayref to lcp
  // todo: put this into a function as it is used multiple times
  const std::vector<size_t>& b = checkpoint_.getBounds();
  IndexVector lowerBounds = { b[0], b[2], b[4], b[6], b[8], b[10] };
  IndexVector upperBounds = { b[1], b[3], b[5], b[7], b[9], b[11] };
  IndexVector lcpSizes = upperBounds - lowerBounds;
  MultiArrayRef<GeneComplex,6> lcpData =
    createMultiArrayRef<GeneComplex,6>( checkpoint_.getData(), lcpSizes );

  // get multiarrayref to dfg
  MultiArrayRef6 dfgData =
      createMultiArrayRef<CombiDataType,6>( *dfg_ );

  // copy data back to lcp
  // note that the last grid points in x,z,v,w dimension are ignored
  const size_t* lcpShape = lcpData.shape();
  for( size_t n=0; n < lcpShape[0]; ++n ){ //n_spec
    for( size_t m=0; m < lcpShape[1]; ++m ){ //w
      for( size_t l=0; l < lcpShape[2]; ++l ){ //v
        for( size_t k=0; k < lcpShape[3]; ++k ){ //z
          for( size_t j=0; j < lcpShape[4]; ++j ){ //y
            for( size_t i=0; i < lcpShape[5]; ++i ){ //x
              lcpData[n][m][l][k][j][i].r = dfgData[n][m][l][k][j][i].real();
              lcpData[n][m][l][k][j][i].i = dfgData[n][m][l][k][j][i].imag();
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
void GeneTask::normalizeDFG(){
  const bool normalizePhase = false;
  const bool normalizeAmplitude = false;

  // 1: normalize phase

  if( normalizePhase ){
    // compute local mean value of dfg
    CombiDataType localMean(0.0);
    std::vector<CombiDataType>& data = dfg_->getElementVector();
    for( size_t i=0; i<data.size(); ++i )
      localMean += data[i];

    // allreduce mean to get global mean value
    CombiDataType globalMean(0.0);
    MPI_Allreduce( &localMean, &globalMean, 1,
                   dfg_->getMPIDatatype(), MPI_SUM,
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
    std::vector<CombiDataType>& data = dfg_->getElementVector();
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
