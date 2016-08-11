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

namespace combigrid
{

GeneTask::GeneTask( DimType dim, LevelVector& l,
                    std::vector<bool>& boundary, real coeff, LoadModel* loadModel,
                    std::string& path, real dt, size_t nsteps,
                    real shat, real kymin, real lx, int ky0_ind,
                    IndexVector p )
    : Task( dim, l, boundary, coeff, loadModel),
      path_( path ),
      dt_( dt ),
      nsteps_( nsteps ),
      shat_( shat ),
      kymin_( kymin ),
      lx_( lx ),
      ky0_ind_( ky0_ind ),
      p_(p)
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
    dfg_(NULL)
{
}

GeneTask::~GeneTask()
{
// TODO Auto-generated destructor stub
}

void
GeneTask::run( CommunicatorType lcomm )
{
// change dir to wdir
if( chdir( path_.c_str() ) ){


  printf( "could not change to directory %s \n", path_.c_str() );
  MPI_Abort( MPI_COMM_WORLD, 1 );
}

char cwd[1024];
getcwd(cwd, sizeof(cwd));

int rank;
MPI_Comm_rank( MPI_COMM_WORLD, &rank );
std::cout << "rank " << rank << " changed dir to " << cwd << std::endl;

}


void GeneTask::init(CommunicatorType lcomm){

}

void
GeneTask::writeLocalCheckpoint( GeneComplex* data, size_t size,
                                std::vector<size_t>& sizes,
                                std::vector<size_t>& bounds )
{
// todo: doing it like this will require two times copying
checkpoint_.writeCheckpoint( data, size, sizes, bounds );


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
GeneTask::saveCheckpoint( const FullGrid<complex>& fg, const char* filename )
{
assert( fg.isGridCreated() );

size_t nspec = 1;
size_t nw = static_cast<size_t>( std::pow( 2.0, fg.getLevels()[4] ) );
size_t nv = static_cast<size_t>( std::pow( 2.0, fg.getLevels()[3] ) );
size_t nz = static_cast<size_t>( std::pow( 2.0, fg.getLevels()[2] ) );
size_t ny = 1;
size_t nx = static_cast<size_t>( std::pow( 2.0, fg.getLevels()[0] ) ) + 1;

// create temporary global GeneGrid
GeneGrid global_tmp( boost::extents[nspec][nw][nv][nz][ny][nx] );

// convert fg to gene grid
FullGrid<complex>& ncfg = const_cast<FullGrid<complex>&>( fg );
CombiGeneConverter::FullGridToGeneGrid( ncfg, global_tmp );

// save to file
std::ofstream cpFile( filename, std::ofstream::binary );

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
res[5] = 1;
cpFile.write( (char*) res, sizeof(int) * 6 );
std::cout << "resolution:" << "\n" << "\t nx(" << sizeof(int) << ") " << nx
    << "\n" << "\t ny(" << sizeof(int) << ") " << ny << "\n" << "\t nz("
    << sizeof(int) << ") " << nz << "\n" << "\t nv(" << sizeof(int) << ") "
    << nv << "\n" << "\t nw(" << sizeof(int) << ") " << nw << "\n"
    << "\t n_spec(" << sizeof(int) << ") " << nspec << std::endl;

// write data
size_t dsize( global_tmp.num_elements() );
cpFile.write( (char*) global_tmp.origin(), sizeof(complex) * dsize );

cpFile.close();
}


void GeneTask::setZero(){
  // todo: implement
}


void GeneTask::initDFG( CommunicatorType comm,
                        std::vector<IndexVector>& decomposition ){
  if( dfg_ != NULL )
    delete dfg_;

  dfg_ = new DistributedFullGrid<CombiDataType>( dim_, l_, comm,
      this->getBoundary(), p_, decomposition );

  /*
  for( int r=0; r<dfg_->getCommunicatorSize(); ++r ){
    if( r == dfg_->getMpiRank() ){
      std::cout << "rank " << r << "\n"
                << "\t lower bounds " << dfg_->getLowerBounds() << "\n"
                << "\t upper bounds " << dfg_->getUpperBounds()
                << std::endl;
    }

    MPI_Barrier( dfg_->getCommunicator() );
  }*/


}


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

  // compute norm of lcp
  {
	real mynorm = 0.0;
	GeneComplex* data = checkpoint_.getData();
	for( size_t i=0; i<checkpoint_.getSize(); ++i ){
	  CombiDataType tmp( data[i].r, data[i].i );
	  mynorm += std::abs( tmp );
	}

	real norm;
	MPI_Allreduce( &mynorm, &norm, 1, MPI_DOUBLE, MPI_SUM,
				 theMPISystem()->getLocalComm() );

	MASTER_EXCLUSIVE_SECTION{
	  std::cout << "norm before lcp copy " << norm << std::endl;
	}
  }

  {
	  real mynorm = 0.0;

	  for( size_t n=0; n < lcpShape[0]; ++n ){ //n_spec
		for( size_t m=0; m < lcpShape[1]; ++m ){ //w
		  for( size_t l=0; l < lcpShape[2]; ++l ){ //v
			for( size_t k=0; k < lcpShape[3]; ++k ){ //z
			  for( size_t j=0; j < lcpShape[4]; ++j ){ //y
				for( size_t i=0; i < lcpShape[5]; ++i ){ //x
				  CombiDataType tmp;
				  tmp.real( lcpData[n][m][l][k][j][i].r );
				  tmp.imag( lcpData[n][m][l][k][j][i].i );
				  mynorm += std::abs( tmp );
				}
			  }
			}
		  }
		}
	  }

	real norm;
	MPI_Allreduce( &mynorm, &norm, 1, MPI_DOUBLE, MPI_SUM,
			 	 	 theMPISystem()->getLocalComm() );

	MASTER_EXCLUSIVE_SECTION{
		std::cout << "norm of lcpData " << norm << std::endl;
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

  {
	  // calc norm for debugging
	  real mynorm = 0.0;
	  std::vector<CombiDataType>& data = dfg_->getElementVector();
	  for( size_t i=0; i<data.size(); ++i ){
		  mynorm += std::abs( data[i] );
	  }

	  real norm;
	  MPI_Allreduce( &mynorm, &norm, 1, MPI_DOUBLE, MPI_SUM,
					 theMPISystem()->getLocalComm() );

	  MASTER_EXCLUSIVE_SECTION{
		  std::cout << "norm before boundary adaption " << norm << std::endl;
	  }
  }

  // adapt z-boundary
  adaptBoundaryZ();
}


void GeneTask::adaptBoundaryZ(){
  if( dfg_->getParallelization()[2] > 1 )
    adaptBoundaryZglobal();
  else
    adaptBoundaryZlocal();
}



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

  for( size_t n=0; n < dfgShape[0]; ++n ) //n_spec
    for( size_t m=0; m < dfgShape[1]; ++m ) //w
      for( size_t l=0; l < dfgShape[2]; ++l ) //v
          for( size_t j=0; j < dfgShape[4]; ++j ) //y
            for( size_t i=0; i < nkx; ++i ){ //x
                dfgData[n][m][l][ dfgShape[3]-1 ][j][i] =
                    dfgData[n][m][l][0][j][ (i + xoffset)%nkx ];
            }

  // todo: multiply with factor
}


void GeneTask::getOffsetAndFactor( IndexType& xoffset, CombiDataType& factor ){
  // calculate x offset and factor
  int N = int( round( shat_ * kymin_ * lx_ ) );
  int ky0_ind = 1;
  assert( N == 1);

  xoffset = N;
  factor = std::pow(-1.0,N);
}


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

  IndexType numx = dfg_->getGlobalSizes()[0];
  assert( xoffset >= 0 && xoffset <= numx ); // otherwise index calculations may not work

  for( IndexType xi=lBounds[0]; xi<uBounds[0]; ++xi ){
    // create a valid global index vector which has xi
    IndexVector gix ( lBounds );

    if( sendingProc ){
      gix[0] = xi - xoffset;

      if( gix[0] < 0 )
        gix[0] += numx; // avoids unclear definition of modulo operation for negative numbers

      gix[2] = dfg_->getGlobalSizes()[2] - 1;
    }
    else{
      gix[0] = ( xi + xoffset )%numx;
      //gix[0] = xi + xoffset;
      gix[2] = 0;
    }

    //if( !( gix[0] >= 0 && gix[0] < numx ) )
    //  continue;

    // check if index out of bounds
    assert( gix[0] >= 0 && gix[0] < numx );

    // find process which has this point (receiving process)
    IndexVector scoords( dfg_->getDimension() );
    dfg_->getPartitionCoords( gix, scoords );
    RankType r = dfg_->getRank( scoords );

    // append xi to process' list of x coordinates
    transferIndicesBlock1[r].push_back( xi );

    for( RankType rnk =0; rnk < theMPISystem()->getNumProcs(); ++rnk ){
      if( rnk == theMPISystem()->getLocalRank() ){
        if( sendingProc ){
          std::cout << "rank " << theMPISystem()->getLocalRank()
                    << " send " << xi << " to rank" << r << std::endl;
        } else{
          std::cout << "rank " << theMPISystem()->getLocalRank()
                    << " recv " << xi << " from rank" << r << std::endl;
        }
      }
      MPI_Barrier( theMPISystem()->getLocalComm() );
    }
  }




  std::vector< MPI_Request > requests;

  for( RankType r=0; r<dfg_->getCommunicatorSize(); ++r ){
    if( transferIndicesBlock1[r].size() == 0 )
      continue;

    // set lower bounds of subarray
    IndexVector subarrayLowerBounds = dfg_->getLowerBounds();
    subarrayLowerBounds[0] = transferIndicesBlock1[r][0];

    // set upper bounds of subarray
    IndexVector subarrayUpperBounds = dfg_->getUpperBounds();
    subarrayUpperBounds[0] = transferIndicesBlock1[r].back();

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
      MPI_Isend( dfg_->getData(), 1, mysubarray, r, 0, dfg_->getCommunicator(),
                 &req);

    } else{
      MPI_Irecv( dfg_->getData(), 1, mysubarray, r, 0, dfg_->getCommunicator(),
                 &req);
    }

    requests.push_back(req);
  }

  MPI_Waitall( requests.size(), &requests[0], MPI_STATUSES_IGNORE );

  // todo: multiply with factor
}

} /* namespace combigrid */
