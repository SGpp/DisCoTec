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

namespace combigrid
{

GeneTask::GeneTask( DimType dim, LevelVector& l,
                    std::vector<bool>& boundary, real coeff, LoadModel* loadModel,
                    std::string& path, real dt, size_t nsteps )
    : Task( dim, l, boundary, coeff, loadModel), path_( path ),  dt_( dt ),
      nsteps_( nsteps )
{

// theres only one boundary configuration allowed at the moment
assert( boundary[0] == true );//w
assert( boundary[1] == false );//v
assert( boundary[2] == true );//z
assert( boundary[3] == false );//y
assert( boundary[4] == true );//x
}

GeneTask::GeneTask()
{
}

GeneTask::~GeneTask()
{
// TODO Auto-generated destructor stub
}

void
GeneTask::run( CommunicatorType& lcomm )
{
// change dir to wdir
if( chdir( path_.c_str() ) ){


  printf( "could not change to directory %s \n", path_.c_str() );
  MPI_Abort( MPI_COMM_WORLD, 1 );
}

char cwd[1024];
getcwd(cwd, sizeof(cwd));

std::cout << "changed dir to " << cwd << std::endl;

}


void GeneTask::init(CommunicatorType& lcomm){

}

void
GeneTask::writeLocalCheckpoint( GeneComplex* data, size_t size,
                                std::vector<size_t>& sizes,
                                std::vector<size_t>& bounds )
{
// todo: doing it like this will require two times copying
checkpoint_.writeCheckpoint( data, size, sizes, bounds );
}


/*
 * todo: having to pass the communicator here is rather inconvenient. as each
 * task is in some sense attached to a particular process group, this should
 * be represented in the design as well */
void
GeneTask:: getFullGrid( FullGrid<complex>& fg, RankType localRootID,
                        CommunicatorType& lcomm)
{
// check whether checkpoint exists
assert( checkpoint_.getSize() > 0 );

int lsize;
RankType localID;
MPI_Comm_size( lcomm, &lsize );
MPI_Comm_rank( lcomm, &localID );

// start time gather global cp
double tstart_gathercp = MPI_Wtime();

// collect checkpoint on local root
std::vector<GeneComplex> recvdata;
std::vector<int> sendbounds;
std::vector<int> recvbounds;
const GeneLocalCheckpoint& lcp = this->getLocalCheckpoint();

assert( lcp.getSize() > 0 );

// gather data
GeneComplex* senddata = const_cast<GeneComplex*>( lcp.getData() );

if( localID == localRootID )
  recvdata.resize( static_cast<size_t>( lsize ) * lcp.getSize() );

MPI_Gather( senddata, static_cast<int>( lcp.getSize() ), MPI_DOUBLE_COMPLEX,
            &recvdata[0], static_cast<int>( lcp.getSize() ),
            MPI_DOUBLE_COMPLEX,
            localRootID, lcomm );

// gather bounds
sendbounds.assign( lcp.getBounds().begin(), lcp.getBounds().end() );

if( localID == localRootID )
  recvbounds.resize( lsize * sendbounds.size() );

MPI_Gather( &sendbounds[0], static_cast<int>( sendbounds.size() ), MPI_INT,
            &recvbounds[0], static_cast<int>( sendbounds.size() ), MPI_INT,
            localRootID, lcomm );

if( localID == localRootID ){
  // start time convertCPtoFG
  double tstart_converCPtoFG = MPI_Wtime();

  const std::vector<size_t>& globalSize = this->getLocalCheckpoint().getSizes();

  // create GeneGrid
  GeneGrid global(
      boost::extents[globalSize[0]][globalSize[1]][globalSize[2]][globalSize[3]][globalSize[4]][globalSize[5]] );

  // fill global grid with checkpoint data
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

    // fill global array with data
    // set view to global array
    typedef GeneGrid::index_range range;
    GeneGrid::array_view<6>::type lview = global[boost::indices[range(
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
              for( index i = 0; i < lni0; ++i )
                lview[n][m][l][k][j][i] = local[n][m][l][k][j][i];
  }

  assert( globalSize[0] == 1 ); // only one species at the moment
  assert( globalSize[4] == 1 ); // will only work for this at the moment
  assert( globalSize[5] % 2 == 1 ); //nx0==ni0 must be odd

  fg.createFullGrid();
  CombiGeneConverter::GeneGridToFullGrid( global, fg );
}
}


DistributedFullGrid<complex>& GeneTask::getDistributedFullGrid(){
  //todo: implement!
  assert( false && "not implemented" );
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
const GeneLocalCheckpoint& lcp = this->getLocalCheckpoint();

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

GeneComplex* ldata = const_cast<GeneComplex*>( lcp.getData() );
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

} /* namespace combigrid */
