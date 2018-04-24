#include <assert.h>
#include <fstream>
#include "sgpp/distributedcombigrid/utils/Config.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/fullgrid/MultiArray.hpp"
#include <iostream>
#include <vector>
#include <sys/stat.h>
#include "boost/lexical_cast.hpp"
#include "timing.h"

using namespace combigrid;

void readCheckpoint( const char* ggFileName,
                     std::vector<CombiDataType>& data,
                     IndexVector& resolution );

void readPlotFile( const char* pltFileName,
                   std::vector<CombiDataType>& data,
                   IndexVector& resolution );

void
calcNorms(  std::vector<CombiDataType>& dleft,
            std::vector<CombiDataType>& dright,
            std::vector<real>& Norms );

real l2Norm( std::vector<CombiDataType>& data );

int main( int argc, char** argv ){
  std::cout << argc << "\n";
  assert( argc == 6 );

  // mode ff, fc, cf, cc -> <format of first file><format of second file>
  // f -> combigrid plot file
  // c -> gene checkpoint file
  char* mode = argv[1];

  // files
  std::string filenameLeft( argv[2] );
  std::string filenameRight( argv[3] );
  std::string filenameError( argv[4] );
  std::string prefix( argv[5] );

  std::vector<CombiDataType> data1;
  std::vector<CombiDataType> data2;

  IndexVector res1;
  IndexVector res2;

  // get data and res of first file
  if( mode[0] == 'g' ){
    readCheckpoint( filenameLeft.c_str(), data1, res1 );
  } else if( mode[0] == 'f' ){
    readPlotFile( filenameLeft.c_str(), data1, res1 );
  } else{
    assert( !"wrong parameter" );
  }

  // get data and res of second file
  if( mode[1] == 'g' ){
    readCheckpoint( filenameRight.c_str(), data2, res2 );
  } else if( mode[1] == 'f' ){
    readPlotFile( filenameRight.c_str(), data2, res2 );
  } else{
    assert( !"wrong parameter" );
  }

  // check sizes
  assert( data1.size() == data2.size() );
  assert( res1 == res2 );

  // normalize with l2 norm
  real l2norm1 = l2Norm( data1 );
  real l2norm2 = l2Norm( data2 );
  std::cout << "l2 norm grid 1: " << l2norm1 << " l2 norm grid 2: " << l2norm2 << "\n";
  real tmp1 = 1.0/l2norm1;
  for( auto i=0; i<data1.size(); ++i )
    data1[i] *= tmp1;

  real tmp2 = 1.0/l2norm2;
  for( auto i=0; i<data2.size(); ++i )
    data2[i] *= tmp2;
  /*for( auto i=0; i<data1.size(); ++i ){
    std::cout << data1[i] << " ";
  }
  //std::cout << "\n data2";
  for( auto i=0; i<data2.size(); ++i ){
    std::cout << data2[i] << " ";
  }
  std::cout << "\n";
  */
  // calc l2 norm of absolute values
  real err = 0.0;
  for( auto i=0; i<data1.size(); ++i ){
    real tmp = std::abs( data1[i] ) - std::abs( data2[i] );
    
    if(std::abs(tmp)/std::abs(data1[i]) > 1e-12) std::cout <<" i: " << i << " value: " << std::abs(tmp)/std::abs(data1[i]) << " ";
    err += tmp*tmp;
  }
  std::cout << "\n";
  err = std::sqrt(err);
  // open file in append mode
  std::ofstream ofs( filenameError.c_str(), std::ofstream::app );

  // write prefix and err
  ofs << prefix << " " << err << std::endl;

  std::cout << "error " << err << std::endl;

  return 0;
}


real l2Norm( std::vector<CombiDataType>& data ){
  real nrm = 0.0;
  for( auto d : data )
    nrm += std::abs(d)*std::abs(d);

  return sqrt(nrm);
}


// read gene checkpoint
void readCheckpoint( const char* ggFileName,
                     std::vector<CombiDataType>& data,
                     IndexVector& resolution )
{
  // load gene grid from cp file
  std::cout << "reading GENE checkpoint " << ggFileName << std::endl;

  double tstart = timing();

  //check if file exists
  struct stat buffer;
  assert( stat( ggFileName, &buffer ) == 0 );

  std::ifstream ggFile( ggFileName );

  // read prec flag
  char prec[6];
  ggFile.read( prec, 6 );
  std::cout << "prec:(6)" << prec << std::endl;

  // read time and dt
  double mytime, dt;
  ggFile.read( (char*) &mytime, sizeof(double) );
  ggFile.read( (char*) &dt, sizeof(double) );
  std::cout << "mytime(" << sizeof(double) << "): " << mytime << std::endl;
  std::cout << "dt(" << sizeof(double) << "): " << dt << std::endl;

  // read resolution
  int res[6];
  ggFile.read( (char*) res, sizeof(int) * 6 );
  int ni0 = res[0];
  int nj0 = res[1];
  int nz0 = res[2];
  int nv0 = res[3];
  int nw0 = res[4];
  int n_spec = res[5];
  std::cout << "resolution:" << "\n" << "\t nx(" << sizeof(int) << ") " << ni0
      << "\n" << "\t ny(" << sizeof(int) << ") " << nj0 << "\n" << "\t nz("
      << sizeof(int) << ") " << nz0 << "\n" << "\t nv(" << sizeof(int) << ") "
      << nv0 << "\n" << "\t nw(" << sizeof(int) << ") " << nw0 << "\n"
      << "\t n_spec(" << sizeof(int) << ") " << n_spec << std::endl;
  assert( n_spec == 1 );
  //assert( nj0 == 1 );

  int dsize = ni0 * nj0 * nz0 * nv0 * nw0 * n_spec;
  std::cout << "size: " << dsize << std::endl;

  // read data
  data.resize(dsize);
  ggFile.read( (char*) &data[0], sizeof(CombiDataType) * dsize );

  ggFile.close();

  std::cout << "time to load cp file: " << timing() - tstart << "s"
      << std::endl;

  // set resolution
  resolution.resize(6);
  for( size_t i=0; i<6; ++i )
    resolution[i] = res[i];
  /*for(int i; i < data.size(); i++){
    //if(data[i] != CombiDataType(0))
      std::cout << data[i] << "\n";
  }*/

}


void readPlotFile( const char* pltFileName,
                   std::vector<CombiDataType>& data,
                   IndexVector& resolution ){
  // load gene grid from cp file
  std::cout << "reading plot file " << pltFileName << std::endl;

  double tstart = timing();

  //check if file exists
  struct stat buffer;
  assert( stat( pltFileName, &buffer ) == 0 );

  std::ifstream pltFile( pltFileName );

  // read dim and resolution
  int dim;
  pltFile.read( (char*) &dim, sizeof(int) );
  assert( dim == 6 );

  std::vector<int> res(dim);
  for( size_t i=0; i<dim; ++i )
    pltFile.read( (char*) &res[i], sizeof(int) );

  resolution.assign( res.begin(), res.end() );
  std::cout << "resolution " << resolution << std::endl;


  // calc data size
  size_t dsize = 1;
  for( auto r : resolution )
    dsize *= r;


  // read data
  std::vector<CombiDataType> tmp(dsize);
  pltFile.read( (char*) &tmp[0], sizeof(CombiDataType) * dsize );

  pltFile.close();

  std::cout << "time to load plot file: " << timing() - tstart << "s"
      << std::endl;

  // create multiarray view on tmp
  IndexVector shape( resolution.rbegin(), resolution.rend() );
  MultiArrayRef6 grid = createMultiArrayRef<CombiDataType,6>( &tmp[0], shape );

  // copy tmp to data without boundary points
  // copy data from local checkpoint to dfg
  // note that on the last process in some dimensions dfg is larger than the
  // local checkpoint
  int offsetY = 0;
  if(shape[4]>1){
    offsetY = 1;
  }
  for( size_t n=0; n < shape[0]; ++n ){ //n_spec
    for( size_t m=0; m < shape[1]-1; ++m ){ //w
      for( size_t l=0; l < shape[2]-1; ++l ){ //v
        for( size_t k=0; k < shape[3]-1; ++k ){ //z
          for( size_t j=0; j < shape[4]-offsetY; ++j ){ //y
            for( size_t i=0; i < shape[5]; ++i ){ //x
              data.push_back( grid[n][m][l][k][j][i] );
            }
          }
        }
      }
    }
  }

  // correct resolution
//  resolution[0] -= 1; //x
  if(shape[4]>1){
    resolution[1] -= 1; //y
  }
  resolution[2] -= 1; //z
  resolution[3] -= 1; //v
  resolution[4] -= 1; //w
  /*for(int i; i < data.size(); i++){
    //if(data[i] !=  CombiDataType(0))
      std::cout << data[i] << "\n";
  }*/
}


void
calcNorms( std::vector<CombiDataType>& dleft, std::vector<CombiDataType>& dright,
            std::vector<real>& Norms )
{
  Norms.resize( 9 );

  real el1( 0.0 ), el2( 0.0 ), emax( 0.0 );
  real ll1( 0.0 ), ll2( 0.0 ), lmax( 0.0 );
  real rl1( 0.0 ), rl2( 0.0 ), rmax( 0.0 );

  assert( dleft.size() == dright.size() );

  for( size_t i = 0; i < dleft.size(); ++i ){
    const CombiDataType ei = dleft[i] - dright[i];
    const real eiabs = std::abs( ei );
    const real liabs = std::abs( dleft[i] );
    const real riabs = std::abs( dright[i] );

    el1 += eiabs;
    el2 += eiabs * eiabs;
    if( eiabs > emax )
      emax = eiabs;

    ll1 += liabs;
    ll2 += liabs * liabs;
    if( liabs > lmax )
      lmax = liabs;

    rl1 += riabs;
    rl2 += riabs * riabs;
    if( riabs > rmax )
      rmax = riabs;
  }

  el2 = std::sqrt( el2 );
  ll2 = std::sqrt( ll2 );
  rl2 = std::sqrt( rl2 );

  Norms[0] = el1;
  Norms[1] = el2;
  Norms[2] = emax;
  Norms[3] = ll1;
  Norms[4] = ll2;
  Norms[5] = lmax;
  Norms[6] = rl1;
  Norms[7] = rl2;
  Norms[8] = rmax;
}

