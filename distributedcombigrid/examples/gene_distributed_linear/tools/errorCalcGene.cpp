/*
 * * errorCalc.cpp
 *
 *  Created on: Aug 19, 2014
 *      Author: heenemo
 */
#include <assert.h>
#include <fstream>
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include <iostream>
#include <vector>
#include "../src/GeneGrid.hpp"
#include "../src/CombiGeneConverter.hpp"
#include "timing.h"
#include <sys/stat.h>
#include "boost/lexical_cast.hpp"

using namespace combigrid;

FullGrid<complex>*
readGG( const char* ggFileName );
void
calcErrors( FullGrid<complex>& left, FullGrid<complex>& right,
            std::vector<real>& Norms );

int
main( int argc, char** argv )
{
  assert( argc == 8 );

  char* mode = argv[1];

  // file prefixes
  const char* leftPrefix = argv[2];
  const char* rightPrefix = argv[3];

  // output file
  const char* outfile = argv[4];

  // start, end, step
  int start = atoi( argv[5] );
  int end = atoi( argv[6] );
  int step = atoi( argv[7] );

  int count = (end - start) / step + 1;
  std::vector<int> stepslist(count);
  std::vector< std::vector<real> > normslist(count);

  #pragma omp parallel for
  for( int i = start; i <= end; i += step ){
    // report progress
    std::cout << "calculating step " << i/step << " of " << count << std::endl;

    std::string filenameLeft = leftPrefix
        + boost::lexical_cast<std::string>( i ) + ".dat";
    std::string filenameRight = rightPrefix
        + boost::lexical_cast<std::string>( i ) + ".dat";

    FullGrid<complex>* fg_left;
    FullGrid<complex>* fg_right;

    if( mode[0] == 'g' )
      fg_left = readGG( filenameLeft.c_str() );
    else if( mode[0] == 'f' ){
      fg_left = new FullGrid<complex>( filenameLeft.c_str() );
    }
    else{
      assert( !"wrong parameter" );
    }

    if( mode[1] == 'g' ){
      fg_right = readGG( filenameRight.c_str() );
    }
    else if( mode[1] == 'f' ){
      fg_right = new FullGrid<complex>( filenameRight.c_str() );
    }
    else{
      assert( !"wrong parameter" );
    }

    // sanity checks
    assert( fg_left->getLevels() == fg_right->getLevels() );
    assert( fg_left->returnBoundaryFlags() == fg_right->returnBoundaryFlags() );

    std::vector<real> norms(9);
    calcErrors( *fg_left, *fg_right, norms );

    // err norms are relative to left
    norms[0] = norms[0] / norms[3];
    norms[1] = norms[1] / norms[4];
    norms[2] = norms[2] / norms[5];

    stepslist[ i/step - 1 ] = i;
    normslist[ i/step - 1 ] = norms;

    delete fg_left;
    delete fg_right;
   }

  // write to file
  std::ofstream ofs( outfile );
  for( size_t i = 0; i < stepslist.size(); ++i ){
    ofs << stepslist[i];

    for( size_t j = 0; j < normslist[i].size(); ++j )
      ofs << " " << normslist[i][j];

    ofs << std::endl;
  }

  return 0;
}

FullGrid<complex>*
readGG( const char* ggFileName )
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
  assert( nj0 == 1 );

  int dsize = ni0 * nj0 * nz0 * nv0 * nw0 * n_spec;
  std::cout << "size: " << dsize << std::endl;

  // create gene grid
  GeneGrid gg( boost::extents[n_spec][nw0][nv0][nz0][nj0][ni0] );

  // read data
  std::vector<GeneComplex> ggData( dsize );
  ggFile.read( (char*) gg.origin(), sizeof(complex) * dsize );

  ggFile.close();

  // compute levelvector
  std::vector<int> levels( 5 );
  levels[0] = static_cast<int>( std::log( static_cast<double>( ni0 ) - 1.0 )
      / std::log( 2.0 ) ); // no boundary
  levels[1] = 1; //no boundary
  levels[2] = static_cast<int>( std::log( static_cast<double>( nz0 ) )
      / std::log( 2.0 ) ); // boundary
  levels[3] = static_cast<int>( std::log( static_cast<double>( nv0 ) )
      / std::log( 2.0 ) ); // no boundary
  levels[4] = static_cast<int>( std::log( static_cast<double>( nw0 ) )
      / std::log( 2.0 ) );

  std::cout << "levels: " << LevelVector( levels.begin(), levels.end() )
      << std::endl;

  std::vector<BoundaryType> boundary( 5 );
  boundary[0] = true;
  boundary[1] = false;
  boundary[2] = true;
  boundary[3] = false;
  boundary[4] = true;

  std::cout << "time to load cp file: " << timing() - tstart << "s"
      << std::endl;
  tstart = timing();

  // convert gg to fg
  FullGrid<complex>* fg = new FullGrid<complex>( 5, levels, boundary );
  CombiGeneConverter::GeneGridToFullGrid( gg, *fg );

  std::cout << "time to convert cp to fg: " << timing() - tstart << "s"
      << std::endl;

  return fg;
}

void
calcErrors( FullGrid<complex>& left, FullGrid<complex>& right,
            std::vector<real>& Norms )
{
  Norms.resize( 9 );

  real el1( 0.0 ), el2( 0.0 ), emax( 0.0 );
  real ll1( 0.0 ), ll2( 0.0 ), lmax( 0.0 );
  real rl1( 0.0 ), rl2( 0.0 ), rmax( 0.0 );

  std::vector<complex>& dleft = left.getElementVector();
  std::vector<complex>& dright = right.getElementVector();

  assert( dleft.size() == dright.size() );

  for( size_t i = 0; i < dleft.size(); ++i ){
    const complex ei = dleft[i] - dright[i];
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
