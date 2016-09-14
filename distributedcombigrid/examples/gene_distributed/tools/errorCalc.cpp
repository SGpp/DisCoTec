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
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include <iostream>
#include <vector>
#include <sys/stat.h>
#include "boost/lexical_cast.hpp"

using namespace combigrid;

void
calcErrors( FullGrid<CombiDataType>& left, FullGrid<CombiDataType>& right,
            std::vector<real>& Norms );

int
main( int argc, char** argv )
{
  assert( argc == 5 );

  // file prefixes
  std::string filenameLeft( argv[1] );
  std::string filenameRight( argv[2] );
  std::string filenameError( argv[3] );
  std::string prefix( argv[4] );

  FullGrid<CombiDataType> fg_left( filenameLeft.c_str() );
  FullGrid<CombiDataType> fg_right( filenameRight.c_str() );

  // sanity checks
  assert( fg_left.getLevels() == fg_right.getLevels() );
  assert( fg_left.returnBoundaryFlags() == fg_right.returnBoundaryFlags() );

  std::vector<real> norms(9);
  calcErrors( fg_left, fg_right, norms );

  // err norms are relative to left
  norms[0] = norms[0] / norms[3];
  norms[1] = norms[1] / norms[4];
  norms[2] = norms[2] / norms[5];

  // open file in append mode
  std::ofstream ofs( filenameError.c_str(), std::ofstream::app );

  // write prefix
  ofs << prefix << " ";

  // write norms
  for( size_t j = 0; j < norms.size(); ++j )
      ofs << " " << norms[j];

  ofs << std::endl;

  return 0;
}



void
calcErrors( FullGrid<CombiDataType>& left, FullGrid<CombiDataType>& right,
            std::vector<real>& Norms )
{
  Norms.resize( 9 );

  real el1( 0.0 ), el2( 0.0 ), emax( 0.0 );
  real ll1( 0.0 ), ll2( 0.0 ), lmax( 0.0 );
  real rl1( 0.0 ), rl2( 0.0 ), rmax( 0.0 );

  std::vector<CombiDataType>& dleft = left.getElementVector();
  std::vector<CombiDataType>& dright = right.getElementVector();

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
