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
#include "sgpp/distributedcombigrid/fullgrid/MultiArray.hpp"
#include <iostream>
#include <vector>
#include <sys/stat.h>
#include "boost/lexical_cast.hpp"
#include "boost/multi_array.hpp"

using namespace combigrid;


void testBoundaryZ(FullGrid<CombiDataType>& fg, std::string filePrefix ){
  const double tol = 1e-15;

  MultiArrayRef6 fgData = createMultiArrayRef<CombiDataType,6>( fg );

  // loop over all dimensions except z
  bool check = true;

  const size_t* fgShape = fgData.shape();
  for( size_t n=0; n < fgShape[0]; ++n ) //n_spec
    for( size_t m=0; m < fgShape[1]; ++m ) //w
      for( size_t l=0; l < fgShape[2]; ++l ) //v
          for( size_t j=0; j < fgShape[4]; ++j ) //y
            for( size_t i=0; i < fgShape[5]; ++i ){ //x
              if( std::abs( fgData[n][m][l][ fgShape[3]-1 ][j][i] -
                    fgData[n][m][l][0][j][ (i + 1)%fgShape[5] ] )
                    > tol ){
                std::cout << "left: " << fgData[n][m][l][ fgShape[3]-1 ][j][i]
                          << " right: " << fgData[n][m][l][0][j][ (i + 1)%fgShape[5] ]
                          << std::endl;

                  std::cout << "i = " << i
                            << " j = " << j
                            << " l = " << l
                            << " m = " << m
                            << " n = " << n
                            << std::endl;

                  check = false;
                  //assert( check );
              }
            }
}


int
main( int argc, char** argv )
{
  assert( argc == 3 );

  // file prefixes
  std::string filenameLeft( argv[1] );
  std::string prefix( argv[2] );

  FullGrid<CombiDataType> fg( filenameLeft.c_str() );

  testBoundaryZ( fg, prefix );

  return 0;
}
