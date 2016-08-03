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
#include "boost/multi_array.hpp"

using namespace combigrid;


void testBoundaryZ(FullGrid<CombiDataType>& fg, std::string filePrefix ){
  DimType dim = fg.getDimension();

  const IndexVector fgSizes = fg.getSizes();

  typedef boost::multi_array_ref<complex, 6> DFGGridRef;
  DFGGridRef fgData( fg.getData(),
                      boost::extents[ fgSizes[0] ][ fgSizes[1] ][ fgSizes[2] ]
                                    [ fgSizes[3] ][ fgSizes[4] ][ fgSizes[5] ] );

  // loop over all dimensions except z
  // x, y, z, v, w, spec
  bool check = true;

  for( size_t xi=0; xi<fgSizes[0]-1; ++xi )
    for( size_t yi=0; yi<fgSizes[1]; ++yi )
      for( size_t vi=0; vi<fgSizes[3]; ++vi )
        for( size_t wi=0; wi<fgSizes[4]; ++wi )
          for( size_t ni=0; ni<fgSizes[5]; ++ni ){
           if( fgData[xi+1][yi][0][vi][wi][ni]
               != fgData[xi][yi][fgSizes[2]-1][vi][wi][ni] ){
               std::cout << "left: "
                         << fgData[xi+1][yi][0][vi][wi][ni]
                         << " right: "
                         << fgData[xi][yi][fgSizes[2]-1][vi][wi][ni]
                         << std::endl;

               std::cout << "xi = " << xi
                         << " yi = " << yi
                         << " vi = " << vi
                         << " wi = " << wi
                         << " ni = " << ni
                         << std::endl;

               check = false;
               assert( check );
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
