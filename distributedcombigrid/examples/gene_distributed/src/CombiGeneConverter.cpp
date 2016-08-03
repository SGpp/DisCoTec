/*
 * CombiGeneConverter.cpp
 *
 *  Created on: Jan 14, 2014
 *      Author: heenemo
 */

#include "CombiGeneConverter.hpp"
#include <boost/multi_array.hpp>
#include <cmath>
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include <cassert>

namespace combigrid {

CombiGeneConverter::CombiGeneConverter() {
  // TODO Auto-generated constructor stub

}

CombiGeneConverter::~CombiGeneConverter() {
  // TODO Auto-generated destructor stub
}

void CombiGeneConverter::GeneGridToFullGrid(  GeneGrid& gg,
                                          FullGrid<complex> &fg ){
  // check if gg and fg compatible
  assert( fg.getDimension() == 5 );

  // check boundary
  const std::vector<bool>& boundary = fg.returnBoundaryFlags();
  assert( boundary[0] == true );
  assert( boundary[1] == false );
  assert( boundary[2] == true );
  assert( boundary[3] == false );
  assert( boundary[4] == true );

  // check sizes
  size_t nspec = gg.shape()[0];
  assert( nspec == 1 );

  size_t nw    = gg.shape()[1];
  size_t nv    = gg.shape()[2];
  size_t nz    = gg.shape()[3];
  size_t ny    = gg.shape()[4];
  size_t nx    = gg.shape()[5];

  const LevelVector& levels = fg.getLevels();

  assert( levels[0] == static_cast<int>( std::log( static_cast<double>(nx) - 1.0 )
                        / std::log( 2.0 ) ) ); // boundary
  assert( levels[1] == 1 );
  assert( ny == 1 );
  assert( levels[2] == static_cast<int>( std::log( static_cast<double>(nz) )
                          / std::log( 2.0 ) ) ); // boundary
  assert( levels[3] == static_cast<int>( std::log( static_cast<double>(nv) )
                          / std::log( 2.0 ) ) ); // no boundary
  assert( levels[4] == static_cast<int>( std::log( static_cast<double>(nw) )
                            / std::log( 2.0 ) ) );  // boundary





  // create new array of nz0 extent +1
  const size_t* oshape = gg.shape();
  // shape: n_spec nw0 nv0 nz0 nky0 nx0
  GeneGrid tmp0(boost::extents[ oshape[0]   ][ oshape[1] ][ oshape[2] ]
                                     [ oshape[3] ][ oshape[4] ][ oshape[5] ]);

  // rearrange x data from 0 1 3 2 4 -3 -2 -1 to -3,-2,-1,0,1,2,3
  // nx must be an odd number
  assert( oshape[5]%2 == 1 );
  size_t truncInd = oshape[5] / 2 + 1;
  for( size_t n=0; n<oshape[0]; ++n )
    for( size_t m=0; m<oshape[1]; ++m )
      for( size_t l=0; l<oshape[2]; ++l )
        for( size_t k=0; k<oshape[3]; ++k )
          for( size_t j=0; j<oshape[4]; ++j )
            for( size_t i=0; i < truncInd - 1; ++i )
              tmp0[n][m][l][k][j][i] = gg[n][m][l][k][j][ truncInd + i ];

  for( size_t n=0; n<oshape[0]; ++n )
    for( size_t m=0; m<oshape[1]; ++m )
      for( size_t l=0; l<oshape[2]; ++l )
        for( size_t k=0; k<oshape[3]; ++k )
          for( size_t j=0; j<oshape[4]; ++j )
            for( size_t i=truncInd - 1; i < oshape[5]; ++i )
              tmp0[n][m][l][k][j][i] = gg[n][m][l][k][j][ i - truncInd + 1 ];

  // print to file
  //write to file
  //FILE *fp;
  //fp=fopen("./cp_after_x", "wb");
  size_t dsize = oshape[0]*oshape[1]*oshape[2]*oshape[3]*oshape[4]*oshape[5];
  //fwrite(tmp0.data(),sizeof(GeneComplex),dsize,fp);

  GeneGrid tmp(boost::extents[ oshape[0]     ][ oshape[1] ][ oshape[2] ]
                                     [ oshape[3] + 1 ][ oshape[4] ][ oshape[5] ]);

  // adapt z
  const size_t* nshape = tmp.shape();

  // copy data to new array
  for( size_t n=0; n<nshape[0]; ++n )
    for( size_t m=0; m<nshape[1]; ++m )
      for( size_t l=0; l<nshape[2]; ++l )
        for( size_t k=0; k<oshape[3]; ++k )
          for( size_t j=0; j<nshape[4]; ++j )
            for( size_t i=0; i<nshape[5]; ++i )
              tmp[n][m][l][k][j][i] = tmp0[n][m][l][k][j][i];

  // fill new z layer with data
  for( size_t n=0; n<nshape[0]; ++n )
    for( size_t m=0; m<nshape[1]; ++m )
      for( size_t l=0; l<nshape[2]; ++l )
          for( size_t j=0; j<nshape[4]; ++j )
            for( size_t i=0; i<nshape[5]; ++i )
              tmp[n][m][l][ nshape[3]-1 ][j][i] = tmp0[n][m][l][ oshape[3]-1 ][j][i];

  //write to file
  //FILE *fp2;
  //fp2=fopen("./cp_after_z1", "wb");
  dsize = nshape[0]*nshape[1]*nshape[2]*nshape[3]*nshape[4]*nshape[5];
  //fwrite(tmp.data(),sizeof(GeneComplex),dsize,fp2);

  real shat = 0.7960; // todo get shat from parfile
  real kymin = 0.3000; // todo get kymin from parfile
  real lx = 4.18760; // todo get lx from parfile
  int N = int( round( shat * kymin * lx ) );
  int ky0_ind = 1; // todo get ky0_ind from parfile

  int kx_offset = static_cast<int>(nshape[5]) / 2;

  // todo indexing wrong when nky > 1
  for( int j = ky0_ind; j < static_cast<int>(nshape[4]) + ky0_ind; ++j ){
    for( int i = -kx_offset; i < static_cast<int>(nshape[5]) -1 - kx_offset; ++i ){
      int kx_ = (i + N * j);
      real factor = std::pow(-1.0,N*j);
      //std::cout << "tmp" << nshape[3]-1<<","<<j-1<<","<<i+kx_offset<< " "
      //          << "gg" << 0<<","<<j-1<<","<<kx_+kx_offset<<" "<<factor<<std::endl;

      for( size_t n=0; n<nshape[0]; ++n )
        for( size_t m=0; m<nshape[1]; ++m )
          for( size_t l=0; l<nshape[2]; ++l ){

            tmp[n][m][l][ nshape[3]-1 ][ j-1 ][ i+kx_offset ].r
              = tmp0[n][m][l][0][ j-1 ][ kx_ + kx_offset ].r * factor;
            tmp[n][m][l][ nshape[3]-1 ][ j-1 ][ i+kx_offset ].i
              = tmp0[n][m][l][0][ j-1 ][ kx_ + kx_offset ].i * factor;
          }
    }
  }

  //write to file
  //FILE *fp3;
  //fp3=fopen("./cp_after_z2", "wb");
  dsize = nshape[0]*nshape[1]*nshape[2]*nshape[3]*nshape[4]*nshape[5];
  //fwrite(tmp.data(),sizeof(GeneComplex),dsize,fp3);

  // throw away first entry of v
  GeneGrid tmp2(boost::extents[ nshape[0] ][ nshape[1] ][ nshape[2]-1 ]
                                      [ nshape[3] ][ nshape[4] ][ nshape[5]   ]);
  const size_t* n2shape = tmp2.shape();

  // copy data to new array
  for( size_t n=0; n<n2shape[0]; ++n )
    for( size_t m=0; m<n2shape[1]; ++m )
      for( size_t l=0; l<n2shape[2]; ++l )
        for( size_t k=0; k<n2shape[3]; ++k )
          for( size_t j=0; j<n2shape[4]; ++j )
            for( size_t i=0; i<n2shape[5]; ++i )
              tmp2[n][m][l][k][j][i] = tmp[n][m][l+1][k][j][i];

  //write to file
  //FILE *fp4;
  //fp4=fopen("./cp_after_v", "wb");
  dsize = n2shape[0]*n2shape[1]*n2shape[2]*n2shape[3]*n2shape[4]*n2shape[5];
  //fwrite(tmp2.data(),sizeof(GeneComplex),dsize,fp4);

  // enlarge w
  GeneGrid tmp3(boost::extents[ n2shape[0] ][ n2shape[1]+1 ][ n2shape[2] ]
                                      [ n2shape[3] ][ n2shape[4] ][ n2shape[5] ]);
  const size_t* n3shape = tmp3.shape();

  // copy data to new array
  for( size_t n=0; n<n2shape[0]; ++n )
    for( size_t m=0; m<n2shape[1]; ++m )
      for( size_t l=0; l<n2shape[2]; ++l )
        for( size_t k=0; k<n2shape[3]; ++k )
          for( size_t j=0; j<n2shape[4]; ++j )
            for( size_t i=0; i<n2shape[5]; ++i )
              tmp3[n][m][l][k][j][i] = tmp2[n][m][l][k][j][i];

  // set to new w values to zero
  for( size_t n=0; n<n2shape[0]; ++n )
    for( size_t l=0; l<n2shape[2]; ++l )
      for( size_t k=0; k<n2shape[3]; ++k )
        for( size_t j=0; j<n2shape[4]; ++j )
          for( size_t i=0; i<n2shape[5]; ++i ){
            tmp3[n][ n3shape[1]-1 ][l][k][j][i].r = 0.0;
            tmp3[n][ n3shape[1]-1 ][l][k][j][i].i = 0.0;
          }

  //write to file
  //FILE *fp5;
  //fp5=fopen("./cp_after_w", "wb");
  dsize = tmp3.num_elements();
  //fwrite(tmp3.data(),sizeof(GeneComplex),dsize,fp5);

  // create fullgrid
  fg.createFullGrid();

  // check if size is compatible
  //std::cout << "fg size = " << fg.getNrElements()
  //          << " tmp3 size = " << dsize << std::endl;
  assert( static_cast<size_t>( fg.getNrElements() ) == dsize );

  // copy values from multiarray to fullgrid
  std::vector<complex>& fgData = fg.getElementVector();
  GeneComplex* ggData = tmp3.data();
  for( size_t i = 0; i < static_cast<size_t>( fg.getNrElements() ); ++i ){
    fgData[i].real( ggData[i].r );
    fgData[i].imag( ggData[i].i );
  }

  //std::cout << "converted gene grid to ct fullgrid" << std::endl;
}


void CombiGeneConverter::FullGridToGeneGrid( FullGrid<complex> &fg,
                                                    GeneGrid& gg ){
  size_t nspec = gg.shape()[0];
  size_t nw    = gg.shape()[1];
  size_t nv    = gg.shape()[2];
  size_t nz    = gg.shape()[3];
  size_t ny    = gg.shape()[4];
  size_t nx    = gg.shape()[5];

  // check if sizes fit
  assert( static_cast<size_t>( std::pow(2.0,fg.getLevels()[0]) ) == nx - 1 );
  assert( 1 == ny );
  assert( static_cast<size_t>( std::pow(2.0,fg.getLevels()[2]) ) == nz );
  assert( static_cast<size_t>( std::pow(2.0,fg.getLevels()[3]) ) == nv );
  assert( static_cast<size_t>( std::pow(2.0,fg.getLevels()[4]) ) == nw );

  // nx must be an odd number
  assert( nx%2 == 1 );

  // create 6-dim view to fg data
  typedef boost::multi_array_ref<complex, 6> array_type_ref;
  complex* data = &fg.getElementVector()[0];
  array_type_ref fg_data( data, boost::extents[nspec][nw+1][nv-1][nz+1 ][ny  ][nx  ] );

  // copy data to gene grid
  // rearrange x data from -3,-2,-1,0,1,2,3 to 0 1 3 4 -3 -2 -1
  size_t offset = nx / 2;
  for( size_t n=0; n<nspec; ++n )
    for( size_t m=0; m<nw; ++m )
      for( size_t l=1; l<nv; ++l )
        for( size_t k=0; k<nz; ++k )
          for( size_t j=0; j<ny; ++j )
            for( size_t i=0; i < offset + 1; ++i ){
              gg[n][m][l][k][j][i].r = fg_data[n][m][l-1][k][j][ offset + i ].real();
              gg[n][m][l][k][j][i].i = fg_data[n][m][l-1][k][j][ offset + i ].imag();
            }

  for( size_t n=0; n<nspec; ++n )
    for( size_t m=0; m<nw; ++m )
      for( size_t l=1; l<nv; ++l )
        for( size_t k=0; k<nz; ++k )
          for( size_t j=0; j<ny; ++j )
            for( size_t i=offset+1; i < nx; ++i ){
              gg[n][m][l][k][j][i].r = fg_data[n][m][l-1][k][j][ i - offset - 1 ].real();
              gg[n][m][l][k][j][i].i = fg_data[n][m][l-1][k][j][ i - offset - 1 ].imag();
            }

  // fill v0 with zero
  for( size_t n=0; n<nspec; ++n )
    for( size_t m=0; m<nw; ++m )
        for( size_t k=0; k<nz; ++k )
          for( size_t j=0; j<ny; ++j )
            for( size_t i=0; i<nx; ++i ){
              gg[n][m][0][k][j][i].r = 0.0;
              gg[n][m][0][k][j][i].i = 0.0;
            }
}

} /* namespace combigrid */
