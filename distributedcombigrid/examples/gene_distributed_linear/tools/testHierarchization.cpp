#include <mpi.h>

#include <sgpp/distributedcombigrid/fullgrid/FullGrid.hpp>
#include <sgpp/distributedcombigrid/hierarchization/Hierarchization.hpp>

using namespace combigrid;

void testSerial();
void testDistributed();
real computeNorm( FullGrid<complex>& fg1, FullGrid<complex>& fg2, int nrm );


int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  // testSerial();

  testDistributed();

  MPI_Finalize();

  return 0;
}





void testSerial(){
  // read Full Grid from file
  FullGrid<complex> fg( "test_fg.dat" );
  FullGrid<complex> fg_ref( "test_fg.dat" );

  // hierarchize
  Hierarchization::hierarchize( fg );

  // dehierarchize
  Hierarchization::dehierarchize( fg );

  // compare to reference
  std::cout << "L1 Norm = " << computeNorm( fg, fg_ref, 1 );
}


real computeNorm( FullGrid<complex>& fg1, FullGrid<complex>& fg2, int nrm ){
  std::vector<CombiDataType>& dleft = fg1.getElementVector();
  std::vector<CombiDataType>& dright = fg2.getElementVector();

  assert( dleft.size() == dright.size() );

  real l1(0), l2(0), linf(0);

  for( size_t i = 0; i < dleft.size(); ++i ){
    const CombiDataType ei = dleft[i] - dright[i];

    const real eiabs = std::abs( ei );

    l1 += eiabs;
    l2 += eiabs * eiabs;
    if( eiabs > linf )
      linf = eiabs;
  }

  l2 = std::sqrt( l2 );

  if( nrm == 1 )
    return l1;

  if( nrm == 2 )
    return l2;

  if( nrm == 0 )
    return linf;

  return l2;
}



