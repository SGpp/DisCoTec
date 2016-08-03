#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <mpi.h>
#include <iostream>

static bool CHECK_NPROCS( int nprocs ){
    int size;
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    return size == nprocs;
}

BOOST_AUTO_TEST_CASE( bla4 ){
  if( CHECK_NPROCS(4) ){
    std::cout << "I'm with 4" << std::endl;

    BOOST_CHECK( true );
  }
}


BOOST_AUTO_TEST_CASE( bla8 ){
  if( CHECK_NPROCS(8) ){
    std::cout << "I'm with 8" << std::endl;

    BOOST_CHECK( true );
  }
}
