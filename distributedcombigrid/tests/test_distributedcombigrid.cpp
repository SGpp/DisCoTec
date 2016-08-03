#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SGppDistributedCombigridModule
#include <boost/test/unit_test.hpp>
#include <mpi.h>

struct MpiOnOff {
  MpiOnOff()   { MPI_Init( NULL, NULL); }
  ~MpiOnOff()  { MPI_Finalize(); }
};

BOOST_GLOBAL_FIXTURE( MpiOnOff );

