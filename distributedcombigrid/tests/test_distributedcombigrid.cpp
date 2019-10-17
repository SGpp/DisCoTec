#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SGppDistributedCombigridModule
#include <mpi.h>
#include <boost/test/unit_test.hpp>

struct MpiOnOff {
  MpiOnOff() { MPI_Init(NULL, NULL); }
  ~MpiOnOff() { MPI_Finalize(); }
};

BOOST_GLOBAL_FIXTURE(MpiOnOff);
