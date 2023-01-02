#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SGppDistributedCombigridModule
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include <boost/test/unit_test.hpp>

struct MpiOnOff {
  MpiOnOff() { MPI_Init(NULL, NULL); }
  ~MpiOnOff() { MPI_Finalize(); }
};

BOOST_GLOBAL_FIXTURE(MpiOnOff);
