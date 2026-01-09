#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SGppDistributedCombigridModule
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <boost/test/unit_test.hpp>
#ifdef DISCOTEC_USE_PALIWA
#include <Kokkos_Core.hpp>
#include <ddc/ddc.hpp>
#endif // DISCOTEC_USE_PALIWA

#include "mpi/MPISystem.hpp"

using namespace combigrid;

BOOST_GLOBAL_FIXTURE(MpiOnOff);
#ifdef DISCOTEC_USE_PALIWA
using kokkos_scope_guard = Kokkos::ScopeGuard;
using ddc_scope_guard = ddc::ScopeGuard;
BOOST_GLOBAL_FIXTURE(kokkos_scope_guard);
BOOST_GLOBAL_FIXTURE(ddc_scope_guard);
#endif // DISCOTEC_USE_PALIWA
