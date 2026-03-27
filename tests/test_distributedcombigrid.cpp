#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SGppDistributedCombigridModule
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <boost/test/unit_test.hpp>

#include "discotec/mpi/MPISystem.hpp"

#ifdef DISCOTEC_USE_PALIWA
#include <Kokkos_Core.hpp>
#include <ddc/ddc.hpp>

struct KokkosDDCScopeGuard {
  KokkosDDCScopeGuard() {
    kokkosGuard_ = std::make_unique<Kokkos::ScopeGuard>();
    ddcGuard_ = std::make_unique<ddc::ScopeGuard>();
  }
  ~KokkosDDCScopeGuard() {
    ddcGuard_.reset();
    kokkosGuard_.reset();
  }
  std::unique_ptr<Kokkos::ScopeGuard> kokkosGuard_;
  std::unique_ptr<ddc::ScopeGuard> ddcGuard_;
};
#endif  // DISCOTEC_USE_PALIWA

using namespace combigrid;

BOOST_GLOBAL_FIXTURE(MpiOnOff);
#ifdef DISCOTEC_USE_PALIWA
BOOST_GLOBAL_FIXTURE(KokkosDDCScopeGuard);
#endif
