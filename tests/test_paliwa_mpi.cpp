#define BOOST_TEST_DYN_LINK
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <boost/test/unit_test.hpp>
#include <cmath>
#include <complex>
#include <type_traits>
#include <cstdint>
#include <vector>

#include "discotec/fullgrid/DistributedFullGrid.hpp"
#include "discotec/hierarchization/DistributedHierarchization.hpp"
#include "test_helper.hpp"

#ifdef DISCOTEC_USE_PALIWA
using namespace combigrid;

namespace {
template <DimType DIM, typename Value = double>
void comparePaliwa(const LevelVector& levels, const LevelVector& lmin,
                   const std::vector<int>& procs, bool reverseRanks = false) {
  MPI_Comm comm = TestHelper::getComm(procs);
  if (comm == MPI_COMM_NULL) return;
  MPI_Comm reversed = MPI_COMM_NULL;
  if (reverseRanks) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    MPI_Comm ordered;
    MPI_Comm_split(comm, 0, size - 1 - rank, &ordered);
    std::vector<int> periods(DIM, 1);
    MPI_Cart_create(ordered, DIM, procs.data(), periods.data(), 0, &reversed);
    MPI_Comm_free(&ordered);
    comm = reversed;
  }
  {
    size_t count = 1;
    for (auto l : levels) count *= size_t{1} << l;
    int ranks;
    MPI_Comm_size(comm, &ranks);
    std::vector<Value> storage(count / static_cast<size_t>(ranks));
    DistributedFullGrid<Value, DIM> grid(DIM, levels, comm, std::vector<BoundaryType>(DIM, 1),
                                           storage.data(), procs, false);
    std::vector<Value> original(storage.size()), reference(storage.size());
    for (IndexType i = 0; i < grid.getNrLocalElements(); ++i) {
      // Irregular, non-separable data exposes axis swaps that a roundtrip misses.
      uint64_t x = static_cast<uint64_t>(grid.getGlobalLinearIndex(i)) + 1234;
      x += 0x9e3779b97f4a7c15ULL;
      x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
      x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
      x ^= x >> 31;
      original[static_cast<size_t>(i)] = 1.0 + static_cast<double>(x >> 11) * 0x1.0p-53;
      if constexpr (std::is_same_v<Value, std::complex<double>>) {
        original[static_cast<size_t>(i)].imag(static_cast<double>(x & 0xffff) / 65536.0);
      }
    }
    const std::vector<bool> dims(DIM, true);
    for (auto basis : {BasisFunctionType::HAT_PERIODIC, BasisFunctionType::BIORTHOGONAL_PERIODIC,
                       BasisFunctionType::FULLWEIGHTING_PERIODIC}) {
      const std::vector<BasisFunctionType> bases(DIM, basis);
      std::copy(original.begin(), original.end(), storage.begin());
      DistributedHierarchization::hierarchize(grid, dims, bases, lmin,
                                               HierarchizationBackend::DISCOTEC);
      std::copy(storage.begin(), storage.end(), reference.begin());
      std::copy(original.begin(), original.end(), storage.begin());
      DistributedHierarchization::hierarchize(grid, dims, bases, lmin,
                                               HierarchizationBackend::PALIWA);
      double coefficientError = 0, inverseError = 0;
      for (size_t i = 0; i < storage.size(); ++i) {
        BOOST_REQUIRE(std::isfinite(std::abs(storage[i])));
        coefficientError = std::max(coefficientError, std::abs(storage[i] - reference[i]));
      }
      // Start inverse from independently generated native coefficients.
      std::copy(reference.begin(), reference.end(), storage.begin());
      DistributedHierarchization::dehierarchize(grid, dims, bases, lmin,
                                                 HierarchizationBackend::PALIWA);
      for (size_t i = 0; i < storage.size(); ++i) {
        BOOST_REQUIRE(std::isfinite(std::abs(storage[i])));
        inverseError = std::max(inverseError, std::abs(storage[i] - original[i]));
      }
      MPI_Allreduce(MPI_IN_PLACE, &coefficientError, 1, MPI_DOUBLE, MPI_MAX, comm);
      MPI_Allreduce(MPI_IN_PLACE, &inverseError, 1, MPI_DOUBLE, MPI_MAX, comm);
      BOOST_TEST(coefficientError < 1e-11);
      BOOST_TEST(inverseError < 1e-11);
    }
  }
  if (reversed != MPI_COMM_NULL) MPI_Comm_free(&reversed);
}
}  // namespace

BOOST_FIXTURE_TEST_SUITE(paliwa_mpi, TestHelper::BarrierAtEnd)

BOOST_AUTO_TEST_CASE(anisotropic_coefficients_serial) {
  comparePaliwa<2>({5, 3}, {0, 0}, {1, 1});
  comparePaliwa<3>({4, 3, 2}, {2, 1, 0}, {1, 1, 1});
}

#ifdef PALIWA_WITH_MPI
BOOST_AUTO_TEST_CASE(slabs_two_ranks) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(2));
  comparePaliwa<2>({5, 3}, {0, 0}, {2, 1});
  comparePaliwa<2>({5, 3}, {5, 1}, {2, 1});
  comparePaliwa<2>({5, 3}, {2, 1}, {1, 2});
  comparePaliwa<3>({4, 3, 2}, {1, 0, 1}, {1, 2, 1});
}

BOOST_AUTO_TEST_CASE(slabs_four_ranks) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(4));
  comparePaliwa<2>({5, 3}, {0, 0}, {4, 1});
  comparePaliwa<2>({5, 3}, {2, 1}, {1, 4});
  comparePaliwa<3>({4, 3, 2}, {1, 0, 1}, {1, 4, 1});
}

BOOST_AUTO_TEST_CASE(complex_four_ranks) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(4));
  comparePaliwa<2, std::complex<double>>({5, 3}, {1, 0}, {2, 2});
}

BOOST_AUTO_TEST_CASE(one_point_per_rank) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(4));
  comparePaliwa<2>({2, 3}, {0, 0}, {4, 1});
}

BOOST_AUTO_TEST_CASE(cartesian_four_ranks_reordered) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(4));
  comparePaliwa<2>({5, 3}, {1, 0}, {2, 2}, true);
  comparePaliwa<3>({4, 3, 2}, {0, 1, 0}, {2, 1, 2}, true);
}

BOOST_AUTO_TEST_CASE(reject_nonuniform_partition) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(2));
  const std::vector<int> procs{2, 1};
  const auto comm = TestHelper::getComm(procs);
  if (comm == MPI_COMM_NULL) return;
  // All ranks see the same decomposition and must reject before communication.
  std::vector<double> storage(16 * 8);
  DistributedFullGrid<double, 2> grid(2, {4, 3}, comm, {1, 1}, storage.data(), procs,
                                      false, {{0, 3}, {0}});
  const std::vector<BasisFunctionType> bases(2, BasisFunctionType::HAT_PERIODIC);
  BOOST_CHECK_THROW(DistributedHierarchization::hierarchize(
                        grid, {true, true}, bases, {0, 0}, HierarchizationBackend::PALIWA),
                    std::runtime_error);
  BOOST_CHECK_THROW(DistributedHierarchization::dehierarchize(
                        grid, {true, true}, bases, {0, 0}, HierarchizationBackend::PALIWA),
                    std::runtime_error);
}
#endif
BOOST_AUTO_TEST_SUITE_END()
#endif
