#define BOOST_TEST_DYN_LINK
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <boost/test/unit_test.hpp>
#include <cstdio>
#include <fstream>

#include "discotec/io/ParameterIO.hpp"
#include "discotec/mpi/MPISystem.hpp"
#include "test_helper.hpp"

using namespace combigrid;

BOOST_AUTO_TEST_SUITE(parameterio)

// Helper: rank 0 writes a temp file, all ranks barrier
std::string writeTempCtparam(const std::string& content) {
  std::string path = "test_ctparam_tmp.ini";
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    std::ofstream f(path);
    f << content;
    f.close();
  }
  MPI_Barrier(MPI_COMM_WORLD);
  return path;
}

void cleanup(const std::string& path) {
  MPI_Barrier(MPI_COMM_WORLD);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) std::remove(path.c_str());
}

BOOST_AUTO_TEST_CASE(test_readParameterFile_basic) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(1));

  auto path = writeTempCtparam(
      "[manager]\nngroup = 2\nnprocs = 4\n\n"
      "[ct]\ndim = 2\nlmin = 3 3\nlmax = 5 5\np = 2 2\nncombi = 10\n\n"
      "[application]\ndt = 0.01\n");

  auto [ngroup, nprocs, cfg] = readParameterFile(path, MPI_COMM_WORLD);

  BOOST_TEST(ngroup == 2u);
  BOOST_TEST(nprocs == 4u);
  BOOST_TEST(cfg.get<DimType>("ct.dim") == 2);
  BOOST_TEST(cfg.get<size_t>("ct.ncombi") == 10u);
  BOOST_TEST(cfg.get<double>("application.dt") == 0.01);

  cleanup(path);
}

BOOST_AUTO_TEST_CASE(test_readParameterFile_broadcast) {
  // All 9 ranks must see the same data after broadcast
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(9));

  auto path = writeTempCtparam(
      "[manager]\nngroup = 1\nnprocs = 9\n\n"
      "[ct]\ndim = 3\nlmin = 2 3 4\nlmax = 5 6 7\np = 1 3 3\nncombi = 42\n"
      "basis = biorthogonal_periodic\n\n"
      "[application]\ndt = 0.125\n");

  auto [ngroup, nprocs, cfg] = readParameterFile(path, MPI_COMM_WORLD);

  BOOST_TEST(ngroup == 1u);
  BOOST_TEST(nprocs == 9u);
  BOOST_TEST(cfg.get<DimType>("ct.dim") == 3);
  BOOST_TEST(cfg.get<size_t>("ct.ncombi") == 42u);
  BOOST_TEST(cfg.get<std::string>("ct.basis") == "biorthogonal_periodic");
  BOOST_TEST(cfg.get<double>("application.dt") == 0.125);

  DimType dim = 3;
  LevelVector lmin(dim), lmax(dim);
  cfg.get<std::string>("ct.lmin") >> lmin;
  cfg.get<std::string>("ct.lmax") >> lmax;
  BOOST_TEST(lmin[0] == 2);
  BOOST_TEST(lmin[1] == 3);
  BOOST_TEST(lmin[2] == 4);
  BOOST_TEST(lmax[0] == 5);
  BOOST_TEST(lmax[1] == 6);
  BOOST_TEST(lmax[2] == 7);

  cleanup(path);
}

BOOST_AUTO_TEST_CASE(test_readParameterFile_defaults) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(1));

  auto path = writeTempCtparam(
      "[manager]\nngroup = 1\nnprocs = 1\n\n"
      "[ct]\ndim = 3\nlmin = 2 2 2\nlmax = 4 4 4\nncombi = 5\n");

  auto [ngroup, nprocs, cfg] = readParameterFile(path, MPI_COMM_WORLD);

  BOOST_TEST(cfg.get<std::string>("ct.basis", "hat_periodic") == "hat_periodic");
  BOOST_TEST(cfg.get<uint32_t>("ct.chunkSize", 128) == 128);
  BOOST_TEST(!cfg.get_child_optional("ct.p"));
  BOOST_TEST(!cfg.get_child_optional("ct.boundary"));

  cleanup(path);
}

// buildCombiParametersFromConfig needs theMPISystem initialized with
// initWorldReusable(), which conflicts with the test harness's global
// MpiOnOff fixture.  These tests run as a standalone executable
// (e.g. the combi_workers_only example) rather than in the shared test
// harness.  Here we only test the validation logic which doesn't need
// the full MPI system.
BOOST_AUTO_TEST_CASE(test_parallelization_validation) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(1));

  // p = 2 2 → product 4, but nprocs = 1 → should throw
  auto path = writeTempCtparam(
      "[manager]\nngroup = 1\nnprocs = 1\n\n"
      "[ct]\ndim = 2\nlmin = 2 2\nlmax = 3 3\np = 2 2\nncombi = 5\n");

  auto [ngroup, nprocs, cfg] = readParameterFile(path, MPI_COMM_WORLD);

  // The validation happens inside buildCombiParametersFromConfig before
  // it touches theMPISystem, so this should throw even without init.
  BOOST_CHECK_THROW(buildCombiParametersFromConfig(cfg), std::invalid_argument);

  cleanup(path);
}

BOOST_AUTO_TEST_SUITE_END()
