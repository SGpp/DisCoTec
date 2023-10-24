#define BOOST_TEST_DYN_LINK

#include <boost/serialization/export.hpp>
#include <boost/test/unit_test.hpp>

#include "io/MPICompression.hpp"
#include "mpi/MPISystem.hpp"
#include "test_helper.hpp"
#include "utils/MonteCarlo.hpp"

using namespace combigrid;

void checkCompression() {
  auto commSize = getCommSize(MPI_COMM_WORLD);
  int halfTheRanks = commSize / 2;
  CommunicatorType comm = TestHelper::getComm(halfTheRanks * 2);
  if (comm == MPI_COMM_NULL) {
    BOOST_TEST_CHECKPOINT("drop out of test comm");
    return;
  }

  combigrid::Stats::initialize();
  theMPISystem()->initWorldReusable(comm, 2, halfTheRanks, false);
  auto rank = theMPISystem()->getLocalRank();
  auto localComm = theMPISystem()->getLocalComm();

  // generate random vector of doubles
  size_t numValues =
      1000000 + static_cast<size_t>(montecarlo::getRandomNumber(-rank - 100, rank + 1));
  const std::vector<double> originalValues =
      std::move(montecarlo::getRandomCoordinates(1, numValues)[0]);

  std::vector<char> compressedString;
  mpiio::compressBufferToLZ4Frame(originalValues.data(), numValues, localComm, compressedString);

  // assume size increase for non-compressible, random data
  BOOST_CHECK_GT(compressedString.size(), numValues * sizeof(double));

  std::vector<double> decompressedValues(originalValues.size());

  mpiio::decompressLZ4FrameToBuffer(std::move(compressedString), localComm, decompressedValues);

  BOOST_CHECK_EQUAL(originalValues.size(), decompressedValues.size());
  BOOST_CHECK_EQUAL_COLLECTIONS(originalValues.begin(), originalValues.end(),
                                decompressedValues.begin(), decompressedValues.end());

  combigrid::Stats::finalize();
  MPI_Barrier(comm);
  BOOST_CHECK(!TestHelper::testStrayMessages(localComm));
  BOOST_CHECK(!TestHelper::testStrayMessages(comm));
}

BOOST_FIXTURE_TEST_SUITE(io, TestHelper::BarrierAtEnd, *boost::unit_test::timeout(60))

BOOST_AUTO_TEST_CASE(test_1, *boost::unit_test::timeout(60)) {
#ifdef DISCOTEC_USE_LZ4
  checkCompression();
#endif
  MPI_Barrier(MPI_COMM_WORLD);
  BOOST_CHECK(!TestHelper::testStrayMessages());
}

BOOST_AUTO_TEST_SUITE_END()
