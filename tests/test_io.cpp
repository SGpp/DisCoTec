#define BOOST_TEST_DYN_LINK

#include <boost/serialization/export.hpp>
#include <boost/test/unit_test.hpp>

#include "io/MPICompression.hpp"
#include "io/MPIInputOutput.hpp"
#include "mpi/MPISystem.hpp"
#include "test_helper.hpp"
#include "utils/MonteCarlo.hpp"

using namespace combigrid;

void checkCompressionFrame(size_t numValues = 1000000) {
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

  for (bool randomValues : std::vector<bool>({false, true})) {
    numValues = numValues + static_cast<size_t>(montecarlo::getRandomNumber(-rank - 100, rank + 1));
    std::vector<double> originalValues;
    if (randomValues) {
      // generate random vector of doubles
      originalValues = std::move(montecarlo::getRandomCoordinates(1, numValues)[0]);
    } else {
      originalValues.resize(numValues, 1234.56);
    }

    std::vector<char> compressedString;
    mpiio::compressBufferToLZ4Frame(originalValues.data(), numValues, compressedString);

    if (randomValues) {
      // assume size increase for non-compressible, random data
      BOOST_CHECK_GT(compressedString.size(), numValues * sizeof(double));
    } else {
      // assume size decrease for compressible, same-value data
      BOOST_CHECK_LT(compressedString.size(), numValues * sizeof(double));
    }

    std::vector<double> decompressedValues(originalValues.size());

    mpiio::decompressLZ4FrameToBuffer(std::move(compressedString), decompressedValues);

    BOOST_CHECK_EQUAL(originalValues.size(), decompressedValues.size());
    BOOST_CHECK_EQUAL_COLLECTIONS(originalValues.begin(), originalValues.end(),
                                  decompressedValues.begin(), decompressedValues.end());
  }
  combigrid::Stats::finalize();
  MPI_Barrier(comm);
  BOOST_CHECK(!TestHelper::testStrayMessages(localComm));
  BOOST_CHECK(!TestHelper::testStrayMessages(comm));
}

void checkCompressionWithHeader(size_t numValues = 1000000) {
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

  for (bool randomValues : std::vector<bool>({false, true})) {
    numValues = numValues + static_cast<size_t>(montecarlo::getRandomNumber(-rank - 100, rank + 1));
    std::vector<double> originalValues;
    if (randomValues) {
      // generate random vector of doubles
      originalValues = std::move(montecarlo::getRandomCoordinates(1, numValues)[0]);
    } else {
      originalValues.resize(numValues, 1234.56);
    }

    std::string filename =
        "compression_testfile_" + std::to_string(theMPISystem()->getProcessGroupNumber());
    mpiio::writeCompressedValuesConsecutive(
        originalValues.data(), originalValues.size(), filename, localComm);

    std::vector<double> decompressedValues(originalValues.size());

    auto numValuesObtained = mpiio::readCompressedValuesConsecutive(
        decompressedValues.data(), decompressedValues.size(), filename, localComm);

    BOOST_CHECK_EQUAL(originalValues.size(), numValuesObtained);
    BOOST_CHECK_EQUAL_COLLECTIONS(originalValues.begin(), originalValues.end(),
                                  decompressedValues.begin(), decompressedValues.end());

    numValuesObtained = mpiio::readReduceCompressedValuesConsecutive(
        decompressedValues.data(), decompressedValues.size(), filename, localComm,
        std::plus<double>{});

    // double the original values for comparison
    for (auto& val : originalValues) {
      val *= 2.0;
    }
    BOOST_CHECK_EQUAL(originalValues.size(), numValuesObtained);
    BOOST_CHECK_EQUAL_COLLECTIONS(originalValues.begin(), originalValues.end(),
                                  decompressedValues.begin(), decompressedValues.end());
  }
  combigrid::Stats::finalize();
  MPI_Barrier(comm);
  BOOST_CHECK(!TestHelper::testStrayMessages(localComm));
  BOOST_CHECK(!TestHelper::testStrayMessages(comm));
}

BOOST_FIXTURE_TEST_SUITE(io, TestHelper::BarrierAtEnd, *boost::unit_test::timeout(60))

BOOST_AUTO_TEST_CASE(test_1, *boost::unit_test::timeout(60)) {
#ifdef DISCOTEC_USE_LZ4
  checkCompressionFrame();
#endif
}

BOOST_AUTO_TEST_CASE(test_2, *boost::unit_test::timeout(60)) {
#ifdef DISCOTEC_USE_LZ4

  BOOST_CHECK_NO_THROW(
      try { checkCompressionWithHeader(); } catch (std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        throw e;
      });

#endif
}

BOOST_AUTO_TEST_SUITE_END()
