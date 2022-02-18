#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include "sgpp/distributedcombigrid/third_level/NetworkUtils.hpp"
#include "test_helper.hpp"

static std::string host             = "localhost";
static std::vector<double> testData = std::vector<double>({0., 1., 2., 3., 4.});

using namespace combigrid;

void testBinarySendClient(MPI_Comm comm, unsigned short port) {
  // wait until server is set up
  MPI_Barrier(comm);

  // connect to server
  ClientSocket client(host, port);
  BOOST_TEST_CHECKPOINT("before client.init()");
  client.init();
  BOOST_TEST_CHECKPOINT("after client.init()");

  // send data
  client.sendallBinary(testData.data(), testData.size());
  BOOST_TEST_CHECKPOINT("client sent all data");
}

void testBinarySendRecvServer(MPI_Comm comm, unsigned short port) {
  // setup server
  ServerSocket server(port);
  BOOST_TEST_CHECKPOINT("before server.init()");
  server.init();
  BOOST_TEST_CHECKPOINT("after server.init()");

  // server is now ready to accept client
  MPI_Barrier(comm);
  std::shared_ptr<ClientSocket> client(server.acceptClient());
  BOOST_CHECK(client != nullptr);

  // receive data
  std::vector<double> recvData(testData.size());
  client->recvallBinaryAndCorrectInPlace(recvData.data(), recvData.size());
  BOOST_TEST_CHECKPOINT("server received all data");

  // check if received data has been successfully received
  std::cout << std::endl << "-- Receive Test:" << std::endl;
  for (size_t i = 0; i < recvData.size(); ++i) {
    std::cout << "received: " << recvData[i] << " expected: " << testData[i] << std::endl;
    BOOST_CHECK(recvData[i] == testData[i]);
  }
}

void testBinarySendReduceServer(MPI_Comm comm, unsigned short port) {
  // init server
  ServerSocket server(port);
  server.init();

  // server is now ready to accept client
  MPI_Barrier(comm);
  std::shared_ptr<ClientSocket> client(server.acceptClient());

  // reduce data
  std::vector<double> recvData(testData);
  client->recvallBinaryAndReduceInPlace<double>(recvData.data(),
      recvData.size(),
      [](const double& a, const double& b) -> double {return a + b;});

  // check if received data has been successfully reduced
  std::cout << std::endl << "-- Reduce Test:" << std::endl;
  for (size_t i = 0; i < recvData.size(); ++i) {
    std::cout << "received: " << recvData[i] << " expected: " << 2*testData[i] << std::endl;
    BOOST_CHECK(recvData[i] == 2*testData[i]);
  }
}


BOOST_FIXTURE_TEST_SUITE(networkutils, TestHelper::BarrierAtEnd, *boost::unit_test::timeout(60))

BOOST_AUTO_TEST_CASE(testBinarySendRecv) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(2));
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // create communicator with first two procs only
  MPI_Comm newComm = TestHelper::getComm(2);
  if (newComm == MPI_COMM_NULL)
    return;

  for (unsigned short port : {11111}) {
    BOOST_TEST_CHECKPOINT(std::to_string(port));
    MPI_Comm_rank(newComm, &rank);
    if (rank == 0)
      testBinarySendRecvServer(newComm, port);
    if (rank == 1)
      testBinarySendClient(newComm, port);
  }
}

BOOST_AUTO_TEST_CASE(testBinarySendReduce) {
  BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(2));
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // create communicator with first two procs only
  MPI_Comm newComm = TestHelper::getComm(2);
  if (newComm == MPI_COMM_NULL)
    return;

  for (unsigned short port : {11111}) {
    MPI_Comm_rank(newComm, &rank);
    if (rank == 0)
      testBinarySendReduceServer(newComm, port);
    if (rank == 1)
      testBinarySendClient(newComm, port);
  }
}

BOOST_AUTO_TEST_SUITE_END()
