#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "sgpp/distributedcombigrid/utils/Stats.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "test_helper.hpp"


void checkStats(int size){
    BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(size));
    BOOST_CHECK( size % 2 == 1);

    combigrid::CommunicatorType comm = TestHelper::getComm(size);

    if(comm == MPI_COMM_NULL){
        return;
    }
    
    BOOST_TEST_WARN(false);
    combigrid::Stats::initialize();
    BOOST_TEST_WARN(true);
    combigrid::Stats::finalize();
    MPI_Barrier(comm);
    combigrid::Stats::write("./Stats_output", comm);
    MPI_Barrier(comm);
}

BOOST_AUTO_TEST_SUITE(stats)

BOOST_AUTO_TEST_CASE(test_3, * boost::unit_test::timeout(20)) {
    checkStats(3);
}

BOOST_AUTO_TEST_CASE(test_1, * boost::unit_test::timeout(20)) {
    checkStats(9);
}

BOOST_AUTO_TEST_SUITE_END()
