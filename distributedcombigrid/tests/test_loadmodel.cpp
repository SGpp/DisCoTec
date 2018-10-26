#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <chrono>
#include <thread>

#include "sgpp/distributedcombigrid/loadmodel/LearningLoadModel.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/utils/Stats.hpp"
#include "test_helper.hpp"

using namespace combigrid;

void testDataSave(int size){

    BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(size));

    MPI_Barrier(MPI_COMM_WORLD);

    CommunicatorType comm = TestHelper::getComm(size);

    if (comm == MPI_COMM_NULL){ return; }

    combigrid::Stats::initialize();

    size_t ngroup = size;
    size_t nprocs = 1;

    LearningLoadModel llmodel;

    Stats::Event e = Stats::Event();

    LevelVector l = { getCommRank(comm), 1, 4, 4, 3, 1 };
    e.end = std::chrono::high_resolution_clock::now();

    llmodel.addDataPoint(l, e, nprocs);

    BOOST_TEST(llmodel.eval(l) == (e.end - e.start).count());

    for (auto& file: llmodel.files_){
        file.second->remove_file();
    }
    combigrid::Stats::finalize();
}

void testCSVFile(){
    int rank = getCommRank(MPI_COMM_WORLD);

    csvfile file("./test" + std::to_string(rank) + ".csv");
    file << rank << std::to_string(rank * -1.0) << endrow;

    std::vector<long int> col = file.readColumn<long int>(0); 
    std::vector<long int> col_float = file.readColumn<long int>(1); 
    BOOST_TEST(col.front() == rank);    
    BOOST_TEST(col_float.front() == (rank * -1.0));
    file.remove_file();
}

BOOST_AUTO_TEST_SUITE(loadmodel)

BOOST_AUTO_TEST_CASE(test_1) {
    testDataSave(1);
}

BOOST_AUTO_TEST_CASE(test_9, * boost::unit_test::timeout(20)) {
    testDataSave(9);
}

BOOST_AUTO_TEST_CASE(test_c1, * boost::unit_test::timeout(20)) {
    testCSVFile();
}

BOOST_AUTO_TEST_SUITE_END()
