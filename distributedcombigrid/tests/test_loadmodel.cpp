#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <chrono>
#include <thread>

#include "sgpp/distributedcombigrid/loadmodel/LearningLoadModel.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/utils/Stats.hpp"
#include "test_helper.hpp"

using namespace combigrid;

std::ifstream::pos_type getFileSize(std::string filename){
    std::ifstream in(filename.c_str(), std::ifstream::ate | std::ifstream::binary);
    return in.tellg(); 
}

void testDataSave(int size){

    BOOST_REQUIRE(TestHelper::checkNumMPIProcsAvailable(size));

    MPI_Barrier(MPI_COMM_WORLD);

    CommunicatorType comm = TestHelper::getComm(size);

    if (comm == MPI_COMM_NULL){ return; }

    combigrid::Stats::initialize();

    size_t ngroup = size;
    size_t nprocs = 64000;
    std::string fileName;
    Stats::Event e = Stats::Event();
    {
        LevelVector l = { getCommRank(comm), 1, 4, 4, 3, 1 };
        durationsFile::DurationsWriteFile writefile(l);
        fileName = getFilename(l);

        e.end = std::chrono::high_resolution_clock::now();

        writefile.write(e, nprocs);
    }

    BOOST_TEST(getFileSize(fileName) > 0);

    if(getCommRank(comm) == 0){
        std::vector<LevelVector> lvectorvector;
        for (int i=0; i < size; ++i){
            LevelVector lv = { i, 1, 4, 4, 3, 1 };
            lvectorvector.push_back( lv );
        }
        LearningLoadModel llmodel(lvectorvector);

        LevelVector l = { 0, 1, 4, 4, 3, 1 };
        BOOST_TEST(llmodel.eval(l) == (e.end - e.start).count());
    }

    combigrid::Stats::finalize();
}


BOOST_AUTO_TEST_SUITE(loadmodel)

BOOST_AUTO_TEST_CASE(test_1) {
    testDataSave(1);
}

BOOST_AUTO_TEST_CASE(test_9, * boost::unit_test::timeout(20)) {
    testDataSave(9);
}

BOOST_AUTO_TEST_SUITE_END()
