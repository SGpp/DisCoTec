#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <chrono>
#include <thread>
#include <limits>
#include <random>

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

    BOOST_REQUIRE(true); //if things go wrong weirdly, see where things go wrong

    size_t nprocs = 64000;

    std::random_device r; 
    // Choose a random duration
    std::default_random_engine e1(r());
    std::uniform_int_distribution<long> uniform_dist(1/*, default value for upper limit is std::numeric_limits<long>::max()*/);
    int d = uniform_dist(e1);
    
    Stats::Event e = Stats::Event();
    LevelVector l = { getCommRank(comm), 1, 4, 4, 3, 1 };
    std::string fileName = getFilename(l);
    remove( fileName.c_str() );

    {
        durationsFile::DurationsWriteFile writefile(l);
        
        e.end = e.start + std::chrono::milliseconds(d);
        // std::this_thread::sleep_for(std::chrono::seconds(1));
        // e.end = std::chrono::high_resolution_clock::now();
        
        writefile.write(e, nprocs); 
        writefile.write(e, nprocs); 
    }
    {
        durationsFile::DurationsWriteFile writefile(l);
        
        // e.end = e.start + std::chrono::milliseconds(100);
        // e.end = e.start + std::chrono::milliseconds(10000000);
        std::this_thread::sleep_for(std::chrono::seconds(1));
        e.end = std::chrono::high_resolution_clock::now();
        
        writefile.write(e, nprocs); 
        writefile.write(e, nprocs); 
        // durationsFile::durationInformation di = {(long int) 1000000, nprocs};
        // writefile.write(&di);
    }
    std::chrono::milliseconds msec = std::chrono::duration_cast<std::chrono::milliseconds>(e.end - e.start);

    BOOST_TEST(getFileSize(fileName) > 0);

    {
        durationsFile::DurationsReadFile readfile(l);
        std::vector<durationInformation> result = readfile.readFromBeginning(2);
        BOOST_TEST(result.back().duration == msec.count());
    }

    if(getCommRank(comm) == 0){
        std::vector<LevelVector> lvectorvector;
        for (int i=0; i < size; ++i){
            LevelVector lv = { i, 1, 4, 4, 3, 1 };
            lvectorvector.push_back( lv );
        }
        LearningLoadModel llmodel(lvectorvector);
        llmodel.setNumberOfEntriesExpected(2);

        LevelVector l = { 0, 1, 4, 4, 3, 1 };
        BOOST_TEST(llmodel.eval(l) == msec.count());
    }
    BOOST_REQUIRE(true); 
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
