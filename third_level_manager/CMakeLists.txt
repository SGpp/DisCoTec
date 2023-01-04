cmake_minimum_required(VERSION 3.24.2)

project("DisCoTec Tests"
        LANGUAGES CXX
        DESCRIPTION "DisCoTec manager for distributed combitechnic")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(thirdLevelManager Params.cpp Stats.cpp System.cpp thirdLevelManager.cpp)
add_executable(test_uftp Params.cpp Stats.cpp System.cpp test_uftp.cpp)
find_package(MPI REQUIRED)
target_link_libraries(thirdLevelManager PUBLIC MPI::MPI_CXX)
target_link_libraries(test_uftp PUBLIC MPI::MPI_CXX)

target_link_libraries(thirdLevelManager PUBLIC discotec)
target_link_libraries(test_uftp PUBLIC discotec)
find_package(Boost REQUIRED COMPONENTS serialization program_options)
target_link_libraries(thirdLevelManager PUBLIC Boost::boost)
target_link_libraries(test_uftp PUBLIC Boost::boost)