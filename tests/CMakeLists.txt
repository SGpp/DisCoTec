cmake_minimum_required(VERSION 3.24.2)

project("DisCoTec Tests"
        LANGUAGES CXX
        DESCRIPTION "Tests for DisCoTec")

if (NOT TARGET discotec)
    add_subdirectory(../src discotec)
endif ()

include(CTest)

find_package(MPI REQUIRED)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB UNIT_TESTS_SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_executable(test_distributedcombigrid_boost ${UNIT_TESTS_SRC_FILES})

find_package(Boost REQUIRED COMPONENTS unit_test_framework headers serialization)

target_link_libraries(test_distributedcombigrid_boost PUBLIC discotec Boost::boost Boost::unit_test_framework Boost::serialization MPI::MPI_CXX)
target_include_directories(test_distributedcombigrid_boost PUBLIC Boost::headers)
set_property(TARGET test_distributedcombigrid_boost PROPERTY CROSSCOMPILING_EMULATOR '${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 9')

find_package(HDF5 REQUIRED)
target_link_libraries(test_distributedcombigrid_boost PUBLIC HDF5::HDF5)

if(DISCOTEC_USE_HIGHFIVE)
    find_package(HighFive REQUIRED)
    target_link_libraries(test_distributedcombigrid_boost PUBLIC HighFive)
    target_compile_definitions(test_distributedcombigrid_boost PUBLIC HAVE_HIGHFIVE)
endif()

foreach (test_file ${UNIT_TESTS_SRC_FILES})
    get_filename_component(test_name ${test_file} NAME_WE)
    string(REPLACE "test_" "" test_name ${test_name})
    if (NOT ${test_name} STREQUAL "distributedcombigrid")
        add_test(NAME ${test_name} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 9 ${CMAKE_CURRENT_SOURCE_DIR}/test_distributedcombigrid_boost --run_test=${test_name} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    endif ()
endforeach ()

install(TARGETS test_distributedcombigrid_boost DESTINATION bin/tests)