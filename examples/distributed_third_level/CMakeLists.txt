# Copyright (C) 2008-today The SG++ Project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

cmake_minimum_required(VERSION 3.24.2)

project("DisCoTec distributed third level example"
        LANGUAGES CXX
        DESCRIPTION "Configuration for the DisCoTec third level distributed example")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

if (NOT TARGET discotec)
    add_subdirectory(../../src discotec)
endif ()

find_package(MPI REQUIRED)

find_package(Boost REQUIRED)

add_executable(distributed_third_level_combi_example combi_example.cpp)
target_include_directories(distributed_third_level_combi_example PRIVATE ${MPI_CXX_INCLUDE_DIRS} ../../src)
target_compile_features(distributed_third_level_combi_example PRIVATE cxx_std_17)
target_link_libraries(distributed_third_level_combi_example PRIVATE MPI::MPI_CXX discotec Boost::boost)

add_executable(distributed_third_level_workers_only combi_example_worker_only.cpp)
target_include_directories(distributed_third_level_workers_only PRIVATE ${MPI_CXX_INCLUDE_DIRS} ../../src)
target_compile_features(distributed_third_level_workers_only PRIVATE cxx_std_17)
target_link_libraries(distributed_third_level_workers_only PRIVATE MPI::MPI_CXX discotec Boost::boost)
if(DISCOTEC_USE_HIGHFIVE)
    enable_language(C)
    find_package(HDF5 REQUIRED)
    find_package(HighFive REQUIRED)
    target_include_directories(discotec PRIVATE ${HighFive_INCLUDE_DIRS})
    target_link_libraries(distributed_third_level_workers_only PUBLIC HighFive)
    target_compile_definitions(distributed_third_level_workers_only PUBLIC DISCOTEC_USE_HIGHFIVE)
endif()

#TODO: only once for all examples?
add_executable(distributed_third_level_get_mc_variance get_mc_variance.cpp)
target_compile_features(distributed_third_level_get_mc_variance PRIVATE cxx_std_17)
target_include_directories(distributed_third_level_get_mc_variance PRIVATE ${MPI_CXX_INCLUDE_DIRS} ../../src)
target_link_libraries(distributed_third_level_get_mc_variance PRIVATE discotec Boost::boost)

#TODO: only once for all examples?
add_executable(distributed_third_level_mpimemory mpimemory.cpp)
target_compile_features(distributed_third_level_mpimemory PRIVATE cxx_std_17)
target_include_directories(distributed_third_level_mpimemory PRIVATE ${MPI_CXX_INCLUDE_DIRS} ../../src)
target_link_libraries(distributed_third_level_mpimemory PRIVATE discotec Boost::boost)

install(TARGETS distributed_third_level_combi_example distributed_third_level_get_mc_variance distributed_third_level_mpimemory DESTINATION examples/distributed_third_level)
install(FILES ctparam ctparam_tl_system0 ctparam_tl_system1 ctscheme_tl_system_0.json ctscheme_tl_system_1.json run.sh DESTINATION examples/distributed_third_level)
