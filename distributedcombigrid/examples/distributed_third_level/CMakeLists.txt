# Copyright (C) 2008-today The SG++ Project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

cmake_minimum_required(VERSION 3.18)

project("DisCoTec distributed third level example"
        LANGUAGES CXX
        DESCRIPTION "TODO")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

if (NOT TARGET discotec)
    add_subdirectory(../../src/sgpp/distributedcombigrid discotec)
endif ()

find_package(MPI REQUIRED)


add_executable(distributed_third_level_combi_example combi_example.cpp)
target_include_directories(distributed_third_level_combi_example PRIVATE ${MPI_CXX_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/../../src)
target_link_libraries(distributed_third_level_combi_example PRIVATE MPI::MPI_CXX discotec)


add_executable(distributed_third_level_combi_example_manager_only combi_example_manager_only.cpp)
target_include_directories(distributed_third_level_combi_example_manager_only PRIVATE ${MPI_CXX_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/../../src)
target_link_libraries(distributed_third_level_combi_example_manager_only PRIVATE discotec)

#TODO: only once for all examples?
add_executable(distributed_third_level_get_mc_variance get_mc_variance.cpp)
target_include_directories(distributed_third_level_get_mc_variance PRIVATE ${MPI_CXX_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/../../src)
target_link_libraries(distributed_third_level_get_mc_variance PRIVATE discotec)

#TODO: only once for all examples?
add_executable(distributed_third_level_mpimemory mpimemory.cpp)
target_include_directories(distributed_third_level_mpimemory PRIVATE ${MPI_CXX_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/../../src)
target_link_libraries(distributed_third_level_mpimemory PRIVATE discotec)
