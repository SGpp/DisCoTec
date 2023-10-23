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

add_executable(subspace_writer subspace_writer.cpp)
target_include_directories(subspace_writer PRIVATE ${MPI_CXX_INCLUDE_DIRS} ../../../src)
target_compile_features(subspace_writer PRIVATE cxx_std_17)
target_link_libraries(subspace_writer PRIVATE MPI::MPI_CXX discotec Boost::boost)

install(TARGETS subspace_writer DESTINATION tools/subspace_writer)
