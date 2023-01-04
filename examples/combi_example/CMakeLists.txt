# Copyright (C) 2008-today The SG++ Project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

cmake_minimum_required(VERSION 3.18)

project("DisCoTec distributed combi example"
        LANGUAGES CXX
        DESCRIPTION "TODO")

if (NOT TARGET discotec)
    add_subdirectory(../../src discotec)
endif ()
find_package(Boost REQUIRED)
find_package(HighFive REQUIRED)
find_package(HDF5 REQUIRED)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_executable(combi_example combi_example.cpp)
target_compile_features(combi_example PRIVATE cxx_std_17)
target_link_libraries(combi_example discotec Boost::boost HighFive::HighFive HDF5::HDF5)
