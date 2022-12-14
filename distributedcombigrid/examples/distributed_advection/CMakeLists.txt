# Copyright (C) 2008-today The SG++ Project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

cmake_minimum_required(VERSION 3.18)

project("DisCoTec distributed advection example"
        LANGUAGES CXX
        DESCRIPTION "TODO")

if (NOT TARGET discotec)
    add_subdirectory(../../src/sgpp/distributedcombigrid discotec)
endif ()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(distributed_combi_example combi_example.cpp)
target_link_libraries(distributed_combi_example discotec)

#TODO: Boost
#TODO: get_mc_variance auslagern ?