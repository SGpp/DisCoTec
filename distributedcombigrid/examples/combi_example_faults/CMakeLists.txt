# Copyright (C) 2008-today The SG++ Project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

cmake_minimum_required(VERSION 3.18)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(combi_example_faults combi_example_faults.cpp)
target_link_libraries(combi_example_faults discotec)

#TODO: Boost
# TODO: test_ft_sim_*.cpp
# TODO: MPi