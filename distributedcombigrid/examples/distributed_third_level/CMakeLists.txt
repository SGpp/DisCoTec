# Copyright (C) 2008-today The SG++ Project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

cmake_minimum_required(VERSION 3.18)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(distributed_third_level combi_example.cpp combi_example_manager_only.cpp get_mc_variance.cpp mpimemory.cpp)
target_link_libraries(distributed_third_level discotec)

#TODO: Boost